import string
import collections
import subprocess
import itertools
import edlib
import os
import multiprocessing

from _utils import quality_trim

def pairwise_align(idx):
    i,j = idx
    if i == j:
        return(-1,-1)
    shorter_primer,longer_primer = (primers[i],primers[j]) \
                                   if len(primers[i]) < len(primers[j]) \
                                      else (primers[j],primers[i])
    editdist_cutoff = int(0.1*len(shorter_primer))
    temp = edlib.align(shorter_primer,longer_primer,mode="HW",k=editdist_cutoff) # align in infix mode,
    if temp["editDistance"] != -1:
        return (shorter_primer,longer_primer)
    else:
        return (-1,-1)

class PrimerDataStruct(object):
    '''
    '''
    def __init__(self,*args,**kwargs):
        '''
        '''
        self.primer_file = kwargs["primer_file"]
        self.cdhit_est   = kwargs["cdhit_est"]
        self.k = kwargs["k"]
        self.ncpu = kwargs["ncpu"]
        self.seqtype = kwargs["seqtype"]

        self._primer_col = 4 if self.seqtype == "rna" else 3 # rna or dna primer file parsing
        
    def _cluster_primer_seqs(self):
        ''' Pairwise comparison of primers , store information
        for other closley matching primers for each primer
        '''
        global primers
        primers = []
        self._primer_clusters = collections.defaultdict(set)
        with open(self.primer_file,"r") as IN:
            for line in IN:
                primer = line.strip("\n\r").split("\t")[self._primer_col]
                if primer not in primers:
                    primers.append(primer)
                    
        pairwise_indices = []
        for i in range(0,len(primers)):
            for j in range(0,len(primers)):
                if i == j:
                    continue
                pairwise_indices.append((i,j))

        p = multiprocessing.Pool(self.ncpu)
        for res in p.map(pairwise_align,pairwise_indices):
            if res[0] == -1:
                continue
            shorter_primer,longer_primer = res
            self._primer_clusters[shorter_primer].add(longer_primer)
            self._primer_clusters[longer_primer].add(shorter_primer)
        p.close()
        p.join()                    
            
    def _create_primer_search_datastruct(self):
        ''' Create k-mer index on the primers
        '''
        self._primer_kmers = collections.defaultdict(set)
        self._primer_info = collections.defaultdict(list) 

        primer_id=1
        with open(self.primer_file,"r") as IN:
            for line in IN:
                contents = line.strip("\n\r").split("\t")
                primer   = contents[-1]
                if primer in self._primer_info:
                    temp = [str(primer_id),contents,[]]
                    self._primer_info[primer].append(temp) # store exact matching primers in the same value bucket
                else:
                    temp = [str(primer_id),contents,[]] # empty list for clustered primers from cd-hit
                    self._primer_info[primer] = [temp]

                    # create k-mer index
                    k = 8
                    kmers = set("".join(itertools.islice(primer,i,i+k)) for i in range(len(primer)+1-k))
                    for kmer in kmers:
                        self._primer_kmers[kmer].add(primer)

                primer_id+=1

        # add close primer sequences for each primer based on pairwise clustering
        for p1 in self._primer_clusters:
            for p2 in self._primer_clusters[p]:
                self._primer_info[p1][0][1].append(p2)

    @property
    def primer_search_datastruct(self):
        '''
        '''
        self._cluster_primer_seqs()
        self._create_primer_search_datastruct()
        return (self._primer_info,self._primer_kmers)

class Trimmer(object):
    '''
    '''
    def __init__(self,*args,**kwargs):
        ''' Class constructor
        '''        
        # user must provide
        self.is_nextseq                = kwargs["is_nextseq"]
        self.max_mismatch_rate_primer  = kwargs["max_mismatch_rate_primer"]
        self.max_mismatch_rate_overlap = kwargs["max_mismatch_rate_overlap"]
        self.synthetic_oligo_len       = kwargs["synthetic_oligo_len"]
        self.overlap_check_len         = kwargs["overlap_check_len"]
        self.primer3_R1                = kwargs["primer3_R1"]
        self.primer3_R2                = kwargs["primer3_R2"]
        self.tagname_umi               = kwargs["tagname_umi"]
        self.tagname_primer            = kwargs["tagname_primer"]
        self.tagname_primer_error      = kwargs["tagname_primer_error"]        
        self.trim_custom_seq_adapter   = kwargs["trim_custom_seq_adapter"]
        
        # user can provide, if not defaults used
        self.trim_polyA = True if "polyA_trim" in kwargs else False
        self.check_primer_side = True if "check_primer_side" in kwargs else False
        self.flip_R1_R2 = True if "flip_R1_R2" in kwargs else False

        # user can overide these defaults if needed
        self._k = 8
        self._r = 30

        # other constants
        self._revcomp_table             = bytes.maketrans(b"ACTG", b"TGAC")
        self._padding                   = 5
        self._custom_sequencing_adapter = b"AATGTACAGTATTGCGTTTTG"
        
    # kmer size
    @property
    def k(self):
        return self._k
    @k.setter
    def k(self,val):
        assert isinstance(val,int),"Please specify an integer value for k"
        self._k = val
        
    # how much of read sequence to use
    @property
    def r(self):
        return self._r        
    @r.setter
    def r(self,val):
        assert isinstance(val,int),"Please specify an integer value for k"
        self._r = val

    def custom_sequencing_adapter_check(self,r1_seq):
        ''' Check for custom sequencing adapter on r1
        :param bytes r1_seq: R1 sequence
        :rtype int
        :returns end pos of adapter, -1 if not found
        '''
        alignment = edlib.align(self._custom_sequencing_adapter,
                                r1_seq[0:len(self._custom_sequencing_adapter)+3],
                                mode="SHW",task="locations")
        if float(alignment["editDistance"])/len(self._custom_sequencing_adapter) <= 0.18:
            return alignment["locations"][-1][1]
        else:
            return -1

    def primer_trim(self,primer_datastruct,r1_seq):
        ''' Trim primer to said bases if present on read, account for certain mismatch
        :param tuple primer_datastruct: (primer_info,primer_kmer)
        :param str   r1_seq: r1 read sequence
        '''
        # initialize variables
        best_editdist       = None
        best_score          = None
        best_primer         = None
        best_primer_len     = None
        candidates          = set()
        
        # unpack argument
        primer_info, primer_kmer = primer_datastruct
        
        temp = r1_seq[0:self._r].decode("ascii")
        for oligo in set("".join(itertools.islice(temp,i,i+self._k))
                         for i in range(self._r+1-self._k)):
            if oligo in primer_kmer:
                for c in primer_kmer[oligo]:
                    candidates.add(c)
                    for similar_primer in primer_info[c][0][1]: # iterate over all similar primers to this primer
                        candidates.add(similar_primer)
        
        if len(candidates) == 0: # no hits in the index, exhaustive search over all primers
            candidates = primer_info
        num_candidates = len(candidates)
        
        for primer in candidates:
            primer_len = len(primer)
            if r1_seq[0:primer_len] == primer and num_candidates == 1: # exact match and no other primers to check
                return (primer_len,primer,0,0)
            
            alignment = edlib.align(primer,r1_seq[0:len(primer)+self._padding],mode="SHW")
            editdist  = alignment["editDistance"]
            score     = float(editdist)/primer_len

            if best_score is None or score < best_score:
                best_align      = alignment
                best_primer     = primer
                best_score      = score
                best_editdist   = editdist
                best_primer_len = primer_len
            elif score == best_score:
                if best_primer_len < primer_len:
                    best_align      = alignment
                    best_primer     = primer
                    best_score      = score
                    best_editdist   = editdist
                    best_primer_len = primer_len

        assert best_score is not None
            
        if best_score <= self.max_mismatch_rate_primer:
            return (best_align["locations"][-1][1], best_primer, best_editdist, best_score) # return 0 based position where primer ends on the read
        else:
            return (-1,None,-1,-1)
        

    '''                   Overlap checks

         Primer  
         ======--------------------------=====>        R1
               xxxxxx
                                   xxxxxx 
           <===--------------------------========      R2
                                      UMI (synthetic side)

           ===   : the synthetic part we want to trim
           xxx   : the regions checked for overlap of R1 and R2
           ---   : endogenous read sequence

    '''    
        
    def synthetic_side_check(self,r1_seq,r2_seq):
        ''' Check for overlap between R1 and R2 on the synthetic/UMI side
        :param bytes/str r1_seq: R1 read sequence
        :param bytes/str r2_seq: R2 read sequence
        :rtype tuple
        :returns (start,end) pos on R1 (0-based coordinates) 
                 (-1,-1) if no alignment within mismatch rate
        '''  
        query = self.revcomp(r2_seq[self.synthetic_oligo_len:self.synthetic_oligo_len+self.overlap_check_len])
        assert len(query) != 0,"{r1_seq}\t{r2_seq}".format(r1_seq=r1_seq.decode("ascii"),r2_seq=r2_seq.decode("ascii"))  # edlib hangs forever on empty strings     
        alignment = edlib.align(query,r1_seq,mode="HW",task="locations")

        if float(alignment["editDistance"])/self.overlap_check_len <= 0.12:
            return alignment["locations"][-1]
        else:
            return (-1,-1)
        
    def primer_side_check(self,primer_end,r1_seq,r2_seq):
        '''
        '''
        query = self.revcomp(r1_seq[primer_end+1:primer_end+1+self.overlap_check_len])
        assert len(query) != 0,"{r1_seq}\t{primer_end}\t{r2_seq}".format(r1_seq=r1_seq.decode("ascii"),primer_end=primer_end,r2_seq=r2_seq.decode("ascii"))        
        alignment = edlib.align(query,r2_seq,mode="HW",task="locations")

        
        if float(alignment["editDistance"])/self.overlap_check_len <= 0.12:
            return alignment["locations"][-1]
        else:
            return (-1,-1)        
        
    def quality_trim_(self,qual_string,seq_string,cutoff,base=33):
        return quality_trim(qual_string.decode("ascii"),seq_string.decode("ascii"),cutoff,self.is_nextseq,base)

    def revcomp(self,seq):
        ''' Reverse complement a sequence
        '''
        return seq.translate(self._revcomp_table)[::-1]
