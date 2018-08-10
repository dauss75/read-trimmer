import string
import collections
import subprocess
import itertools
import edlib
import os

from _utils import quality_trim

class PrimerDataStruct(object):
    '''
    '''
    def __init__(self,*args,**kwargs):
        '''
        '''
        self.primer_file = kwargs["primer_file"]
        self.cdhit_est   = kwargs["cdhit_est"]
        self.k = kwargs["k"]

        self._cdhit_fasta      =   self.primer_file + ".fasta"
        self._cdhit_out_prefix =   self.primer_file + ".clusters.temp"
        self._cdhit_out        =   self.primer_file + ".clusters.temp.clstr"
        self._cdhit_out_parsed =   self.primer_file + ".clusters"

    def _parse_cdhit(self):
        ''' Parse cd-hit output to usable format
        '''
        primers = []
        i=0

        with open(self.primer_file,"r") as IN:
            for line in IN:
                primers.append(line.strip("\n\r").split("\t")[-1])

        cluster_info = collections.defaultdict(list)
        with open(self._cdhit_out,"r") as IN:
            for line in IN:
                if line.startswith(">"):
                    cluster_num = int(line.strip("\n").split(" ")[1])
                    new_cluster = True
                else:
                    temp = line.strip("\n").split("\t")[1].split(",")[1].strip().split("...")[0].strip(">")
                    primer = primers[int(temp)]
                    cluster_info[str(cluster_num)].append(primer)

        with open(self._cdhit_out_parsed,"w") as OUT:
            for cluster in cluster_info:
                if len(cluster_info[cluster]) > 1:
                    OUT.write(",".join(cluster_info[cluster])+"\n")
        

    def _cluster_primer_seqs(self):
        ''' Cluster similar primer sequences for use later in primer trimming
        '''
        # create a fasta of the primers
        cmd1 = (
            """ i=0; while read chrom pos strand primer; do echo ">"$i; echo $primer|tr -d '\r'; """
            """ i=$(($i+1)); done < {primerfile} > {fasta}; """.format(
                primerfile = self.primer_file,fasta=self._cdhit_fasta)
        )
        subprocess.check_call(cmd1,shell=True)
        # run cd-hit
        cmd2 = "{cdhit_est} -i {fasta} -o {out_prefix}".format(
            cdhit_est = self.cdhit_est,
            fasta = self._cdhit_fasta,
            out_prefix = self._cdhit_out_prefix)
        subprocess.check_call(cmd2,shell=True)
        # parse output to usable format
        self._parse_cdhit()

    def _cleanup(self):
        ''' Clean-up temp files
        '''
        for f in [self._cdhit_fasta,self._cdhit_out_prefix,self._cdhit_out]:
            os.remove(f)
            
    def _create_primer_search_datastruct(self):
        '''
        '''
        self._primer_kmers = collections.defaultdict(set)
        self._primer_info = collections.defaultdict(list) 

        primer_id=1
        with open(self.primer_file,"r") as IN:
            for line in IN:
                contents = line.strip("\n\r").split("\t")
                primer   = contents[-1]
                if primer in self._primer_info:
                    temp = [str(primer_id),[]]
                    self._primer_info[primer].append(temp) # store exact matching primers in the same value bucket
                else:
                    temp = [str(primer_id),[]] # empty list for clustered primers from cd-hit
                    self._primer_info[primer] = [temp]

                    # create k-mer index
                    k = 8
                    kmers = set("".join(itertools.islice(primer,i,i+k)) for i in range(len(primer)+1-k))
                    for kmer in kmers:
                        self._primer_kmers[kmer].add(primer)

                primer_id+=1

        # add close primer sequences for each primer based on cd-hit results
        with open(self._cdhit_out_parsed,"r") as IN:
            for line in IN:
                temp = line.strip("\n").split(",")
                
                for p1,p2 in itertools.product(temp,repeat=2):                    
                    assert p1 in self._primer_info and p2 in self._primer_info,"Primer(s) from cd-hit not in PrimerFile !!"
                    if p1 == p2: # same primer sequences
                        continue

                    self._primer_info[p1][0][1].append(p2)
                    self._primer_info[p2][0][1].append(p1)                

    @property
    def primer_search_datastruct(self):
        '''
        '''
        self._cluster_primer_seqs()
        self._create_primer_search_datastruct()
        self._cleanup()
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
                                r1_seq[0:len(self._custom_sequencing_adapter)+4],
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
