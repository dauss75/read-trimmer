import multiprocessing
import gzip
from functools import partial

import edlib
import gc

from trimmer import PrimerDataStruct, Trimmer
from _utils import two_fastq_heads


class QiaSeqTrimmer(Trimmer):
    '''
    '''
    def __init__(self,*args,**kwargs):
        super(QiaSeqTrimmer,self).__init__(*args,**kwargs)
        
        # r1,r2 info 
        self._r1_info = None
        self._r2_info = None        

    # boolean variables requested from this class
    @property
    def is_r1_primer_trimmed(self):
        return self._is_r1_primer_trimmed
    @property
    def is_r1_syn_trimmed(self):
        return self._is_r1_syn_trimmed
    @property
    def is_r2_primer_trimmed(self):
        return self._is_r2_primer_trimmed
    @property
    def is_r1_r2_overlap(self):
        return self._is_r1_r2_overlap
    @property
    def is_too_short(self):
        return self._is_too_short
    
    # R1 info
    @property
    def r1_info(self):
        return self._r1_info
    @r1_info.setter
    def r1_info(self,info_tuple):
        self._r1_info = info_tuple
    # R2 info
    @property
    def r2_info(self):
        return self._r2_info
    @r2_info.setter
    def r2_info(self,info_tuple):
        self._r2_info = info_tuple

    def reformat_readid(self,read_id,umi,primer_id):
        ''' update read id with umi and primer info
        '''
        idx = read_id.find(b" ") if read_id.find(b" ") != -1 else len(read_id)
        primer_id = primer_id.encode("ascii")
        umi_str = b":".join([self.tagname_umi,b"Z",umi])
        primer_str = b":".join([self.tagname_primer,b"Z",primer_id])
        return b" ".join([read_id[0:idx],umi_str,primer_str])
        
    def qiaseq_trim(self,primer_datastruct):
        """ Trim QiaSeq DNA/RNA reads
        """
        # init some bools
        self._is_r1_primer_trimmed = False
        self._is_r1_syn_trimmed    = False
        self._is_r2_primer_trimmed = False
        self._is_r1_r2_overlap     = False
        self._is_too_short         = False
        self.is_qual_trim_r1       = False
        self.is_qual_trim_r2       = False
        
        primer_side_overlap = False
        syn_side_overlap    = False
        primer_id = "-1"
        
        # unpack R1,R2 info
        r1_id,r1_seq,r1_qual = self._r1_info
        r2_id,r2_seq,r2_qual = self._r2_info

        # get umi
        if not r2_seq.startswith(b"N"):
            umi = r2_seq[0:12]
            synthetic_oligo_len = self.synthetic_oligo_len
        else:
            umi = r2_seq[1:13]
            synthetic_oligo_len = self.synthetic_oligo_len + 1

        # quality trimming        
        r1_qual_end = self.quality_trim_(r1_qual,r1_seq,14)
        r2_qual_end = self.quality_trim_(r2_qual,r2_seq,14)
        temp = len(r1_seq)  - r1_qual_end
        self.r1_qual_trim_len = temp
        if temp > 0:
            self.is_qual_trim_r1 = True
        temp = len(r2_seq)  - r2_qual_end
        if temp > 0:
            self.is_qual_trim_r2 = True
        self.r2_qual_trim_len = temp
        
        # update r1,r2
        r1_seq  = r1_seq[0:r1_qual_end]
        r1_qual = r1_qual[0:r1_qual_end]
        r2_seq  = r2_seq[0:r2_qual_end]
        r2_qual = r2_qual[0:r2_qual_end]        
                    
        if len(r1_seq) < 50 or len(r2_seq) < 50: # skip reads too short after qual trimming
            self._is_too_short = True
            return
            
        # track original qual trimmed lengths
        r1_len = len(r1_seq)
        r2_len = len(r2_seq)

        # trim primer on R1
        r1_primer_end_pos,primer,editdist,score = self.primer_trim(primer_datastruct,r1_seq)

        if r1_primer_end_pos == -1:
            r1_id = self.reformat_readid(r1_id,umi,primer_id)
            r2_id = self.reformat_readid(r2_id,umi,primer_id)
            self._r1_info = (r1_id,r1_seq,r1_qual)
            self._r2_info = (r2_id,r2_seq,r2_qual)
            return
            
        primer_id = str(primer_datastruct[0][primer][0][0])

        # look for overlap on synthetic side , i.e. align endo seq from R2 on R1
        syn_side_overlap_start,syn_side_overlap_end = self.synthetic_side_check(r1_seq,r2_seq)
        if syn_side_overlap_start != -1: # synthetic side overlap was successful
            syn_side_overlap = True

        if self.check_primer_side or not syn_side_overlap:
            # look for overlap on primer side , i.e. align endo seq from R1 on R2
            primer_side_overlap_start,primer_side_overlap_end = self.primer_side_check(r1_primer_end_pos,r1_seq,r2_seq)
            if primer_side_overlap_start != -1:
                primer_side_overlap = True
                
        if syn_side_overlap:
            r1_trim_end = syn_side_overlap_end + 1
            if self.check_primer_side and primer_side_overlap : # use primer side coordinates for trimming
                r2_trim_end = primer_side_overlap_end + 1 + 8 # python will truncate to string end pos if we overflow the length
            else: # use syn side coordinates
                r2_trim_end = synthetic_oligo_len + (syn_side_overlap_end - r1_primer_end_pos) + 1 + 8
        else:
            if primer_side_overlap:
                # use primer side coordinates for trimming            
                r2_trim_end = primer_side_overlap_end + 1 + 8 # python will truncate to string end pos if we overflow the length
                r1_trim_end = r1_primer_end_pos + 1 + primer_side_overlap_end + 1 #- self.overlap_check_len
            else: # no overlap on syn or primer side
                r1_trim_end = r1_len
                r2_trim_end = r2_len
                
        # update r1,r2 qual and sequence
        r1_trim_start = r1_primer_end_pos - 8 + 1
        r2_trim_start = synthetic_oligo_len
        
        r1_seq  = r1_seq[r1_trim_start:r1_trim_end]
        r1_qual = r1_qual[r1_trim_start:r1_trim_end]
        r2_seq  = r2_seq[r2_trim_start:r2_trim_end]
        r2_qual = r2_qual[r2_trim_start:r2_trim_end]
            
        # update read ids
        r1_id = self.reformat_readid(r1_id,umi,primer_id)
        r2_id = self.reformat_readid(r2_id,umi,primer_id)
        
        # update read info tuple
        self._r1_info = (r1_id,r1_seq,r1_qual)
        self._r2_info = (r2_id,r2_seq,r2_qual)
        
        # update bools
        self._is_r1_primer_trimmed = True
        if r1_trim_end < r1_len:
            self._is_r1_syn_trimmed = True
        if r2_trim_end < r2_len:
            self._is_r2_primer_trimmed = True        
        self._is_r1_r2_overlap = primer_side_overlap or syn_side_overlap
        
        assert len(r1_seq) != 0 or len(r2_seq) != 0 or len(r1_qual) != 0 or len(r2_qual) != 0
        
def wrapper_func(buffers):
    '''
    '''
    trim_obj = QiaSeqTrimmer(is_nextseq=False,max_mismatch_rate_primer=0.12,
                             max_mismatch_rate_overlap=0.12,
                             synthetic_oligo_len = 23,
                             overlap_check_len = 25,
                             check_primer_side = True,
                             tagname_umi = b"mi",
                             tagname_primer = b"pr")
                                     
    buff_r1,buff_r2 = buffers
    r1_lines = buff_r1.split(b"\n")
    r2_lines = buff_r2.split(b"\n")
    
    # init counters
    num_too_short         = 0
    num_reads             = 0
    num_r1_primer_trimmed = 0
    num_r1_syn_trimmed    = 0
    num_r2_primer_trimmed = 0
    num_r1_r2_overlap     = 0
    
    num_qual_trim_bases_r1 = 0
    num_qual_trim_bases_r2 = 0
    num_qual_trim_r1 = 0
    num_qual_trim_r2 = 0
    
    out = []
    out_lines = []
    
    i = 1
    for line in zip(r1_lines,r2_lines):
        if line[0] == "": # last element is empty because of the split("\n") above
            continue        
        if i % 4 == 1: # header
            r1_readid,r2_readid = line
        elif i % 4 == 2: # seq
            r1_seq,r2_seq = line
        elif i % 4 == 3: # placeholder
            pass
        elif i % 4 == 0: # qual
            num_reads+=1

            r1_qual,r2_qual = line
            # have R1 and R2 ready to process now
            trim_obj.r1_info = (r1_readid,r1_seq,r1_qual)
            trim_obj.r2_info = (r2_readid,r2_seq,r2_qual)

            trim_obj.qiaseq_trim(primer_datastruct)
            
            if trim_obj.is_qual_trim_r1:
                num_qual_trim_bases_r1 += trim_obj.r1_qual_trim_len
                num_qual_trim_r1 += 1
            if trim_obj.is_qual_trim_r2:
                num_qual_trim_bases_r2 += trim_obj.r2_qual_trim_len
                num_qual_trim_r2 += 1
            
            if trim_obj.is_too_short:
                num_too_short+=1
                i+=1
                continue
            
            # retrieve trimmed sequences
            trimmed_r1_info = trim_obj.r1_info
            trimmed_r2_info = trim_obj.r2_info
            trimmed_r1_lines = b"\n".join([trimmed_r1_info[0],trimmed_r1_info[1],b"+",trimmed_r1_info[2]])
            trimmed_r2_lines = b"\n".join([trimmed_r2_info[0],trimmed_r2_info[1],b"+",trimmed_r2_info[2]])

            if trim_obj.is_r1_primer_trimmed:
                num_r1_primer_trimmed+=1
            if trim_obj.is_r1_syn_trimmed:
                num_r1_syn_trimmed+=1                
            if trim_obj.is_r2_primer_trimmed:
                num_r2_primer_trimmed+=1
            if trim_obj.is_r1_r2_overlap:
                num_r1_r2_overlap+=1

            out_lines.append((trimmed_r1_lines,trimmed_r2_lines))
            
        i+=1
        
    metrics = (num_r1_primer_trimmed,num_r1_syn_trimmed,num_r2_primer_trimmed,num_r1_r2_overlap,num_too_short,num_reads,num_qual_trim_bases_r1,num_qual_trim_bases_r2,num_qual_trim_r1,num_qual_trim_r2)

    out = [out_lines, metrics]
            
    return out

            
def iterate_fastq(f,f2,ncpu,buffer_size=4*1024**2):
    ''' Copied from cutadapt, added logic to yield
    a list of buffers equal to the number of CPUs
    '''
    buf1 = bytearray(buffer_size)
    buf2 = bytearray(buffer_size)
    
    # Read one byte to make sure we are processing FASTQ
    start1 = f.readinto(memoryview(buf1)[0:1])
    start2 = f2.readinto(memoryview(buf2)[0:1])

    if (start1 == 1 and buf1[0:1] != b'@') or (start2 == 1 and buf2[0:1] != b'@'):
        raise Exception('Paired-end data must be in FASTQ format when using multiple cores')

    to_yield = []
    nchunks = 0
    while True:
        bufend1 = f.readinto(memoryview(buf1)[start1:]) + start1
        bufend2 = f2.readinto(memoryview(buf2)[start2:]) + start2
        if start1 == bufend1 and start2 == bufend2:
            break

        end1, end2 = two_fastq_heads(buf1, buf2, bufend1, bufend2)
        assert end1 <= bufend1
        assert end2 <= bufend2

        if end1 > 0 or end2 > 0:
            nchunks+=1
            to_yield.append((memoryview(buf1)[0:end1].tobytes(), memoryview(buf2)[0:end2].tobytes()))
            
        start1 = bufend1 - end1
        assert start1 >= 0
        buf1[0:start1] = buf1[end1:bufend1]
        start2 = bufend2 - end2
        assert start2 >= 0
        buf2[0:start2] = buf2[end2:bufend2]
        
        if nchunks == ncpu:
            yield to_yield
            to_yield = []
            nchunks = 0

    if start1 > 0 or start2 > 0:
        to_yield.append((memoryview(buf1)[0:start1].tobytes(), memoryview(buf2)[0:start2].tobytes()))
    if len(to_yield) > 0:
        yield to_yield

def init(l):
    global lock
    lock = l

def main():
    '''
    '''
    global primer_datastruct
    OUT_R1 = open("NEB_R1.trimmed.fastq","wb")
    OUT_R2 = open("NEB_R2.trimmed.fastq","wb")
    OUT_metrics = open("NEB_trimming.metrics.txt","w")
    
    primer_datastruct = PrimerDataStruct(k=8,primer_file="DHS-101Z.primers.txt").primer_search_datastruct
    l = multiprocessing.Lock()

#    f  = gzip.open("NEB_S2_L001_R1_001.fastq.gz")
#    f2 = gzip.open("NEB_S2_L001_R2_001.fastq.gz")
    f =  open("NEB_R1.fastq","rb")
    f2 = open("NEB_R2.fastq","rb")

    num_r1_primer_trimmed = 0
    num_r1_syn_trimmed    = 0
    num_r2_primer_trimmed = 0
    num_r1_r2_overlap     = 0
    num_too_short         = 0
    total_reads           = 0
    
    num_qual_trim_r1_bases = 0
    num_qual_trim_r2_bases = 0
    num_qual_trim_r1 = 0
    num_qual_trim_r2 = 0
    
    nchunk = 1
    
    p = multiprocessing.Pool(16,initializer=init, initargs=(l,))
    for chunks in iterate_fastq(f,f2,16):
        res = p.map(wrapper_func,chunks)

        for r1_r2_lines,counter_tup in res:

            # unpack counters
            num_r1_primer_trimmed+=counter_tup[0]
            num_r1_syn_trimmed+=counter_tup[1]
            num_r2_primer_trimmed+=counter_tup[2]
            num_r1_r2_overlap+=counter_tup[3]
            num_too_short+=counter_tup[4]
            total_reads+=counter_tup[5]
            num_qual_trim_r1_bases+=counter_tup[6]
            num_qual_trim_r2_bases+=counter_tup[7]        
            num_qual_trim_r1+=counter_tup[8]
            num_qual_trim_r2+=counter_tup[9]


            for i in range(len(r1_r2_lines)):
                trimmed_r1_lines  =  r1_r2_lines[i][0]
                trimmed_r2_lines  =  r1_r2_lines[i][1]

                OUT_R1.write(trimmed_r1_lines)
                OUT_R1.write(b"\n")
                OUT_R2.write(trimmed_r2_lines)
                OUT_R2.write(b"\n")

            nchunk+=1
            print("Processed : {n} reads".format(n=total_reads))
            
    p.close()
    p.join()

    
    out_metrics = [
        "Total read fragments : {tot_reads}",
        "Num R1 reads primer trimmed : {num_r1_pr_trimmed}",
        "Num R2 reads primer trimmed : {num_r2_pr_trimmed}",
        "Num R1 reads synthetic side trimmed : {num_syn_trimmed}",
        "Num read fragments overlapping : {num_overlap}",
        "Num read fragments dropped too short : {too_short}",
        "Avg num bases qual trimmed R1 : {qual_trim_r1}",
        "Avg num bases qual trimmed R2 : {qual_trim_r2}",
        "Num reads qual trimmed R1 : {num_qual_trim_r1}",
        "Num reads qual trimmed R2 : {num_qual_trim_r2}",        
    ]
    out_metrics_lines = "\n".join(out_metrics).format(tot_reads = total_reads, num_r1_pr_trimmed = num_r1_primer_trimmed,
                                                      num_r2_pr_trimmed = num_r2_primer_trimmed,
                                                      num_syn_trimmed = num_r1_syn_trimmed,
                                                      num_overlap = num_r1_r2_overlap,
                                                      too_short = num_too_short,
                                                      qual_trim_r1 = float(num_qual_trim_r1_bases)/(num_qual_trim_r1),
                                                      qual_trim_r2 = float(num_qual_trim_r2_bases)/(num_qual_trim_r2),
                                                      num_qual_trim_r1 = num_qual_trim_r1,
                                                      num_qual_trim_r2 = num_qual_trim_r2)
                                                      
    OUT_metrics.write(out_metrics_lines)
    OUT_metrics.write("\n")
    f.close()
    f2.close()
    OUT_R1.close()
    OUT_R2.close()
    OUT_metrics.close()
        
    
if __name__ == '__main__':
    main()
