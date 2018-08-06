import multiprocessing
import gzip
from functools import partial

import edlib

from fastq_parser import open_by_magic,iterate_fastq
import trimmer

def wrapper_func(buffers):
    '''
    '''
    trim_obj = trimmer.Trimmer(instrument="MiSeq",max_mismatch_rate_primer=0.12,
                               max_mismatch_rate_overlap=0.12,
                               synthetic_oligo_len = 23,
                               overlap_check_len = 25
    )
    buff_R1,buff_R2 = buffers
    R1 = buff_R1.split(b"\n")
    R2 = buff_R2.split(b"\n")
    i = 1
    j = 0
    for line in zip(R1,R2):
        if line[0] == "": # last element because of the split above
            continue
        read = None        
        if i % 4 == 1: # header
            R1_readid,R2_readid = line            
            assert R1_readid.split(b" ")[0] == R2_readid.split(b" ")[0]
        elif i % 4 == 2: # seq
            R1_seq,R2_seq = line
        elif i % 4 == 3: # placeholder
            pass
        elif i % 4 == 0: # qual
            R1_qual,R2_qual = line

            R1_primer_end,primer,editdist,score = trim_obj.primer_trim(primer_datastruct,R1_seq)
            R1_qual_3prime_end = trim_obj.quality_trim_(R1_qual,R1_seq,21,False)
            R2_qual_3prime_end = trim_obj.quality_trim_(R2_qual,R2_seq,21,False)
            syn_side_overlap_start,syn_side_overlap_end = trim_obj.synthetic_side_check(R1_seq,R2_seq)
            #print(start,end)            
            j+=1
        i+=1    

    return j
        

def init(l):
    global lock
    lock = l
    
def main():
    '''
    '''
    global primer_datastruct
    primer_datastruct = trimmer.PrimerDataStruct(k=8,primer_file="DHS-101Z.primers.txt").primer_search_datastruct
    l = multiprocessing.Lock()
    p = multiprocessing.Pool(32,initializer=init, initargs=(l,))

#    f  = gzip.open("NEB_S2_L001_R1_001.fastq.gz")
#    f2 = gzip.open("NEB_S2_L001_R2_001.fastq.gz")
    f =  open("NEB_R1.fastq","rb")
    f2 = open("NEB_R2.fastq","rb")

    N=0
    
    i = 1
    chunks = []
    for num_reads in p.imap(wrapper_func,iterate_fastq(f,f2),chunksize=32):
        N+=num_reads

    p.close()
    p.join()

    f.close()
    f2.close()
        
    
if __name__ == '__main__':
    main()
