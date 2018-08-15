import argparse
import os
import multiprocessing

import qiaseq_trimmer

def init_parser():
    '''
    '''
    global parser
    parser = argparse.ArgumentParser(description = "qiaseq trimmer: read trimming using paired end overlap")
    parser.add_argument("--r1", required = True,
                        help = "Input R1 fastq file")
    parser.add_argument("--r2", required = True,
                        help = "Input R2 fastq file")
    parser.add_argument("--out_r1", required = True,
                        help = "Output file for R1 trimmed fastq")
    parser.add_argument("--out_r2", required = True,
                        help = "Output file for R2 trimmed fastq")
    parser.add_argument("--out_metrics", required = True,
                        help = "Output file path for trimming metrics")
    parser.add_argument("--primer_file", required = True,
                        help = "Primer file with 3' coordinates" )
    parser.add_argument("--seqtype", default="dna",const="dna",
                        nargs="?",choices=["dna","rna"],
                        help="Sequencing type : dna/rna. Default : %(default)s")
    parser.add_argument("--is-nextseq", action = "store_true",
                        help = "Whether this is a NextSeq sequencing run")
    parser.add_argument("--is_duplex", action = "store_true",
                         help = "Whether this is a duplex sequencing experiment")
    parser.add_argument("--check-primer-side", action = "store_true",
                        help = "User primer side overlap coordinates for trimming R2 3' end")
    parser.add_argument("--tagname_primer", default = "pr",
                        help = "Tag name for Primer ID")
    parser.add_argument("--tagname_umi", default = "mi",
                        help = "Tag name for UMI sequence")
    parser.add_argument("--primer3_bases_R1", required = False, default = 8,
                        type = int, help = "Number of bases on 3' end of primer to keep on R1")
    parser.add_argument("--primer3_bases_R2", required = False, default = 8,
                        type = int, help = "Number of bases on 3' end of primer to keep on R2")
    parser.add_argument("--synthetic_oligo_len", required = True,
                        type = int, help = "Length of synthetic region on R2 5' end")
    parser.add_argument("--overlap_check_len", default = 25,
                        type = int, help = "Sequence length for overlap check")
    parser.add_argument("--ncpu", default = multiprocessing.cpu_count(),
                        type = int, help = "Number of CPUs to use")
    parser.add_argument("--max_mismatch_rate_overlap", default = 0.12,
                        type = float, help = "Mismatch rate to tolerate for overlap check")
    parser.add_argument("--max_mismatch_rate_primer", default = 0.12,
                        type = float, help = "Mismatch rate to tolerate for primer identification")
    parser.add_argument("--cdhit_est", help = "Path to cd-hit-est program",
                        default = os.path.join(os.path.dirname(
                            os.path.dirname(os.path.abspath(__file__))),"lib/cd-hit-est"))

def main(args):
    '''
    '''
    qiaseq_trimmer.main(args)
    
if __name__ == "__main__":
    init_parser()
    args = parser.parse_args()
    main(args)

    
    
    
    
    
    
    
    
    
    
    

    
