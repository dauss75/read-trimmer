import sys
import os
from argparse import Namespace


sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"trimmer/"))
from qiaseq_trimmer import QiaSeqTrimmer

def helper_return_qiaseq_obj(args):
    ''' Helper function to initiliaze Trimmer Object
    '''
    return QiaSeqTrimmer(is_nextseq                =   args.is_nextseq,
                         is_duplex                 =   args.is_duplex,
                         seqtype                   =   args.seqtype,
                         max_mismatch_rate_primer  =   args.max_mismatch_rate_primer,
                         max_mismatch_rate_overlap =   args.max_mismatch_rate_overlap,
                         custom_seq_adapter        =   args.custom_seq_adapter,
                         umi_len                   =   args.umi_len,
                         common_seq_len            =   args.common_seq_len,
                         overlap_check_len         =   args.overlap_check_len,
                         min_primer_side_len       =   args.min_primer_side_len,
                         min_umi_side_len          =   args.min_umi_side_len,
                         check_primer_side         =   args.check_primer_side,
                         umi_filter_min_bq         =   args.umi_filter_min_bq,
                         umi_filter_max_lowQ_bases =   args.umi_filter_max_lowQ_bases,
                         umi_filter_max_Ns         =   args.umi_filter_max_Ns,
                         trim_custom_seq_adapter   =   args.trim_custom_seq_adapter,
                         primer3_R1                =   args.primer3_R1,
                         primer3_R2                =   args.primer3_R2,
                         poly_tail_primer_side     =   args.poly_tail_primer_side,
                         poly_tail_umi_side        =   args.poly_tail_umi_side,
                         tagname_umi               =   args.tagname_umi,
                         tagname_duplex            =   args.tagname_duplex,
                         tagname_primer            =   args.tagname_primer,
                         tagname_primer_error      =   args.tagname_primer_error,
                         tag_seperator             =   args.tag_seperator,
                         no_tagnames               =   args.no_tagnames)


def helper_return_args():
    ''' Helper function to initiliaze arguments
    '''
    return Namespace(is_nextseq                =   False,
                     is_duplex                 =   False,
                     seqtype                   =   "dna",
                     max_mismatch_rate_primer  =   0.12,
                     max_mismatch_rate_overlap =   0.12,
                     custom_seq_adapter        =   b"AATGTACAGTATTGCGTTTTG",
                     umi_len                   =   12,
                     common_seq_len            =   11,
                     overlap_check_len         =   25,
                     min_primer_side_len       =   50,
                     min_umi_side_len          =   50,
                     check_primer_side         =   True,
                     umi_filter_min_bq         =   20,
                     umi_filter_max_lowQ_bases =   2,
                     umi_filter_max_Ns         =   2,
                     trim_custom_seq_adapter   =   False,
                     primer3_R1                =   8,
                     primer3_R2                =   8,
                     poly_tail_primer_side     =   "none",
                     poly_tail_umi_side        =   "none",
                     tagname_umi               =   b"mi",
                     tagname_duplex            =   b"DU",
                     tagname_primer            =   b"pr",
                     tagname_primer_error      =   b"pe",
                     tag_seperator             =   b"\t",
                     no_tagnames               =   False)
        
