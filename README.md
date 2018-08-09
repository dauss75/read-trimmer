# read-trimmer
Paired end trimmer for QIASeq DNA/RNA reads 

## Installation
```
pip3 install edlib
pip3 install Cython
git clone https://github.com/qiaseq/read-trimmer.git
cd read-trimmer
python3 setup.py build_ext --inplace
# for cd-hit-est
sudo apt-get install libstdcc++6 
```

## Usage Instructions
```
Run as python3 read-trimmer/trimmer/run.py

usage: run.py [-h] --r1 R1 --r2 R2 --out_r1 OUT_R1 --out_r2 OUT_R2
              --out_metrics OUT_METRICS --primer_file PRIMER_FILE
              [--is-nextseq] [--check-primer-side]
              [--tagname_primer TAGNAME_PRIMER] [--tagname_umi TAGNAME_UMI]
              --synthetic_oligo_len SYNTHETIC_OLIGO_LEN
              [--overlap_check_len OVERLAP_CHECK_LEN] [--ncpu NCPU]
              [--max_mismatch_rate_overlap MAX_MISMATCH_RATE_OVERLAP]
              [--max_mismatch_rate_primer MAX_MISMATCH_RATE_PRIMER]
              [--cdhit_est CDHIT_EST]

qiaseq trimmer: read trimming using paired end overlap

optional arguments:
  -h, --help            show this help message and exit
  --r1 R1               Input R1 fastq file
  --r2 R2               Input R2 fastq file
  --out_r1 OUT_R1       Output file for R1 trimmed fastq
  --out_r2 OUT_R2       Output file for R2 trimmed fastq
  --out_metrics OUT_METRICS
                        Output file path for trimming metrics
  --primer_file PRIMER_FILE
                        Primer file with 3' coordinates
  --is-nextseq          Whether this is a NextSeq sequencing run
  --check-primer-side   User primer side overlap coordinates for trimming R2
                        3' end
  --tagname_primer TAGNAME_PRIMER
                        Tag name for Primer ID
  --tagname_umi TAGNAME_UMI
                        Tag name for UMI sequence
  --synthetic_oligo_len SYNTHETIC_OLIGO_LEN
                        Length of synthetic region on R2 5' end
  --overlap_check_len OVERLAP_CHECK_LEN
                        Sequence length for overlap check
  --ncpu NCPU           Number of CPUs to use
  --max_mismatch_rate_overlap MAX_MISMATCH_RATE_OVERLAP
                        Mismatch rate to tolerate for overlap check
  --max_mismatch_rate_primer MAX_MISMATCH_RATE_PRIMER
                        Mismatch rate to tolerate for primer identification
  --cdhit_est CDHIT_EST
                        Path to cd-hit-est program


```

