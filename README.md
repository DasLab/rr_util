# rr_util
Utilities for generating raw reads datasets for deep neural network training


```
usage: bam2reads.py [-h] -b BAM -s FASTA [-o OUT_TAG] [-n NUM_SEQ_PER_CHUNK]

Takes sorted bam file and produces reads file and index for ML modeling

options:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     bam files to process
  -s FASTA, --fasta FASTA
                        FASTA of reference sequences used to align bam
  -o OUT_TAG, --out_tag OUT_TAG
                        tag used for tag.reads.txt and tags.index.csv
  -n NUM_SEQ_PER_CHUNK, --num_seq_per_chunk NUM_SEQ_PER_CHUNK
                        number of sequences

Outputs are .index.csv with 1-indexed start/end of reads for each ref sequence and .reads.txt
```
