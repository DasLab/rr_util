#!/usr/bin/env python3

import argparse
import os
import time
import shutil
import gzip
import glob

parser = argparse.ArgumentParser(
                    prog = 'bam2reads.py',
                    description = 'Takes sorted bam file and produces reads file and index for ML modeling',
                    epilog = 'Outputs are .index.csv with 1-indexed start/end of reads for each ref sequence and .reads.txt')

parser.add_argument('-b','--bam',nargs='*', required=True, help='bam files to process')
parser.add_argument('-s','--fasta', type=str,required=True, help='FASTA of reference sequences used to align bam')
parser.add_argument('-o','--out_tag',default='',help='tag used for tag.reads.txt and tags.index.csv')
parser.add_argument('-mq','--map_quality',default=10,type=int,help=argparse.SUPPRESS )#help='minimum Bowtie2 MAPQ to consider read')
parser.add_argument('-n','--num_seq_per_chunk', default=10000, type=int, help='number of sequences')

assert(shutil.which('samtools') )

def read_fasta( fasta_file, force = False ):
    if len(fasta_file)>3 and fasta_file[-3:]=='.gz': lines = gzip.open( fasta_file, 'rt' ).readlines()
    else: lines = open( fasta_file ).readlines()
    sequences = []
    headers = []
    header = None
    sequence = ''
    for line in lines:
        if len(line) > 0 and line[0] == '>':
            if header is not None:
                headers.append(header)
                sequences.append(sequence)
            sequence = ''
            header = line[1:].strip('\n')
            continue
        sequence = sequence + line.strip('\n')
    if header is not None:
        headers.append(header)
        sequences.append(sequence)
    assert( len(sequences) == len(headers ) )
    return (sequences,headers)

args = parser.parse_args()

if len(args.bam)>1 and len(args.out_tag)>1:
    print('With multiple bam files, cannot supply --out_tag. Rerun without --out_tag and output files will have tags based on .bam files')

for bam in args.bam:
    assert(bam.find('.bam')>-1)
    assert(os.path.isfile( bam ))
    assert(os.path.isfile( bam+'.bai' ))

(fa_sequences, fa_headers) = read_fasta( args.fasta )
print( 'Read in %d sequences from %s' % (len(fa_sequences),args.fasta) )

for bam in args.bam:
    out_tag = args.out_tag
    if len(out_tag)==0: out_tag = bam.replace('.bam','')
    fid_index = open('%s.index.csv' % out_tag,'w')
    fid_reads = open('%s.reads.txt' % out_tag,'w')

    # rely on awk to handle lines quickly. MAPQ is field 5.
    command = "samtools view %s | awk '$5 >= %d {print $3,$10}'" % (bam,args.map_quality)
    fid_bam = os.popen( command )

    line = fid_bam.readline()
    read_count = 1
    num_reads = 0
    fa_idx = 0 # where we are in sequence file
    header_old = ''
    read_start = 0

    print('fasta_idx,read_start,read_end,num_reads,header,sequence',file=fid_index)
    while line:
        (header,seq) = line.split()
        print(seq,file=fid_reads)
        if header != header_old:
            #num_reads = read_end-read_start+1
            while fa_headers[fa_idx].split()[0] != header:
                print(fa_idx+1,read_start,read_count-1,read_count-read_start, fa_headers[fa_idx],fa_sequences[fa_idx],sep=',',file=fid_index)
                num_reads += read_count-read_start
                read_start = read_count
                fa_idx += 1
            header_old = header
        else:
            pass

        line = fid_bam.readline() # check the next line
        read_count += 1

    # need to closeout index file.
    while fa_idx<len(fa_headers):
        print(fa_idx+1,read_start,read_count-1,fa_headers[fa_idx],fa_sequences[fa_idx],sep=',',file=fid_index)
        num_reads += read_count-read_start
        read_start = read_count
        fa_idx += 1

    assert(num_reads==read_count)
    fid_reads.close()
    fid_index.close()
    print('\nOutputted %d reads to %s' % (read_count,fid_reads.name) )
    print('Outputted information for %d sequences with a total of  %d reads to %s'  % (len(fa_headers),num_reads,fid_index.name) )
