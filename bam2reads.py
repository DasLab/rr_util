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
parser.add_argument('-n','--chunk_size', default=0, type=int, help='split with this number of sequences per chunk')

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

(ref_sequences, ref_headers) = read_fasta( args.fasta )
Nref = len(ref_sequences)
print( 'Read in %d reference sequences from %s' % (Nref,args.fasta) )

chunk_size = args.chunk_size

def update_out_files( out_tag,ref_idx, chunk_size, Nref, fid_index=None, fid_reads=None):
    if fid_index is not None: fid_index.close()
    if fid_reads is not None: fid_reads.close()
    split_tag=''
    if chunk_size>0 and Nref>chunk_size:
        chunk_start = ref_idx+1
        chunk_end   = min( ref_idx + chunk_size, len(ref_headers))
        split_tag='.%07d_%07d' % (chunk_start,chunk_end)
    fid_index = open('%s%s.index.csv' % (out_tag,split_tag),'w')
    fid_reads = open('%s%s.reads.txt' % (out_tag,split_tag),'w')
    return fid_index, fid_reads

for bam in args.bam:
    out_tag = args.out_tag
    if len(out_tag)==0: out_tag = bam.replace('.bam','')

    # rely on awk to handle lines quickly. MAPQ is field 5.
    command = "samtools view %s | awk '$5 >= %d {print $3,$10}'" % (bam,args.map_quality)
    fid_bam = os.popen( command )

    line = fid_bam.readline()
    read_count = 1
    num_reads = 0
    ref_idx = 0 # where we are in sequence file
    header_old = ''
    read_start = 0
    fid_index, fid_reads = update_out_files( out_tag,ref_idx, chunk_size, Nref)

    print('fasta_idx,read_start,read_end,num_reads,header,sequence',file=fid_index)
    while line:
        (header,seq) = line.split()
        print(seq,file=fid_reads)
        if header != header_old:
            #num_reads = read_end-read_start+1
            while ref_headers[ref_idx].split()[0] != header:
                print(ref_idx+1,read_start,read_count-1,read_count-read_start, ref_headers[ref_idx],ref_sequences[ref_idx],sep=',',file=fid_index)
                num_reads += read_count-read_start
                read_start = read_count
                ref_idx += 1
                if chunk_size>0 and (ref_idx+1) % chunk_size == 1: # about to start the next chunk.
                    fid_index, fid_reads = update_out_files( out_tag,ref_idx, chunk_size, Nref, fid_index, fid_reads)

            header_old = header

        line = fid_bam.readline() # check the next line
        read_count += 1

    # need to closeout index file.
    while ref_idx<Nref:
        print(ref_idx+1,read_start,read_count-1,ref_headers[ref_idx],ref_sequences[ref_idx],sep=',',file=fid_index)
        num_reads += read_count-read_start
        read_start = read_count
        ref_idx += 1

    assert(num_reads==read_count)
    fid_reads.close()
    fid_index.close()
    print('\nOutputted %d reads to %s' % (read_count,fid_reads.name) )
    print('Outputted information for %d sequences with a total of %d reads to %s'  % (Nref,num_reads,fid_index.name) )
