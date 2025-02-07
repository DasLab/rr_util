#!/usr/bin/env python3

import argparse
import os
import time
import csv
import io
import shutil

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

def update_out_files( fids, out_tag,ref_idx, chunk_size, Nref ):
    if len(fids)>0:
        for fid in fids[-1]: fid.close()
    split_tag=''
    if chunk_size>0 and Nref>chunk_size:
        chunk_start = ref_idx+1
        chunk_end   = min( ref_idx + chunk_size, len(ref_headers))
        split_tag='.%07d_%07d' % (chunk_start,chunk_end)
    fid_index = open('%s%s.index.csv' % (out_tag,split_tag),'w')
    fid_reads = open('%s%s.reads.txt' % (out_tag,split_tag),'w')
    print('ref_idx,read_start,read_end,num_reads,header,sequence',file=fid_index)
    fids.append( (fid_index, fid_reads) )

def csv_format(value):
    output = io.StringIO()
    writer = csv.writer(output, quoting=csv.QUOTE_ALL)
    writer.writerow([value])
    return output.getvalue().strip()  # Strip removes trailing newline


def output_to_index( ref_idx, read_start, read_count, ref_headers, ref_sequences, fids ):
    fid_index = fids[-1][0]
    print(ref_idx+1,read_start,read_count-1,read_count-read_start,csv_format(ref_headers[ref_idx]),ref_sequences[ref_idx],sep=',',file=fid_index)
    read_start = read_count # on to the next ref sequence, record where we are in the reads
    ref_idx += 1
    if chunk_size>0 and (ref_idx+1) % chunk_size == 1:
        update_out_files( fids, out_tag,ref_idx, chunk_size, Nref)
    return ref_idx, read_start

for bam in args.bam:
    out_tag = args.out_tag
    if len(out_tag)==0: out_tag = bam.replace('.bam','')

    # rely on awk to handle lines quickly. MAPQ is field 5.
    command = "samtools view %s | awk '$5 >= %d {print $3,$10}'" % (bam,args.map_quality)
    fid_bam = os.popen( command )

    line = fid_bam.readline()
    read_start = 1 # 1-indexed
    read_count = 1 # 1-indexed
    ref_idx = 0 # where we are in sequence file, 0 indexed
    header_old = ''
    fids = [] # outfiles
    update_out_files( fids, out_tag,ref_idx, chunk_size, Nref)
    fid_index, fid_reads = fids[-1]

    while line:
        (header,seq) = line.split()
        print(seq,file=fid_reads)
        if header != header_old:
            while ref_headers[ref_idx].split()[0] != header:
                ref_idx, read_start = output_to_index( ref_idx, read_start, read_count, ref_headers, ref_sequences, fids )
            header_old = header
            fid_index, fid_reads = fids[-1]
        line = fid_bam.readline() # check the next line
        read_count += 1

    # need to closeout index file.
    while ref_idx<Nref:
        ref_idx, read_start = output_to_index( ref_idx, read_start, read_count, ref_headers, ref_sequences, fids )

    for fid in fids[-1]: fid.close()

    if len(fids)==1:
        fid_index_tag = fids[-1][0].name
        fid_reads_tag = fids[-1][1].name
    else:
        fid_index_tag = '%d files [%s ... %s]' % ( len(fids),fids[0][0].name,fids[-1][0].name)
        fid_reads_tag = '%d files [%s ... %s]' % ( len(fids),fids[0][1].name,fids[-1][1].name)

    print('\nOutputted %d reads to %s' % (read_count,fid_reads_tag) )
    print('Outputted information for %d reference sequences with a total of %d reads to %s'  % (Nref,read_count,fid_index_tag) )
