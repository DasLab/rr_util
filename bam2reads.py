#!/usr/bin/env python3

import argparse
import os
import time
import csv
import io
import shutil
from read_util import *

parser = argparse.ArgumentParser(
                    prog = 'bam2reads.py',
                    description = 'Takes sorted bam file and produces reads file and index for ML modeling',
                    epilog = 'Outputs are .index.csv with 1-indexed start/end of reads for each ref sequence and .reads.txt')

parser.add_argument('-b','--bam','--cram',nargs='*', required=True, help='bam files to process')
parser.add_argument('-s','--fasta', type=str,required=True, help='FASTA of reference sequences used to align bam')
parser.add_argument('-o','--out_tag',default='',help='tag used for tag.reads.txt and tags.index.csv')
parser.add_argument('-mq','--map_quality',default=10,type=int,help='minimum Bowtie2 MAPQ to consider read')
parser.add_argument('--mutdel_cutoff',type=int,default=10,help='Filter for maximum number of mut/del in read (default 0 means no filter)' )
parser.add_argument('-n','--chunk_size', default=0, type=int, help='split with this number of sequences per chunk')
parser.add_argument('--start_idx', default=0, type=int, help='only do the reference sequences from start_idx onwards [default all]')
parser.add_argument('--end_idx', default=0, type=int, help='only do the reference sequences up to end_idx [default all]')
parser.add_argument('--check_md',action='store_true',help='Run check on md by inferring align_read' )

assert(shutil.which('samtools') )

def read_fasta( fasta_file ):
    """Reads a FASTA file and returns a list of sequences and headers.

    Args:
        fasta_file (str): The path to the FASTA file.  Supports gzipped files (`.gz`).

    Returns:
        tuple: A tuple containing two lists:
            - sequences: A list of strings, where each string is a FASTA sequence.
            - headers: A list of strings, where each string is a FASTA header (the line starting with '>').

    Raises:
        AssertionError: If the number of sequences and headers do not match.

    Notes:
        Handles both standard FASTA files and those compressed with gzip.
    """
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


def update_out_files( fids, out_tag,ref_idx, chunk_size, Nref, start_idx, end_idx ):
    """Updates output files based on chunk size and reference index.

    This function closes previously opened files if any exist, and opens new files for writing
    the next chunk of data.  It handles splitting the output into multiple files based on `chunk_size`.

    Args:
        fids (list): A list of tuples, where each tuple contains file handles for index and reads files.
        out_tag (str): The base name for output files.
        ref_idx (int): The current index of the reference sequence.
        chunk_size (int): The size of each chunk.  If 0, no chunking is performed.
        Nref (int): The total number of reference sequences.

    Returns:
        None. Modifies `fids` in place by appending new file handles.
    """
    if len(fids)>0:
        for fid in fids[-1]: fid.close()
    split_tag=''
    if chunk_size>0 and Nref>chunk_size:
        chunk_start = ref_idx+1
        chunk_end   = min( ref_idx + chunk_size, len(ref_headers))
        split_tag='.%07d_%07d' % (chunk_start,chunk_end)
    if start_idx>0:
        assert( chunk_size == 0 )
        chunk_start = start_idx
        chunk_end = Nref
        if end_idx > 0: chunk_end = end_idx
        split_tag='.%07d_%07d' % (chunk_start,chunk_end)
    if ref_idx < Nref-1:
        fid_index = open('%s%s.index.csv' % (out_tag,split_tag),'w')
        fid_raw_reads = open('%s%s.raw_reads.txt' % (out_tag,split_tag),'w')
        fid_align_reads = open('%s%s.align_reads.txt' % (out_tag,split_tag),'w')
        fid_start_md = open('%s%s.start_md_reads.txt' % (out_tag,split_tag),'w')
        print('ref_idx,read_start,read_end,num_reads,header,sequence',file=fid_index)
        fids.append( (fid_index, fid_raw_reads, fid_align_reads, fid_start_md ) )

def csv_format(value):
    """Formats a value as a CSV string, quoting all fields.

    Args:
        value: The value to format.

    Returns:
        str: A CSV-formatted string representation of the value.
    """
    output = io.StringIO()
    writer = csv.writer(output, quoting=csv.QUOTE_ALL)
    writer.writerow([value])
    return output.getvalue().strip()  # Strip removes trailing newline

def output_to_index( ref_idx, read_start, read_count, ref_headers, ref_sequences, fids, total_reads, start_idx, end_idx ):
    """Writes data to the index file and handles chunk updates.

    Args:
        ref_idx (int): The index of the current reference sequence.
        read_start (int): The starting index of reads.
        read_count (int): The total count of reads.
        ref_headers (list): A list of reference headers.
        ref_sequences (list): A list of reference sequences.
        fids (list): A list of file handles.
        total_reads (int): The running total of reads processed.

    Returns:
        tuple: A tuple containing updated ref_idx, read_start, read_count, and total_reads.

    """
    fid_index = fids[-1][0]
    ok_start_end_idx = start_idx==0 or (ref_idx+1 >= start_idx and (end_idx==0 or ref_idx+1 <=end_idx) )
    if ok_start_end_idx:
        print(ref_idx+1,read_start,read_count,read_count-read_start+1,csv_format(ref_headers[ref_idx]),ref_sequences[ref_idx],sep=',',file=fid_index)
    total_reads += read_count-read_start+1
    ref_idx += 1
    if chunk_size>0 and (ref_idx+1) % chunk_size == 1:
        update_out_files( fids, out_tag,ref_idx, chunk_size, Nref, start_idx, end_idx)
        read_count = 0 # reset
    read_start = read_count+1 # on to the next ref sequence, record where we are in the reads
    return ref_idx, read_start, read_count, total_reads

args = parser.parse_args()

if len(args.bam)>1 and len(args.out_tag)>1:
    print('With multiple bam files, cannot supply --out_tag. Rerun without --out_tag and output files will have tags based on .bam files')

for bam in args.bam:
    assert(bam.find('.bam')>-1 or bam.find('.cram')>-1)
    assert(os.path.isfile( bam ))
    #assert(os.path.isfile( bam+'.bai' ))

(ref_sequences, ref_headers) = read_fasta( args.fasta )
ref_sequences = [sequence.replace('U','T') for sequence in ref_sequences] # need to convert to DNA!
Nref = len(ref_sequences)
print( 'Read in %d reference sequences from %s' % (Nref,args.fasta) )


chunk_size = args.chunk_size
start_idx = args.start_idx
end_idx = args.end_idx
assert( not( chunk_size>0 and start_idx>0) )
md_idx = -1

for bam in args.bam:
    out_tag = args.out_tag
    if len(out_tag)==0: out_tag = bam.replace('.bam','').replace('.cram','')

    # rely on awk to handle lines quickly. MAPQ is field 5.
    command = "samtools view %s | awk '$5 >= %d'" % (bam,args.map_quality)
    fid_bam = os.popen( command )

    line = fid_bam.readline()
    cols = line.split()
    read_count = 0 # 1-indexed
    read_start = 1 # 1-indexed
    total_reads = 0
    ref_idx = 0 # where we are in reference sequence file, 0 indexed
    header = cols[2]

    fids = [] # outfiles
    update_out_files( fids, out_tag,ref_idx, chunk_size, Nref, start_idx, end_idx)
    fid_index, fid_raw_reads, fid_align_reads, fid_start_md = fids[-1]

    while line:
        cols = line.split()
        next_header = cols[2]
        if header != next_header: # switching to a new set of reads
            while ref_headers[ref_idx].split()[0] != next_header: #print index information for last set of reads until we reach next_header
                ref_idx, read_start, read_count, total_reads = output_to_index( ref_idx, read_start, read_count, ref_headers, ref_sequences, fids, total_reads, start_idx, end_idx )
            header = next_header
            fid_index, fid_raw_reads, fid_align_reads, fid_start_md = fids[-1]

        ok_start = start_idx == 0 or ref_idx+1 >= start_idx
        ok_end   = end_idx == 0 or ref_idx+1 <= end_idx

        if ok_start and ok_end:
            ref_seq = ref_sequences[ref_idx]
            total_mutdel, align_read = get_total_mutdel( cols,ref_seq )
            if args.mutdel_cutoff == 0 or total_mutdel <= args.mutdel_cutoff:
                start_md = get_md_convert( cols,ref_seq )
                if args.check_md: check_md_convert( start_md, ref_seq, align_read )
                read = cols[9]
                print(read, file=fid_raw_reads)
                print(align_read, file=fid_align_reads)
                print(','.join(start_md), file=fid_start_md)
                read_count += 1

        line = fid_bam.readline() # check the next line

    # need to close out index file.
    while ref_idx<Nref:
        ref_idx, read_start, read_count, total_reads = output_to_index( ref_idx, read_start, read_count, ref_headers, ref_sequences, fids, total_reads, start_idx, end_idx )

    for fid in fids[-1]: fid.close()

    if len(fids)==1:
        fid_index_tag = fids[-1][0].name
        fid_raw_reads_tag = fids[-1][1].name
        fid_align_reads_tag = fids[-1][2].name
        fid_start_md_tag = fids[-1][3].name
    else:
        fid_index_tag = '%d files [%s ... %s]' % ( len(fids),fids[0][0].name,fids[-1][0].name)
        fid_raw_reads_tag = '%d files [%s ... %s]' % ( len(fids),fids[0][1].name,fids[-1][1].name)
        fid_align_reads_tag = '%d files [%s ... %s]' % ( len(fids),fids[0][2].name,fids[-1][2].name)
        fid_start_md_tag = '%d files [%s ... %s]' % ( len(fids),fids[0][3].name,fids[-1][3].name)

    print('\nOutputted %d raw reads to %s' % (total_reads,fid_raw_reads_tag) )
    print('Outputted %d aligned reads to %s' % (total_reads,fid_align_reads_tag) )
    print('Outputted %d start/md to %s' % (total_reads,fid_start_md_tag) )
    print('Outputted information for %d reference sequences with a total of %d reads to %s'  % (Nref,total_reads,fid_index_tag) )
