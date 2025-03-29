#!/usr/bin/env python3

import argparse
import os
import time
import shutil
import gzip
import glob
import pandas as pd
import tempfile


parser = argparse.ArgumentParser(
                    prog = 'merge_reads.py',
                    description = 'Takes reads files from bam2reads aligned to same sequences and merges in order.',
                    epilog = 'Input/output are .index.csv with 1-indexed start/end of reads for each ref sequence and .reads.txt')

parser.add_argument('-i','--reads_files',nargs='*', required=True, help='names of reads.txt files to merge')
parser.add_argument('-o','--out_tag',default='',help='name of output reads.txt file')
parser.add_argument('-v','--verbose',action="store_true",help='output as each sequence is processed')

args = parser.parse_args()

def parse_reads_file_name( reads_file ):
    cols= reads_file.split('.')
    readcol = next((i for (i,s) in enumerate(cols) if 'reads' in s), None)
    if readcol is None: return cols,None,None,None
    prefix = '.'.join(cols[:readcol])
    suffix = '.'.join(cols[readcol:])
    return cols,readcol,prefix,suffix

df_all = []
if True: # later will process files that have same name but across a bunch of subdirectories, like ubr_merge.py

    reads_files = []
    for read_file in args.reads_files: reads_files += sorted(glob.glob(read_file))

    out_tag = args.out_tag
    prefix = parse_reads_file_name(out_tag)[2]
    if prefix is not None: out_tag = prefix

    if args.verbose: print("Collecting information from index.csv files")
    df0 = None
    for (count,reads_file) in enumerate(reads_files):
        cols,readcol,prefix,suffix = parse_reads_file_name( reads_file )
        assert( readcol is not None )
        index_file = prefix + '.index.csv'
        if args.verbose: print( 'Doing file %s, %d out of %d files' % (index_file,count+1,len(reads_files)) )
        df_all.append( pd.read_csv( index_file, usecols=['num_reads','ref_idx'] ) )
        if count==0: df0 = pd.read_csv( index_file)
        #print(reads_file,len(df_all[-1]))

    num_ref = len(df0)
    for df in df_all: assert( len(df) == num_ref )

    tmp_dir = tempfile.mkdtemp(dir="/tmp")  # Creates a unique directory in /tmp
    if args.verbose: print(f"Temporary directory created: {tmp_dir}")
    tmp_files = []
    suffix = parse_reads_file_name( reads_files[0] )[3]
    for ref_idx in df0['ref_idx']: # cleanup old files
        tmp_file = '%s/%d.%s' % (tmp_dir,ref_idx,suffix)
        tmp_files.append(tmp_file)

    for tmp_file in tmp_files:
        if os.path.isfile(tmp_file): os.remove( tmp_file )

    # Now go through each input reads file, and output blocks that correspond to each reference sequence.
    for (count,(reads_file,df)) in enumerate(zip(reads_files, df_all)):
        fid_in = open(reads_file)
        num_reads_all=df['num_reads']
        num_reads_total=0
        if args.verbose: print( 'Doing file %s, %d out of %d files' % (reads_file,count+1,len(reads_files)) )
        for (ref_idx,num_reads,tmp_file) in zip(df['ref_idx'],df['num_reads'],tmp_files):
            fid_out = open( tmp_file, 'a')
            for n in range(num_reads): fid_out.write( fid_in.readline() )
            num_reads_total+=num_reads
        line = fid_in.readline()
        if line: print("PRINT SHOULD NOT BE A LINE HERE!",reads_file,num_reads_total,line)
        assert( not fid_in.readline() ) # hit end of file?
        fid_out.close()
        fid_in.close()

    df_out=df0
    read_start = 1
    if args.verbose: print('Getting ready to output',index_file)
    for n in range(num_ref):
        num_reads = 0
        if args.verbose and (n+1) % 1000 == 0: print( 'Figuring out total reads for ref sequence %d out of %d' % (n+1, num_ref) )
        for df in df_all: num_reads += df.iloc[n]['num_reads']
        df_out.at[n,'read_start'] = read_start
        df_out.at[n,'read_end']   = read_start + num_reads - 1
        df_out.at[n,'num_reads']  = num_reads
        read_start += num_reads
    index_file=out_tag+'.index.csv'
    df_out.to_csv( index_file, index=False )
    print('Outputted %d reads for %d sequences to %s' % (read_start-1,num_ref,index_file) )

    # Now concatenate files above.
    fid_out = open(out_tag+'.'+suffix,'w')
    for (count,(tmp_file,num_reads)) in enumerate(zip(tmp_files,df_out['num_reads'])):
        fid_in = open( tmp_file )
        if args.verbose and (count+1) % 1000 == 0: print( 'Concatenating %s, file %d out of %d' % (tmp_file, count+1, len(tmp_files)) )
        for n in range( num_reads): fid_out.write( fid_in.readline() )
        assert( not fid_in.readline() ) # hit end of file?
        fid_in.close()
    fid_out.close()

    # get rid of tmp files.
    if args.verbose: print('Deleting',tmp_dir)
    shutil.rmtree(tmp_dir)



