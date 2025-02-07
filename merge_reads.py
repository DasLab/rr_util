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

parser.add_argument('-i','--reads_files',nargs='*', required=True, help=argparse.SUPPRESS) # devel option
parser.add_argument('-o','--out_tag',default='',help=argparse.SUPPRESS) # devel option

args = parser.parse_args()

df_all = []
if True: # later will process files that have same name but across a bunch of subdirectories, like ubr_merge.py

    reads_files = args.reads_files
    if len(reads_files)==1:
        reads_files=sorted(glob.glob(reads_files[0]))
    out_tag = args.out_tag.replace('.reads.txt','')

    for reads_file in reads_files:
        assert( reads_file.find('.reads.txt')>-1)
        index_file = reads_file.replace('.reads.txt','.index.csv')
        df_all.append( pd.read_csv( index_file ) )
    num_ref = len(df_all[0])
    for df in df_all: assert( len(df) == num_ref )

    tmp_dir = tempfile.mkdtemp(dir="/tmp")  # Creates a unique directory in /tmp
    print(f"Temporary directory created: {tmp_dir}")
    tmp_files = []
    for ref_idx in df_all[0]['ref_idx']: # cleanup old files
        tmp_file = '%s/%d.reads.txt' % (tmp_dir,ref_idx)
        tmp_files.append(tmp_file)

    for tmp_file in tmp_files:
        if os.path.isfile(tmp_file): os.remove( tmp_file )

    # Now go through each input reads file, and output blocks that correspond to each reference sequence.
    for (reads_file,df) in zip(reads_files, df_all):
        fid_in = open(reads_file)
        num_reads_all=df['num_reads']
        num_reads_total=0
        for (ref_idx,num_reads,tmp_file) in zip(df['ref_idx'],df['num_reads'],tmp_files):
            fid_out = open( tmp_file, 'a')
            for n in range(num_reads): fid_out.write( fid_in.readline() )
            num_reads_total+=num_reads
        line = fid_in.readline()
        if line: print("PRINT SHOULD NOT BE A LINE HERE!",reads_file,num_reads_total,line)
        assert( not fid_in.readline() ) # hit end of file?
        fid_out.close()
        fid_in.close()

    df_out=df_all[0]
    read_start = 1
    for n in range(num_ref):
        num_reads = 0
        for df in df_all: num_reads += df.iloc[n]['num_reads']
        df_out.at[n,'read_start'] = read_start
        df_out.at[n,'read_end']   = read_start + num_reads - 1
        df_out.at[n,'num_reads']  = num_reads
        read_start += num_reads
    index_file=out_tag+'.index.csv'
    df.to_csv( index_file, index=False )
    print('Outputted %d reads for %d sequences to %s' % (read_start-1,num_ref,index_file) )


    # Now concatenate files above.
    fid_out = open(out_tag+'.reads.txt','w')
    for (tmp_file,num_reads) in zip(tmp_files,df_out['num_reads']):
        fid_in = open( tmp_file )
        for n in range( num_reads): fid_out.write( fid_in.readline() )
        assert( not fid_in.readline() ) # hit end of file?
        fid_in.close()
    fid_out.close()

    # get rid of tmp files.
    shutil.rmtree(tmp_dir)



