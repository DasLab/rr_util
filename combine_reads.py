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
                    prog = 'combine_reads.py',
                    description = 'Takes reads files that have been split according to ref sequence chunks and concatenates them',
                    epilog = 'Input/output are .index.csv with 1-indexed start/end of reads for each ref sequence and .reads.txt')

parser.add_argument('reads_files',nargs='*', help='names of reads.txt files to merge')
parser.add_argument('-o','--out_tag',default='',help='name of output reads.txt file')
parser.add_argument('-v','--verbose',action="store_true",help='output as each sequence is processed')

args = parser.parse_args()

df_all = []

reads_files = args.reads_files
out_tags = set()
gzipped=False

def parse_reads_file_name( reads_file ):
    cols= reads_file.split('.')
    readcol = next((i for (i,s) in enumerate(cols) if 'reads' in s), None)
    prefix = '.'.join(cols[:readcol])
    suffix = '.'.join(cols[readcol:])
    return cols,readcol,prefix,suffix

for reads_file in args.reads_files:

    if reads_file[-3:]=='.gz': gzipped=True
    else: assert( not(gzipped) )

    cols,readcol,prefix,suffix = parse_reads_file_name( reads_file )
    if readcol is None: continue

    split_tag = cols[-1]
    out_tags.add( '.'.join(cols[:(readcol-1)]) )
    if split_tag.find('_')==-1: continue
    cols = split_tag.split('_')
    if len(cols)!=2: continue
    reads_files.append( reads_file )

def get_index_file( reads_file ):
    cols,readcol,prefix,suffix = parse_reads_file_name( reads_file )
    return prefix + '.index.csv'

def get_suffix( reads_file ):
    return parse_reads_file_name(reads_file)[3]

if len(out_tags)>1:
    print( 'more than one out_tag provided!', out_tags )
    exit(0)

if len(args.out_tag)>0:
    out_tag = args.out_tag.replace('.reads.txt','')
else:
    out_tag = os.path.basename( list(out_tags)[0] )

if args.verbose: print("Collecting information from index.csv files")
for (count,reads_file) in enumerate(reads_files):
    assert( reads_file.find('reads.txt')>-1)
    index_file = get_index_file( reads_file )
    if args.verbose: print( 'Doing file %s, %d out of %d files' % (index_file,count+1,len(reads_files)) )
    df_all.append( pd.read_csv( index_file ) )

read_start = 1
if args.verbose: print('Getting ready to output',index_file)
for (df_prev,df) in zip(df_all[0:-1],df_all[1:]):
    offset = df_prev['read_end'].iloc[-1]
    for n in range(len(df)):
        df.at[n,'read_start'] += offset
        df.at[n,'read_end']   += offset

df_out = pd.concat(df_all)
index_file=out_tag+'.index.csv'
df_out.to_csv( index_file, index=False )
print('\nOutputted %d reads for %d sequences to %s' % (df_out['read_end'].iloc[-1],len(df_out),index_file) )

# Now concatenate files above.
cat_reads_file = out_tag + '.'+get_suffix(reads_files[0])
command = 'cat %s > %s' % (' '.join(reads_files), cat_reads_file)
print(command)
os.system(command)
