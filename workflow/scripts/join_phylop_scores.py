#!/usr/bin/env python3

import pandas as pd
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--phylop', type=str)
parser.add_argument('--bed', type=str)
parser.add_argument('--out', type=str)

args = parser.parse_args()

# Import data
phylop_files = args.phylop.split()
bed_file = args.bed

phylop = (pd.concat([pd.read_csv(phylop_file, sep = '\t', header = None,
                                    names = ['chrom', 'start', 'stop','id','phyloP'])
                                        for phylop_file in phylop_files]))
bed = pd.read_csv(bed_file, sep = '\t', header = None, names = ['chrom', 'start', 'stop', 'id'])

# (For indels) Keep the highest phylop score per id only
phylop = phylop.groupby(['chrom', 'start', 'stop','id'], as_index = False).max()

# Join the phylop scores back to the original bed file to bring back chromosomes that were not annotated
joined = bed.merge(phylop, on = ['chrom','start','stop','id'], how = 'left')

# Save output
if joined.shape[0] == bed.shape[0]:
    joined[['chrom', 'start', 'stop','phyloP', 'id']].to_csv(args.out, sep = '\t', index = False, header = False)