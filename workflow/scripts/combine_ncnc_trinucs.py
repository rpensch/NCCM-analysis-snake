#!/usr/bin/env python3

import argparse
import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--counts', type=str)
parser.add_argument('--output', type=str)

args = parser.parse_args()

# Load files
count_files = args.counts.split()

# Load counts
chr_counts = pd.concat([pd.read_csv(c ,sep = '\t', header = None, names = ['context','count']) for c in count_files])

# Drop non-ATCG letters
chr_counts = chr_counts[chr_counts['context'].str.fullmatch(r'[ATCG]{3}')].reset_index(drop=True)

counts = chr_counts.groupby('context', as_index = False)['count'].sum()

n_contexts = counts.shape[0]
if n_contexts == 64:
    counts.to_csv(args.output, sep = '\t', index = False)
else:
    raise Exception(f"Number of trinucleotide contexts is incorrect, {n_contexts} instead of 64.")