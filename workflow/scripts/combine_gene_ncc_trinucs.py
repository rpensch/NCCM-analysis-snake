#!/usr/bin/env python3

import argparse
import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--counts', type=str)
parser.add_argument('--gene_flanks', type=str)
parser.add_argument('--output', type=str)

args = parser.parse_args()

# Load files
count_files = args.counts.split()
chr_counts = pd.concat([pd.read_csv(c, sep = '\t', header = None, 
                names = ['gene','context','count']) for c in count_files])

gene_flanks = pd.read_csv(args.gene_flanks, sep = '\t', header = None, 
                names = ['chrom','start','end','gene'])

# Drop non-ATCG letters
chr_counts = chr_counts[chr_counts['context'].str.fullmatch(r'[ATCG]{3}')].reset_index(drop=True)

# Add genes with no ncc positions back in
all_genes = gene_flanks['gene'].unique().tolist()
context_genes = chr_counts['gene'].unique().tolist()
contexts_list = chr_counts['context'].unique().tolist()

add_genes = []
for gene in all_genes:
    if gene not in context_genes:
        for context in contexts_list:
            add_genes.append({'gene':gene,
                              'context':context,
                              'count':0})

chr_counts = pd.concat([chr_counts,pd.DataFrame(add_genes)])

# Save
chr_counts.to_csv(args.output, sep = '\t', index = False)