#!/usr/bin/env python3

import pandas as pd
import argparse
import warnings

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--nccp', type=str)
parser.add_argument('--ncncp', type=str)
parser.add_argument('--gene_set', type=str)
parser.add_argument('--output', type=str)

args = parser.parse_args()

# Load the results for the separate non-coding regions of each gene
nccp = pd.concat([pd.read_csv(nccp_file, 
                        sep = '\t', header = None, names = ['chrom','start','end','gene','nccp']) 
                        for nccp_file in args.nccp.split()])
nccp['region_length'] = nccp['end'] - nccp['start']

ncncp = pd.concat([pd.read_csv(ncncp_file, 
                        sep = '\t', header = None, names = ['chrom','start','end','gene','ncncp']) 
                        for ncncp_file in args.ncncp.split()])
ncncp['region_length'] = ncncp['end'] - ncncp['start']

# Merge constraint and non-constraint
ncp = nccp.merge(ncncp, on = ['gene','chrom','region_length','start','end'])

# Sanity check that all the regions are the same
if (ncp.shape[0] != nccp.shape[0]) | (ncp.shape[0] != ncncp.shape[0]):
    raise Exception("Not the same regions present.")

# Calculate the overall counts for total non-coding regions of each gene
gene_ncp = (ncp.groupby('gene', as_index = False)
               .agg({'region_length':'sum','nccp':'sum','ncncp':'sum'}))

# Determine how many non-coding positions had a phylop score
gene_ncp['phylop_cov'] = (gene_ncp['nccp'] + gene_ncp['ncncp']) / gene_ncp['region_length']

# merge back to the flank bed
flanks = pd.read_csv(args.gene_set, sep = '\t', header = None, names = ['chrom','lflank','uflank','gene'])
merged = flanks.merge(gene_ncp, on = 'gene')

# Final sanity check
if (merged.shape[0] != flanks.shape[0]) | (merged.shape[0] != gene_ncp.shape[0]):
    warnings.warn("MERGE WARNING: Output does not have the expected number of genes.\n\
                    1. This could mean that some regions are fully coding and have been removed.\n\
                    2. It could also mean some sort of other problem\n\
                    Please make sure everything is okay.")

merged.to_csv(args.output, sep = '\t', index = False)