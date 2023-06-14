#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip

from vcf import load_vcf, get_tumour_sample
from annotation import load_bed, extract_annotations, extract_vaf, get_coding_type

# Define composite matrix columns to keep
basic_columns = ['chrom', 'pos', 'ref', 'alt'] + ['sample'] + ['vaf'] + ['coding_type']
phylop_columns = ['start', 'stop', 'phylop', 'mutation_id']
annotations = "Allele | Annotation | Annotation_Impact | Gene_Name | "\
              "Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | "\
              "Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | "\
              "AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO".lower().split(' | ')

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sample', type=str)
parser.add_argument('--ann', type=str)
parser.add_argument('--phylop', type=str)
parser.add_argument('-o', '--out', type=str)

args = parser.parse_args()

# Load files
ann = load_vcf(args.ann)
ann.columns = [c.lower().replace('#','') for c in ann.columns.tolist()]
phylop = load_bed(args.phylop)
phylop.columns = [c.lower().replace('#','') for c in phylop.columns.tolist()]

# Make mutation id for annotation data
ann['mutation_id'] = (args.sample + ':' + ann['chrom'] + ':' + 
                     ann['pos'].astype(str) + ':' + ann['ref'] + ':' + ann['alt'])

# Remove duplicate mutations
ann = ann.drop_duplicates('mutation_id')
phylop = phylop.drop_duplicates('mutation_id')

# Merge annotations and phylop data based on mutation id
composite = ann.merge(phylop[phylop_columns], on = 'mutation_id', how = 'inner')

# Add a sample columns
composite['sample'] = args.sample

# If there are no mutations (the vcf was empty), save empty and exit
if (ann.shape[0] == 0) & (phylop.shape[0] == 0) & (composite.shape[0] == 0):
    composite[basic_columns + phylop_columns + annotations] = 0
    (composite[basic_columns + phylop_columns + annotations]
        .to_csv(args.out, index = False, compression = 'gzip', sep = '\t', na_rep = 'NaN'))
    exit()

# Expand the annotations to columns
composite[annotations] = composite.apply(lambda x: extract_annotations(x['info']), axis = 1, result_type = 'expand')

# Get the coding type
composite['coding_type'] = composite['annotation'].apply(lambda x: get_coding_type(x))

# If possible, get the VAF
tumour_sample = get_tumour_sample(args.ann)
if type(tumour_sample) == str:
    composite['vaf'] = composite.apply(lambda x: extract_vaf(x['format'], x[tumour_sample.lower()]), axis = 1)
else:
    composite['vaf'] = float('NaN')

# Append to composite matrix
if (ann.shape[0] == phylop.shape[0]) & (ann.shape[0] == composite.shape[0]):
    (composite[basic_columns + phylop_columns + annotations]
        .to_csv(args.out, index = False, compression = 'gzip', sep = '\t', na_rep = 'NaN'))
else:
    raise Exception("The composite matrix does not have the expected number of variants.")