#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip

from vcf import load_vcf, count_spims
from annotation import load_bed

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--samples', type=str)
parser.add_argument('--spm_vcfs', type=str)
parser.add_argument('--sim_vcfs', type=str)
parser.add_argument('--spim_vcfs', type=str)
parser.add_argument('--spm_anns', type=str)
parser.add_argument('--sim_anns', type=str)
parser.add_argument('--spm_phylops', type=str)
parser.add_argument('--sim_phylops', type=str)
parser.add_argument('-o', '--out', type=str)

args = parser.parse_args()

# Make dataframe
files_dict = {'sample':     args.samples.split(),
              'spm_ann':    args.spm_anns.split(),
              'sim_ann':    args.sim_anns.split(),
              'spm_phylop': args.spm_phylops.split(),
              'sim_phylop': args.sim_phylops.split()}

if args.spim_vcfs:
    files_dict['spim_vcf'] = args.spim_vcfs.split()

else:
    files_dict['spm_vcf'] = args.spm_vcfs.split()
    files_dict['sim_vcf'] = args.sim_vcfs.split()
    
            
files = pd.DataFrame(data=files_dict)

# Count SPMs and SIMs
summary = pd.DataFrame(columns = ['sample', 
                                  'spm_vcf', 'sim_vcf',
                                  'spm_ann', 'sim_ann',
                                  'spm_phylop', 'sim_phylop'])

for i, row in files.iterrows():

    if args.spim_vcfs:
        n_vcf_spm, n_vcf_sim = count_spims(load_vcf(row['spim_vcf']))
       
    else:
        n_vcf_spm, n_vcf_sim = (count_spims(load_vcf(row['spm_vcf'])
                                    .append(load_vcf(row['sim_vcf']))))
        

    n_ann_spm, n_ann_sim = (count_spims(load_vcf(row['spm_ann'])
                                .append(load_vcf(row['sim_ann']))))
    n_phylop_spm, n_phylop_sim = (load_bed(row['spm_phylop']).shape[0], 
                                 load_bed(row['sim_phylop']).shape[0])

    summary.loc[i] = [row['sample'], 
                      n_vcf_spm, n_vcf_sim, 
                      n_ann_spm, n_ann_sim, 
                      n_phylop_spm, n_phylop_sim]

# Check if all the file types have the same number of variants
if summary.shape[0] != \
    summary[(summary['spm_vcf'] == summary['spm_ann']) & 
            (summary['spm_vcf'] == summary['spm_phylop'])].shape[0]:

    summary.to_csv(args.out, sep='\t', index=False)
    raise Exception(f"The number of somatic point mutations (SPMs) "\
                    "does not match in all files per sample.")

elif summary.shape[0] != \
    summary[(summary['sim_vcf'] == summary['sim_ann']) & 
            (summary['sim_vcf'] == summary['sim_phylop'])].shape[0]:

    summary.to_csv(args.out, sep='\t', index=False)
    raise Exception("The number of somatic indel mutations (SIMs) "\
                    "does not match in all files per sample.")

else:
    summary.to_csv(args.out, sep='\t', index=False)