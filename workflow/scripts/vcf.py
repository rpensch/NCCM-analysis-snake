#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip

''' Functions to load and process vcf files '''

def get_vcf_header(vcf_path):
     
    ''' Get the vcf header '''
     
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:  
            if line.startswith("#CHROM"):
                vcf_names = [n.strip('\n') for n in line.split('\t')]
                break
    return vcf_names

def get_tumour_sample(vcf_path):
     
    ''' Get which sample / columns is the tumor vs. the normal '''
    
    tumour_sample = float('NaN')
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:  
            if line.startswith('##tumor_sample='):
                tumour_sample = line.replace('##tumor_sample=', '').strip()
                break
     
    return tumour_sample

def load_vcf(vcf_path):
         
    ''' Load vcf file '''
     
    vcf_names = get_vcf_header(vcf_path)
    vcf = pd.read_csv(vcf_path, comment='#', sep='\t', header=None, names=vcf_names)
     
    return vcf

def count_spims(vcf):
     
    ''' Count the number of spms and sims '''
     
    n_spm = vcf[(vcf['ALT'].str.len()==1) & 
                (vcf['REF'].str.len()==1)].shape[0]
     
    n_sim = vcf[(vcf['ALT'].str.len()!=1) | 
                (vcf['REF'].str.len()!=1)].shape[0]
     
    return n_spm, n_sim 