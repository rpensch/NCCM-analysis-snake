#!/usr/bin/env python3

''' Functions to Extract annotation data from vcf file '''

import pandas as pd

def load_bed(bed):
     
    ''' Load BED file '''
     
    bed = pd.read_csv(bed, comment='#', sep='\t', header=None, 
                      names=['chrom','start','stop',
                             'phylop','mutation_id'])
     
    return bed

def extract_annotations(info, ann):
     
    ''' Extract the annotations from the VCF info field. '''
     
    # Split the info field to list
    split_info = info.split(';')
     
    # Get the list entries for ANN, LOF and NMD
    ann_str = [i.replace('ANN=', '') for i in split_info if i.startswith('ANN=')]
     
    # Make sure they are not empty and get the string / if empty make an empty entries string
    if len(ann_str) > 0 :
        ann_str = ann_str[0]
    else:
        ann_str = '|'*(len(ann)-1)
     
    # Select the first, highest impact annotation and split fields
    ann_list = ann_str.split(',')[0].split('|')

    return ann_list

# Define coding and non-coding variants
coding = ['HIGH', 'MODERATE', 'LOW']

noncoding = ['MODIFIER']

def get_coding_type(annotation_impact):
     
    ''' Is the variant coding or non-coding? '''

    if annotation_impact in coding:
        return 'coding'
    elif annotation_impact in noncoding:
        return 'noncoding'
    else:
        import warnings
        warnings.warn(f"WARNING: Unknown variant type: {annotation_impact}.")
        return 'NaN'

def extract_vaf(format, t_info):
     
    ''' Extract the variant allele fraction '''
     
    vaf_index = format.split(':').index('AF')
    vaf = t_info.split(':')[vaf_index]
     
    return vaf
