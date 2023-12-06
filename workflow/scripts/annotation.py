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
coding = ['conservative_inframe_deletion', 'conservative_inframe_insertion', 'disruptive_inframe_deletion',
          'disruptive_inframe_deletion&splice_region_variant', 'disruptive_inframe_insertion', 'frameshift_variant',
          'frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant',
          'frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant',
          'frameshift_variant&splice_region_variant', 'frameshift_variant&stop_gained',
          'frameshift_variant&stop_lost', 'initiator_codon_variant', 'missense_variant', 
          'missense_variant&splice_region_variant', 'splice_region_variant&synonymous_variant', 'start_lost', 
          'stop_gained', 'stop_gained&conservative_inframe_insertion', 'stop_gained&disruptive_inframe_deletion', 
          'stop_gained&splice_region_variant', 'stop_lost', 'stop_retained_variant','synonymous_variant',
          'splice_acceptor_variant', 'splice_acceptor_variant&intron_variant', 
          'splice_acceptor_variant&splice_donor_variant&intron_variant',
          'splice_donor_variant','splice_donor_variant&intron_variant', 'splice_region_variant']

noncoding = ['3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant', '5_prime_UTR_variant',
             'downstream_gene_variant', 'intergenic_region', 'intragenic_variant', 'intron_variant',  
             'non_coding_transcript_exon_variant', 'non_coding_transcript_variant', 'splice_region_variant&intron_variant', 
             'splice_region_variant&non_coding_transcript_exon_variant',
             'upstream_gene_variant']

def get_coding_type(annotation):
     
    ''' Is the variant coding or non-coding? '''

    if annotation in coding:
        return 'coding'
    elif annotation in noncoding:
        return 'noncoding'
    else:
        import warnings
        warnings.warn(f"WARNING: Unknown variant type: {annotation}.")
        return 'NaN'

def extract_vaf(format, t_info):
     
    ''' Extract the variant allele fraction '''
     
    vaf_index = format.split(':').index('AF')
    vaf = t_info.split(':')[vaf_index]
     
    return vaf
