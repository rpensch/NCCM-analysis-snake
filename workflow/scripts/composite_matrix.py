#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip

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

def load_bed(bed):
     
    ''' Load BED file '''
     
    bed = pd.read_csv(bed, comment='#', sep='\t', header=None, 
                      names=['chrom','start','stop',
                             'phylop','mutation_id'])
     
    return bed

def unique_values(my_list):
     
    ''' Extract unique values from list (still sorted) '''
     
    unique = []
    for i in my_list:
        if i not in unique:
            unique.append(i)
    return unique

def extract_annotations(info):
     
    ''' Extract the annotations from the VCF info field. '''
     
    # Split the info field to list
    split_info = info.split(';')
     
    # Get the list entries for ANN, LOF and NMD
    ann_str = [i.replace('ANN=', '').strip('()') for i in split_info if i.startswith('ANN=')]
    lof_str = [i.replace('LOF=', '').strip('()') for i in split_info if i.startswith('LOF=')]
    nmd_str = [i.replace('NMD=', '').strip('()') for i in split_info if i.startswith('NMD=')]
     
    # Make sure they are not empty and get the string / if empty make an empty entries string
    if len(ann_str) > 0 :
        ann_str = ann_str[0]
    else:
        ann_str = '|'*(len(ann)-1)
     
    if len(lof_str) > 0 :
        lof_str = lof_str[0]
    else:
        lof_str = '|'*(len(lofs)-1)
     
    if len(nmd_str) > 0 :
        nmd_str = nmd_str[0]
    else:
        nmd_str = '|'*(len(nmds)-1)
     
    # Join multiple entries of the same field together and make a list of tuples for ANN, LOF, NMD 
    ann_list = [','.join(filter(None, unique_values(l))) for l in list(zip(*[a.split('|') for a in ann_str.split(',')]))]
    lof_list = [','.join(filter(None, unique_values(l))) for l in list(zip(*[a.split('|') for a in lof_str.split(',')]))]
    nmd_list = [','.join(filter(None, unique_values(l))) for l in list(zip(*[a.split('|') for a in nmd_str.split(',')]))]
     
    return [ann_list, lof_list, nmd_list]

def extract_all_annotations_tmp(info, which):
     
    ''' Extract the annotations from the VCF info field without removing duplicates '''
     
    if which == 'ann':
        col = 1
    elif which == 'gene':
        col = 3
     
    # Split the info field to list
    split_info = info.split(';')
     
    # Get the list entries for ANN, LOF and NMD
    ann_str = [i.replace('ANN=', '').strip('()') for i in split_info if i.startswith('ANN=')]
     
    # Make sure they are not empty and get the string / if empty make an empty entries string
    if len(ann_str) > 0 :
        ann_str = ann_str[0]
    else:
        ann_str = '|'*(len(ann)-1)
     
    # Join multiple entries of the same field together and make a list of tuples for ANN, LOF, NMD 
    ann = [','.join(l) for l in list(zip(*[a.split('|') for a in ann_str.split(',')]))][col]
     
    return ann

# Define coding and non-coding variants
coding = ['conservative_inframe_deletion', 'conservative_inframe_insertion', 'disruptive_inframe_deletion',
          'disruptive_inframe_deletion&splice_region_variant', 'disruptive_inframe_insertion', 'frameshift_variant',
          'frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant',
          'frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant',
          'frameshift_variant&splice_region_variant', 'frameshift_variant&stop_gained',
          'frameshift_variant&stop_lost', 'initiator_codon_variant', 'missense_variant', 
          'missense_variant&splice_region_variant', 'splice_region_variant&synonymous_variant', 'start_lost', 
          'stop_gained', 'stop_gained&conservative_inframe_insertion', 'stop_gained&disruptive_inframe_deletion', 
          'stop_gained&splice_region_variant', 'stop_lost', 'stop_retained_variant','synonymous_variant']

noncoding = ['3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant', '5_prime_UTR_variant',
             'downstream_gene_variant', 'intergenic_region', 'intragenic_variant', 'intron_variant',  
             'non_coding_transcript_exon_variant', 'non_coding_transcript_variant', 'splice_acceptor_variant',
             'splice_acceptor_variant&intron_variant', 'splice_acceptor_variant&splice_donor_variant&intron_variant',
             'splice_donor_variant','splice_donor_variant&intron_variant', 'splice_region_variant&intron_variant', 
             'splice_region_variant&non_coding_transcript_exon_variant',
             'upstream_gene_variant']

ambiguous = ['splice_region_variant']

def get_coding_type(annotation):
     
    ''' Is the variant coding or non-coding? '''
     
    ann_list = annotation.split(',') # split annotations and flatten
    coding_types = []
    for a in ann_list:
        if a in coding:
            coding_types.append('coding')
        elif a in noncoding:
            coding_types.append('noncoding')
        elif a in ambiguous:
            coding_types.append('splice_region_variant')
        else:
            import warnings
            warnings.warn(f"WARNING: Unknown variant type: {a}.")
    
    if 'coding' in coding_types:
        return 'coding'
    elif 'noncoding' in coding_types:
        return 'noncoding'
    elif 'splice_region_variant' in coding_types:
        return 'splice_region_variant'
    else:
        return 'NaN'

def select_best_annotation(annotation, gene, coding_type):
     
    ''' Select the most important annotation for easier downstream processing '''
     
    ann_list = annotation.split(',')
    gene_list = gene.split(',')
    all_annotations = coding if coding_type == 'coding' else noncoding
    if coding_type == 'coding':
        all_annotations = coding
    elif coding_type == 'noncoding':
        all_annotations = noncoding
    elif coding_type == 'splice_region_variant':
        all_annotations = ambiguous
    else:
        return ['NaN', gene]

    for a in all_annotations:
        if a in ann_list:
            my_anno = a
            index = ann_list.index(a)
            my_gene = gene_list[index]
            return [my_anno, my_gene]

def extract_vaf(format, t_info):
     
    ''' Extract the variant allele fraction '''
     
    vaf_index = format.split(':').index('AF')
    vaf = t_info.split(':')[vaf_index]
     
    return vaf

# Annotations
annotations = "Allele | Annotation | Annotation_Impact | Gene_Name | "\
              "Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | "\
              "Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | "\
              "AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO".lower().split(' | ')
lofs = ['lof_' + lof for lof in "Gene_Name | Gene_ID | Number_of_transcripts_in_gene | "\
       "Percent_of_transcripts_affected".lower().split(' | ')]
nmds = ['nmd_' + nmd for nmd in "Gene_Name | Gene_ID | Number_of_transcripts_in_gene | "\
       "Percent_of_transcripts_affected".lower().split(' | ')]
# Columns to keep
ann_columns = ['chrom', 'pos', 'ref', 'alt'] + ['sample'] + ['vaf'] + ['coding_type']
phylop_columns = ['start', 'stop', 'phylop', 'mutation_id']

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

# If there are no mutations, save empty and exit
if (ann.shape[0] == 0) & (phylop.shape[0] == 0) & (composite.shape[0] == 0):
    composite[ann_columns + phylop_columns + ['main_annotation', 'main_gene']+ annotations + lofs + nmds] = 0
    (composite[ann_columns + phylop_columns + ['main_annotation', 'main_gene']+ annotations + lofs + nmds]
        .to_csv(args.out, index = False, compression = 'gzip', sep = '\t'))
    exit()

# Expand the annotations to columns
composite[['ann_list', 'lof_list', 'nmd_list']] = composite.apply(lambda x: extract_annotations(x['info']), axis = 1, result_type = 'expand')
composite[annotations] = pd.DataFrame(composite['ann_list'].tolist(), index = composite.index)
composite[lofs] = pd.DataFrame(composite['lof_list'].tolist(), index = composite.index)
composite[nmds] = pd.DataFrame(composite['nmd_list'].tolist(), index = composite.index)

# Get the coding type
composite['coding_type'] = composite['annotation'].apply(lambda x: get_coding_type(x))

# Select my favourite annotations
composite['ann_tmp'] = composite.apply(lambda x: extract_all_annotations_tmp(x['info'], 'ann'), axis = 1)
composite['gene_tmp'] = composite.apply(lambda x: extract_all_annotations_tmp(x['info'], 'gene'), axis = 1)
composite[['main_annotation','main_gene']] = composite.apply(lambda row: select_best_annotation(row['ann_tmp'], row['gene_tmp'], row['coding_type']), axis = 1, result_type = 'expand')

# If possible, get the VAF
tumour_sample = get_tumour_sample(args.ann)
if type(tumour_sample) == str:
    composite['vaf'] = composite.apply(lambda x: extract_vaf(x['format'], x[tumour_sample.lower()]), axis = 1)
else:
    composite['vaf'] = float('NaN')

# Append to composite matrix
if (ann.shape[0] == phylop.shape[0]) & (ann.shape[0] == composite.shape[0]):
    (composite[ann_columns + phylop_columns + ['main_annotation', 'main_gene']+ annotations + lofs + nmds]
        .to_csv(args.out, index = False, compression = 'gzip', sep = '\t'))
else:
    raise Exception("The composite matrix does not have the expected number of variants.")