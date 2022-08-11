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
    ann_list = [','.join(l) for l in list(zip(*[a.split('|') for a in ann_str.split(',')]))]
    lof_list = [','.join(l) for l in list(zip(*[a.split('|') for a in lof_str.split(',')]))]
    nmd_list = [','.join(l) for l in list(zip(*[a.split('|') for a in nmd_str.split(',')]))]
     
    return [ann_list, lof_list, nmd_list]

def extract_vaf(format, t_info):
     
    ''' Extract the variant allele fraction '''
     
    vaf_index = format.split(':').index('AF')
    vaf = t_info.split(':')[vaf_index]
     
    return vaf

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

# Columns to keep
ann_columns = ['chrom', 'pos', 'ref', 'alt', 'mutation_id']
phylop_columns = ['start', 'stop', 'phylop', 'mutation_id']

# Annotations
annotations = "Allele | Annotation | Annotation_Impact | Gene_Name | "\
              "Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | "\
              "Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | "\
              "AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO".lower().split(' | ')
lofs = ['lof_' + lof for lof in "Gene_Name | Gene_ID | Number_of_transcripts_in_gene | "\
       "Percent_of_transcripts_affected".lower().split(' | ')]
nmds = ['nmd_' + nmd for nmd in "Gene_Name | Gene_ID | Number_of_transcripts_in_gene | "\
       "Percent_of_transcripts_affected".lower().split(' | ')]

# Merge annotations and phylop data based on mutation id
composite = ann.merge(phylop[phylop_columns], on = 'mutation_id', how = 'inner')
    
# Expand the annotations to columns
composite[['ann_list', 'lof_list', 'nmd_list']] = composite.apply(lambda x: extract_annotations(x['info']), axis = 1, result_type = 'expand')
composite[annotations] = pd.DataFrame(composite['ann_list'].tolist(), index = composite.index)
composite[lofs] = pd.DataFrame(composite['lof_list'].tolist(), index = composite.index)
composite[nmds] = pd.DataFrame(composite['nmd_list'].tolist(), index = composite.index)

# If possible, get the VAF
tumour_sample = get_tumour_sample(args.ann)
if type(tumour_sample) == str:
    composite['vaf'] = composite.apply(lambda x: extract_vaf(x['format'], x[tumour_sample.lower()]), axis = 1)
else:
    composite['vaf'] = float('NaN')
ann_columns = ann_columns + ['vaf']

# Append to composite matrix
if (ann.shape[0] == phylop.shape[0]) & (ann.shape[0] == composite.shape[0]):
    (composite[ann_columns + phylop_columns + annotations + lofs + nmds]
        .to_csv(args.out, index = False, compression = 'gzip', sep = '\t'))
else:
    raise Exception("The composite matrix does not have the expected number of variants.")