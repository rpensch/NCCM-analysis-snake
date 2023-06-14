# pipe-those-NCCMs

# NCCM-analysis pipeline

## Input 

Filtered variants in gzip compressed vcf file. Either one vcf per sample that includes both somatic point mutations (SPMs) and somatic indel mutations (SIMs). Or two vcfs that are somatic point and indel mutation (SPIM) separated (one SPM file and one SIM file as e.g. mutect2 output).

Define samples and their vcfs in the config file. Table with header that is either:

sample  spim

or 

sample  spm     sim

## Software

- Snakemake (version 7.8.5)
- BEDOPS (version 2.4.39)
- snpEff (version 4.3t)
- vcftools (version 0.1.12)
- python3 (version 3.9.5)

module load bioinfo-tools snakemake/7.8.5 BEDOPS/2.4.39 snpEff/4.3t vcftools/0.1.12 python3/3.9.5

## Other stuff

profile with https://github.com/Snakemake-Profiles/slurm

bigwigtobedgraph from UCSC




snakemake --profile profiles/slurm.uppmax composite_matrix 