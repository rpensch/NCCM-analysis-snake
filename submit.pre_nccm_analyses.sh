#!/bin/bash -l
#SBATCH -J nccm_snake
#SBATCH -p core -n 16
#SBATCH -t 3-00:00:00
#SBATCH -A naiss2024-5-39
#SBATCH --mail-user raphaela.pensch@imbim.uu.se
#SBATCH --mail-type=END

module load bioinfo-tools snakemake/7.8.5 BEDOPS/2.4.39 snpEff/4.3t python3/3.9.5 BEDTools/2.31.1

snakemake --cores 16 phylop_counts
snakemake --cores 16 gene_ncc_trinucs
snakemake --cores 8 gene_qc --config qc_bed=/proj/sens2017503/nobackup/12_pancancer_RP/1_new_flanks_file/e_gene_filtering/2_mappability/k100.umap.lt1.sorted.merged.bed.gz
snakemake --cores 16 gene_qc --config qc_bed=/proj/sens2017503/nobackup/12_pancancer_RP/1_new_flanks_file/e_gene_filtering/1_hg19_vs_hg38/hg38_diff.contig_dropped.sorted.bed

