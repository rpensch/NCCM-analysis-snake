#!/bin/bash
#SBATCH -J nccm_snake
#SBATCH -p core -n 16
#SBATCH -t 1-00:00:00
#SBATCH -A sens2017503
#SBATCH --mail-user raphaela.pensch@imbim.uu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools snakemake/7.8.5 BEDOPS/2.4.39 snpEff/4.3t python3/3.9.5 BEDTools/2.31.1

matrix_path=$1
matrix_name=$(basename $matrix_path)
run_name=${matrix_name%.composite_matrix.tsv.gz}

snakemake --cores 16 nccms \
--config gene_set=/proj/sens2017503/nobackup/12_pancancer_RP/1_new_flanks_file/a_create/hg19_gene_10kb_flanks.bed \
run_name=$run_name.10kbp_flanks \
matrix_path=$matrix_path

snakemake --cores 16 nccms \
--config gene_set=/proj/sens2017503/nobackup/12_pancancer_RP/1_new_flanks_file/a_create/hg19_gene_50kb_flanks.bed \
run_name=$run_name.50kbp_flanks \
matrix_path=$matrix_path

snakemake --cores 16 nccms \
--config gene_set=/proj/sens2017503/nobackup/12_pancancer_RP/1_new_flanks_file/a_create/hg19_gene_5kb_flanks.bed \
run_name=$run_name.5kbp_flanks \
matrix_path=$matrix_path

snakemake --cores 16 nccms \
--config gene_set=/proj/sens2017503/nobackup/12_pancancer_RP/1_new_flanks_file/a_create/hg19_gene_100kb_flanks.bed \
run_name=$run_name.100kbp_flanks \
matrix_path=$matrix_path