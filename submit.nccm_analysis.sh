#!/bin/bash -l
#SBATCH -J nccm_snake
#SBATCH -p core -n 16
#SBATCH -t 3-00:00:00
#SBATCH -A naiss2024-5-39
#SBATCH --mail-user raphaela.pensch@imbim.uu.se
#SBATCH --mail-type=END

module load bioinfo-tools snakemake/7.8.5 BEDOPS/2.4.39 snpEff/4.3t python3/3.9.5 BEDTools/2.29.2

snakemake --cores 16 all
