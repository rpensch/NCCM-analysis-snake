# pipe-those-NCCMs

# NCCM-analysis pipeline

## Input 

Filtered somatic variants in gzip compressed vcf file. For each sample, either provide one vcf that includes both somatic point mutations (SPIMs) and somatic indel mutations (SIMs), or two vcfs that are somatic point (SPM) and indel mutation (SIM) separated (one SPM file and one SIM file as e.g. from Mutect2 output).List the input files in an tsv file with with header that is either `sample spim` or `sample spm sim` - see `config/example.input.tsv`

## Required resources

- PhyloP scores separated by chromosome in gzipped bed files. Columns (no header): `chromosome start end id phyloP`

- A set of genes to test for NCCM enrichment, with the following columns (no header): `gene_name chromosome lower_flank gene_start gene_end upper flank cds_size_bp correction_factor`. The lower and upper flank are 100 Kbp from the gene start and end. The correction factor is `100,000 / (upper - lowe flank - cds_size)`. 

## Config

In the config file `config/config.yaml`, change `vcfs: "config/example.input.tsv"` to the name of your input tsv file. Also, set `spim:` to either  `"separated"` or `"combined"`, depending on whether you provided one or two vcf files.

- `genome`: Genome input for [snpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/)

- `phyloP`: Directory of phyloP scores

- `phyloP_threshold`: The threshold for what is considered constraint - we have used 1.2 for human and 1.3 for dogs (8% of the genome)

- `chrom_list`: A list of chromosomes with phyloP scores

- `gene_set`: The location of the gene_set file as described in `Required resources`

## Software

Tested with the followinng software versions. 

- Snakemake (version 7.8.5)
- BEDOPS (version 2.4.39)
- snpEff (version 4.3t)
- vcftools (version 0.1.12)
- python3 (version 3.9.5)
- BEDTools (version 2.29.2)

On Uppmax `module load bioinfo-tools snakemake/7.8.5 BEDOPS/2.4.39 snpEff/4.3t vcftools/0.1.12 python3/3.9.5 BEDTools/2.29.2`
