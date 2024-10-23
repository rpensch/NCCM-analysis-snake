# NCCM-analysis pipeline

## Input 

Filtered somatic variants in gzip compressed vcf file. For each sample, either provide one vcf that includes both somatic point mutations (SPIMs) and somatic indel mutations (SIMs), or two vcfs that are somatic point (SPM) and indel mutation (SIM) separated (one SPM file and one SIM file as e.g. from Mutect2 output). List the input files in an tsv file with with header that is either `sample spim` or `sample spm sim` - see `config/example.input.tsv`

## Required resources

- PhyloP scores separated by chromosome in gzipped bed files. Columns (no header): `chromosome start end id phyloP`. The files have to be named as follows: chr1.bed.gz, chr2.bed.gz, chr3.bed.gz

- A list of chromosomes for your genome, e.g. `resources/canfam4.chromosomes.txt`

- A set of genes to test for NCCM enrichment, with the following columns (no header): `gene_name chromosome lower_flank gene_start gene_end upper flank cds_size_bp correction_factor`. The lower and upper flank are 100 Kbp from the gene start and end. The correction factor is `100,000 / (upper - lowe flank - cds_size)`. See `resources/canfam4_gene_100kb_flanks_4_NCCM.in`

## Config

In the config file `config/config.yaml`, change `vcfs: "config/example.input.tsv"` to the name of your input tsv file. Also, set `spim:` to either  `"separated"` or `"combined"`, depending on whether you provided one or two vcf files.

- `genome`: Genome input for [snpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/)

- `phyloP`: Path for the directory of phyloP score files

- `phyloP_threshold`: The threshold for what is considered constraint - we have used 1.2 for human and 1.3 for dogs (8% of the genome)

- `chrom_list`: Path for file with list of chromosomes, see above in `Require resources`

- `gene_set`: Path of the gene_set file as described in `Required resources`

## Software

Tested with the followinng software versions. 

- Snakemake (version 7.8.5)
- BEDOPS (version 2.4.39)
- snpEff (version 4.3t)
- vcftools (version 0.1.12)
- python3 (version 3.9.5)
- BEDTools (version 2.29.2)

On Uppmax `module load bioinfo-tools snakemake/7.8.5 BEDOPS/2.4.39 snpEff/4.3t vcftools/0.1.12 python3/3.9.5 BEDTools/2.29.2`

## Output

Main output files:

- `results/composite_matrix/*.composite_matrix.tsv.gz` - includes all annotation data (snpeff + phyloP) for all samples

- `results/nccms/*scan.tsv` - this is the main output. 

    - `nonCoding` for all non-coding mutations
    - `consnonCoding` for non-coding constraint mutations (NCCMs)
    - `consSampels` for the number of unique samples with NCCMs
    - `nccmRate` is "the number of NCCMs per 100 Kbp per gene", basically `consnonCoding * the correction factor from the gene_set file` 
    - `samplesRate` is `consSamples * the correction factor from the gene_set file` 

- `results/nccms/*top_nccm_genes_spim.tsv` - a list of annotated mutations for each gene extracted from the composite matrix. This file can become quite big for large datasets. 


To rank genes for downstream analyses, we use `nccmRate / the number of samples in the cohort`.  