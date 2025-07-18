
# NCCM-analysis pipeline

## Input 

Filtered somatic variants in gzip compressed vcf file. For each sample, either provide one vcf that includes both somatic point mutations (SPIMs) and somatic indel mutations (SIMs), or two vcfs that are somatic point (SPM) and indel mutation (SIM) separated (one SPM file and one SIM file as e.g. from Mutect2 output). List the input files in an tsv file with with header that is either `sample spim` or `sample spm sim` - see `config/example.input.tsv`. Provide absolute file paths.

## Required resources

- PhyloP scores separated by chromosome in gzipped bed files. Columns (no header): `chromosome start end id phyloP`. The files have to be named as follows: chr1.bed.gz, chr2.bed.gz, chr3.bed.gz

- A list of chromosomes for your genome, e.g. `resources/canfam4.chromosomes.txt`

- A set of genes to test for NCCM enrichment in bed format, with the following columns (no header): `chromosome lower_flank upper_flank gene`. See `resources/canfam4_gene_100kb_flanks_4_NCCM_v2.bed`

Make sure the chromosome notation is consistent across all input and resource data (either chr1 or 1)! Including vcf files, the phyloP scores bed file content + the naming of the phylP score bed files, the chromosome list file and the gene flanks file. 

## Config

In the config file `config/config.yaml`, change `vcfs: "config/example.input.tsv"` to the name of your input tsv file. Also, set `spim:` to either  `"separated"` or `"combined"`, depending on whether you provided one or two vcf files.

- `start_from`: `"vcf"` or `"matrix"`

    - depending on `start_from` either of the following two must be set:

    - `vcfs`: path to the input tsv file described above

    - `matrix`: path to a pre-existing matrix

- `genome`: Genome input for [snpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/) - e.g. GRCh37.75

- `phyloP`: Path for the directory of phyloP score files

- `phyloP_threshold`: The threshold for what is considered constraint - we have used 1.2 for human and 1.3 for dogs (8% of the genome)

- `chrom_list`: Path for file with list of chromosomes, see above in *Required resources*

- `gene_set`: Path of the gene_set file as described in *Required resources*

## Software

Tested with the followinng software versions. 

- Snakemake (version 7.8.5)
- BEDOPS (version 2.4.39)
- snpEff (version 4.3t)
- python3 (version 3.9.5)
- BEDTools (version 2.29.2)

On Uppmax do `module load bioinfo-tools snakemake/7.8.5 BEDOPS/2.4.39 snpEff/4.3t python3/3.9.5 BEDTools/2.29.2`

## Run

Test/dry run:

`snakemake -np all`

Run the whole workflow with e.g. 16 cores:

`snakemake --cores 16 all`

Create the matrix, don't run the nccm analysis:

`snakemake --cores 16 composite_matrix`

Run only nccm analysis with a pre-existing matrix:

`snakemake --cores 16 nccms`

## Output

Main output files:

- `results/composite_matrix/*.composite_matrix.tsv.gz` - includes all annotation data (snpeff + phyloP) for all samples. Missing phylop cores are `NaN`. 

- `results/nccms/*scan.tsv` - this is the main output. 
