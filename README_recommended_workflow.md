# NCCM analysis as recommended by ME :) 

## 1. Annotate vcf files and create the matrix

If you are starting from vcf files, you need to first annotate and combine them into a data matrix. 

`snakemake --cores 16 composite_matrix`

## 2. Define your regions of interest

Which regions do you want to test for NCCM enrichment?

Create the `gene_set` file, as described in the README. This can be any regions, just make sure all the required columns are there. 

I'm running two separate analysis:

- Gene-centric with variable flank sizes (100 kbp, 50 kbp, 10 kbp, 5 kbp)

- Window-based with 2 kbp sliding windows across the genome that have 50 % overlap. Create those with bedtools. 

For the gene-centric, this means that you have to run all the analyses for each individual flank file. Most importantly, you need covariate data for each of the flank files.

## 3. Create supporting data for the regression

First, we need to know for each gene set file, how many non-coding constraint and non-coding non-constraint positions there are in each region. Also we need to know how many non-coding positions don't have phyloP annotations available at all. To learn this, you need run the following:

```snakemake --cores 16 phylop_counts```

Second, we want to know what the overlap of the gene set file with regions of poor quality is. I have used UMAP k100 mappability scores from <https://bismap.hoffmanlab.org/> and filtered them for values < 1. 

So run:

```snakemake --cores 16 gene_qc```

For both, don't forget to set `gene_set` and `qc_bed` in the config. You can also set them on the command line with e.g.:

```snakemake --cores 16 phylop_counts --config gene_set=/your/gene/set/path```

Save these files for later, we will need them for the regression.

Third, we need covariate data for all of the regions. Some of them you need to create yourself, such as GC content, etc. But for some, when you want to use the amount of overlap with another dataset as your covariate (e.g. ATAC-peaks), you can also use the function above. Just replace the `qc_bed` path with a path to the bed file you want to intersect with. 

## 4. Get the NCCM counts

For each region in the `gene_set` file, we want to know how many NCCMs and NCNCMs there are. To get that, run as previously:

`snakemake --cores 16 all` 

OR with a pre-existing matrix:

`snakemake --cores 16 nccms`

Remember again: You will need to run this step separately for `gene_set` files of all the flank sizes you want to test, e.g. 100 kbp, 50 kbp, 10 kbp, 5 kbp.

## What I think you need to prepare before you start:

- Vcf files or matrix

- PhyloP scores

- A list of chromosomes

- Gene set files (NEW format!)

- Coding regions file

- Gene qc data: Mainly a bed file with regions of poor mappability. 