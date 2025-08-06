#!/bin/bash

GENE_SET=$1
QC_BED=$2
OUT=$3

bedtools intersect -sorted -wao -a $GENE_SET -b $QC_BED |
awk '
BEGIN { OFS="\t" }
{
    # Calculate section overlap
    gene_name = $4
    overlap[gene_name] += $NF
    

    # Store and calculate gene info
    gene_info[gene_name] = $1"\t"$2"\t"$3
    gene_length[gene_name] = $3 - $2
    
}
END {
    # Iterate through genes again
    for (gene in overlap) {

        # Calculate total bp overlap
        bp_overlap = overlap[gene]
        gene_len = gene_length[gene]

        percent_overlap = (bp_overlap / gene_len) * 100

        # Print the final result
        printf "%s\t%s\t%f\n", gene_info[gene], gene, percent_overlap
    }
}' | sort -k1,1 -k2,2n > $OUT