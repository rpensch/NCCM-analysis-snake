#!/bin/bash
MASTER="$1"
GENE_SET="$2"
CONS="$3"
OUTPUT_FILE="$4"

printf "gene\tnc_count\tncncm_count\tnccm_count\tunique_nc_samples\tunique_ncncm_samples\tunique_nccm_samples\n" > $OUTPUT_FILE

bedtools intersect -a $GENE_SET -b $MASTER -wa -wb | \
    awk -F'\t' -v cons_threshold="$CONS" '
    BEGIN { OFS="\t" }
    {
        # Fields from the variant BED file
        sample_id = $9
        cons_score = $8

        # Fields from the gene BED file
        gene = $4
        
        # Increment total non-coding count for the gene
        noncoding_count[gene]++
		nc_samples[gene][sample_id] = 1
        
        # Check conservation
        if (cons_score != "NaN" && cons_score >= cons_threshold) {
            noncoding_constraint_count[gene]++
            nccm_samples[gene][sample_id] = 1
        }
		else if (cons_score != "NaN" && cons_score < cons_threshold) {
			noncoding_nonconstraint_count[gene]++
            ncncm_samples[gene][sample_id] = 1
		}
    }
    END {
        # awk will only know about genes that had at least one variant.
        # To print all genes (even with 0 counts), we need to read the gene set again.
        while ((getline line < "'$GENE_SET'") > 0) {
            split(line, parts, " ")
            gene = parts[4]

            # Get counts, defaulting to 0
            nc_count = noncoding_count[gene] ? noncoding_count[gene] : 0
            nccm_count = noncoding_constraint_count[gene] ? noncoding_constraint_count[gene] : 0
			ncncm_count = noncoding_nonconstraint_count[gene] ? noncoding_nonconstraint_count[gene] : 0
            unique_nc_samples = length(nc_samples[gene])
			unique_nccm_samples = length(nccm_samples[gene])
			unique_ncncm_samples = length(ncncm_samples[gene])

            # Print final result
			printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", gene, nc_count, ncncm_count, nccm_count, unique_nc_samples, unique_ncncm_samples, unique_nccm_samples
        }
    }
' >> "$OUTPUT_FILE"