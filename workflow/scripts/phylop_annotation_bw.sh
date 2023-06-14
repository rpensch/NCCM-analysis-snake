#!/bin/bash

# Name: phyloP_annotation_bw.sh
# by Raphaela Pensch

# Run: 
# phyloP_annotation_bw.sh variant_data_dir_abs_path file_list phylop_location target_dir out_suffix

# Annotate variants in bed files with phyloP scores

#-------------------------------------------------

# BIGWIG_RS - Gerp schema in bigwig format

# Define path variables
INPUT="$1"
OUTPUT="$2"
BIGWIG_RS="$3"


# Annotate SPM files
while read chr start end mutation_id ; do

    scripts/bigWigToBedGraph $BIGWIG_RS -chrom=$chr -start=$start -end=$end phylop.temp

    # Chose the highest phylop score (for indels that have multiple)
    PHYLOP=$(cut -f4 phylop.temp | sort -nr | head -1 )
    if [ -z "$PHYLOP" ] ; then
        PHYLOP="NaN"
    fi
    printf "$chr\t$start\t$end\t$PHYLOP\t$mutation_id\n" >> $OUTPUT

done < $INPUT

exit 0