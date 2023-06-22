#!/bin/bash

# name: 1.1.NCCM_analysis.sh
# Script has two parts:
# PART 1 - Submit NCCM_scan.sh
# In this section, counts of non-coding alterations in the flanking regions of genes (n'hood)
# are tallied and then normalized for length
#-----------------------------------------------------------------------------------------
# PART 2
# Genes which have NCCMs that are extreme (i.e, >=2/kbp) are extracted, 
# followed by the retrieval of their NCCM coordinates
# written by: ShardaS
# adapted: RP
#-----------------------------------------------------------------------------------------
# REQUIRED INPUT FILES - $MASTER
# Composite S[P/I]M matrix file with the following columns:
# Sample , Chrom , Pos,  Constraint score (GERP/PhyloP), CODING/NONCODING
# Give column location of constraint and coding/non-coding as input
# Run: 1.1.NCCM_analysis master_matrix gene_set_with_flanks out_prefix \
# constraint_threshold constr_column type_column (n_parallel_processes)
####################################################################################################################

# Where is the data?
MASTER="$1"
FULL_GENE_SET="$2"
PREF="$3"
# Constraint score threshold
CONS=$4
# In column of master:
CONS_COL=$5
# COding/ non-coding type in column:
TYPE_COL=$6
# Number of parallel processes
PROCESSES=$7

# First split the gene data into as many parts as there are processes
SETS=${FULL_GENE_SET%.in}
split --additional-suffix=.in --numeric-suffixes=1 -n l/$PROCESSES $FULL_GENE_SET ${SETS##*/}-

# Submit the scan as parallel processes
for i in ${SETS##*/}-* ; do 
    workflow/scripts/1.2.NCCM_scan.sh $MASTER $i $PREF-${i##*-} $CONS $CONS_COL $TYPE_COL &
done
wait

# Put together the scans 
printf "GENE\tnonCoding\tconsnonCoding\tconsSamples\tnccmRate\tsamplesRate\n" > $PREF.scan.tsv
cat $PREF-*.scan.tsv >> $PREF.scan.tsv

# get the top nccm_gene + their NCCM rates
#-------------------------------------------

awk -F'\t' 'NR==1{print;next}$5>=2{print | "sort -k5nr"}' $PREF.scan.tsv \
> $PREF.top_nccm_genes.list

echo "NCCM top list ready."

# get the top nccm genes' flanking coordinates
#----------------------------------------------
while read gene nonCoding consnonCoding consSamples nccmRate samplesRate 
do
	awk -F'\t' -v OFS='\t' -v gene=$gene '$1==gene' $FULL_GENE_SET >> $PREF.top_nccm_gene_100kb_flanks.in
done < $PREF.top_nccm_genes.list | tail -n +2

# iterate through the top gene list and get their get the S[P/I]M coordinates
#-----------------------------------------------------------------------------
while read gene chr dflank start stop uflank tlength norm_length
do
	zcat $MASTER | awk -v chr="$chr" -v dflank="$dflank" -v uflank="$uflank" -v type="$TYPE_COL" \
	'$1==chr && $2>=dflank && $2<=uflank && $type=="noncoding"' | 
	awk -v cons="$CONS" -v cons_col="$CONS_COL" '$cons_col!="NaN" && $cons_col>=cons' | cat -n | sed "s/^[ ]*/$gene.NCCM_/g" \
	>> $PREF.top_nccm_genes_spim.tsv
	
done < $PREF.top_nccm_gene_100kb_flanks.in

echo "All done :)"

exit 0