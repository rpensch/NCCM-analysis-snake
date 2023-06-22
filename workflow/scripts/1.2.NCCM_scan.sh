#!/bin/bash

# Where is the data?
MASTER="$1"
GENE_SET="$2"
PREF="$3"
# Conservation score threshold
CONS=$4
# In column of master:
CONS_COL=$5
# COding/ non-coding type in column:
TYPE_COL=$6

while read gene chr dflank start stop uflank tlength norm_length ; do

	# non-coding counts for both raw as well as constrained posi
	#-------------------------------------------------------------
	nonCoding=`zcat $MASTER | awk -v chr="$chr" -v dflank="$dflank" -v uflank="$uflank" -v type="$TYPE_COL" \
	'$1==chr && $2>=dflank && $2<=uflank && $type=="noncoding"'|
	wc -l`

	consnonCoding=`zcat $MASTER | awk -v chr="$chr" -v dflank="$dflank" -v uflank="$uflank" -v type="$TYPE_COL" \
	'$1==chr && $2>=dflank && $2<=uflank && $type=="noncoding"' | 
	awk -v cons="$CONS" -v cons_col="$CONS_COL" '$cons_col!="NaN" && $cons_col>=cons' | wc -l`

	consSamples=`zcat $MASTER | awk -v chr="$chr" -v dflank="$dflank" -v uflank="$uflank" -v type="$TYPE_COL" \
	'$1==chr && $2>=dflank && $2<=uflank && $type=="noncoding"' | 
	awk -v cons="$CONS" -v cons_col="$CONS_COL" '$cons_col!="NaN" && $cons_col>=cons' | cut -f5 | sort | uniq | wc -l`

	# compute the rates of NCCCM with $consnonCoding * $norm_length
	# CODE --> value=`echo "scale=4; $p * $q" | bc`
	#------------------------------------------------------
	nccmRate=`echo "scale=20; $consnonCoding * $norm_length" | bc`
    samplesRate=`echo "scale=20; $consSamples * $norm_length" | bc`
	
	printf "$gene\t$nonCoding\t$consnonCoding\t$consSamples\t$nccmRate\t$samplesRate\n" \
	>> $PREF.scan.tsv

done < $GENE_SET

if [ -f $PREF.scan.tsv ] ; then
    echo "Submission $PREF: Scan completed."
    exit 0
else 
    echo "Submission $PREF: Scan failed. exiting" >&2
    exit 1
fi