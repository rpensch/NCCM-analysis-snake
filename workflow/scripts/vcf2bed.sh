#!/bin/bash

sample=$1
spm=$2
sim=$3
spm_out=$4
sim_out=$5

# Convert spm vcfs to bed files
zcat $spm | awk -v OFS='\t' -v sample="$sample" '{print $0,sample":"$1":"$2":"$4":"$5}' |
convert2bed -i vcf --snvs | awk -v OFS='\t' '{print $1,$2,$3,$(NF)}' | sort -k1,1 -k2,2n \
>> $spm_out

# Convert sim vcfs to bed files
zcat $sim | awk -v OFS='\t' -v sample="$sample" '{print $0,sample":"$1":"$2":"$4":"$5}' |
convert2bed -i vcf --insertions | awk -v OFS='\t' '{print $1,$2,$3,$(NF)}' | sort -k1,1 -k2,2n \
>> $sim_out.tmp
zcat $sim | awk -v OFS='\t' -v sample="$sample" '{print $0,sample":"$1":"$2":"$4":"$5}' |
convert2bed -i vcf --deletions |  awk -v OFS='\t' '{print $1,$2,$3,$(NF)}' | sort -k1,1 -k2,2n \
>> $sim_out.tmp
# Handle MNPs manually
zgrep -v '^#' $sim | 
awk -v OFS='\t' -v sample="$sample" '{ if (length($4)==length($5)) { print $0,sample":"$1":"$2":"$4":"$5} }' |
while read chrom pos id ref alt qual filter info format tumor normal mutation_id ; do
    
    len=`echo $ref | awk '{print length($1)}'`
    printf "$chrom\t$(($pos-1))\t$(($pos-1+$len))\t$mutation_id\n" >> $sim_out.tmp

done

sort -k1,1 -k2,2n $sim_out.tmp > $sim_out
rm $sim_out.tmp

n_spm=$(awk 'END{print NR}' $spm_out)
n_sim=$(awk 'END{print NR}' $spm_out)

if [ $n_spm -eq 0 ] ; then 
    echo "0 somatic point mutations (SPMs) in file."
    exit 1
elif [ $n_sim -eq 0 ] ; then 
    echo "0 somatic indel mutations (SIMs) in file."
    exit 1
else
    exit 0
fi