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
>> $sim_out
zcat $sim | awk -v OFS='\t' -v sample="$sample" '{print $0,sample":"$1":"$2":"$4":"$5}' |
convert2bed -i vcf --deletions |  awk -v OFS='\t' '{print $1,$2,$3,$(NF)}' | sort -k1,1 -k2,2n \
>> $sim_out

n_spm=$(cat $spm_out | wc -l)
n_sim=$(cat $sim_out | wc -l)

if [ $n_spm -eq 0 ] ; then 
    echo "0 somatic point mutations (SPMs) in file."
    exit 1
elif [ $n_sim -eq 0 ] ; then 
    echo "0 somatic indel mutations (SIMs) in file."
    exit 1
else
    exit 0
fi