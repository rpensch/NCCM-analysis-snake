# Convert the vcf files to bed format
rule vcf2bedspm:
    input: "resources/input/{sample}.spm.vcf.gz"
    output: "results/input_beds/{sample}.spm.sorted.bed"
    shell:
        """
        zcat {input} | awk -v FS='\t' -v OFS='\t' -v sample='{wildcards.sample}' 'BEGIN {{ }} /^#/ {{print; next}}{{print $1,$2,$3,$4,$5,$6,$7,sample":"$1":"$2":"$4":"$5}}' | 
        convert2bed -i vcf --snvs | awk -v OFS='\t' '{{print $1,$2,$3,$(NF)}}' | sort -k1,1 -k2,2n \
        >> {output}
        """

rule vcf2bedsim:
    input: "resources/input/{sample}.sim.vcf.gz"
    output: "results/input_beds/{sample}.sim.sorted.bed"
    shell:
        """
        zcat {input} | awk -v FS='\t' -v OFS='\t' -v sample='{wildcards.sample}' 'BEGIN {{ }} /^#/ {{print; next}}{{print $1,$2,$3,$4,$5,$6,$7,sample":"$1":"$2":"$4":"$5}}' | 
        convert2bed -i vcf --insertions | awk -v OFS='\t' '{{print $1,$2,$3,$(NF)}}' | sort -k1,1 -k2,2n \
        >> {output}.tmp
        zcat {input} | awk -v FS='\t' -v OFS='\t' -v sample='{wildcards.sample}' 'BEGIN {{ }} /^#/ {{print; next}}{{print $1,$2,$3,$4,$5,$6,$7,sample":"$1":"$2":"$4":"$5}}' | 
        convert2bed -i vcf --deletions |  awk -v OFS='\t' '{{print $1,$2,$3,$(NF)}}' | sort -k1,1 -k2,2n \
        >> {output}.tmp
        # Handle MNPs manually
        zgrep -v '^#' {input} | 
        awk -v FS='\t' -v OFS='\t' -v sample='{wildcards.sample}' '{{ if (length($4)==length($5)) {{print $1,$2,$3,$4,$5,$6,$7,sample":"$1":"$2":"$4":"$5}} }}' | 
        while read chrom pos id ref alt qual filter info format tumor normal mutation_id ; do 
            len=`echo $ref | awk '{{print length($1)}}'` 
            printf "$chrom\t$(($pos-1))\t$(($pos-1+$len))\t$mutation_id\n" >> {output}.tmp
            done
        sort -k1,1 -k2,2n {output}.tmp > {output}
        rm {output}.tmp
        """

# Annotate the bed files with phyloP scores
if config["phyloP_format"] == "bw":

    rule phyloP:
        input:
            bed = "results/input_beds/{sample}.{spim}.sorted.bed",
            phyloP = config["phyloP"]
        output:
            "results/annotation/phyloP/{sample}.{spim}.phyloP.bed"
        shadow: "full"
        shell:
            """
            mkdir -p results/annotation/ ; mkdir -p results/annotation/phyloP
            workflow/scripts/phylop_annotation_bw.sh {input.bed} {output} {input.phyloP}
            """

elif config["phyloP_format"] == "bed":

    with open(config['chrom_list']) as f:
        chromosomes = f.read().splitlines()

    rule split_chroms:
        input:  "results/input_beds/{sample}.{spim}.sorted.bed"
        output: temp("results/input_beds/{sample}.{spim}.sorted.{chrom}.bed")
        shell: 
            """
            awk -v chr='{wildcards.chrom}' -v OFS='\t' '$1==chr' {input} > {output}
            """

    rule phyloP:
        input: 
            bed = "results/input_beds/{sample}.{spim}.sorted.{chrom}.bed",
            phyloP = config['phyloP'] + "/" + "{chrom}.bed.gz"
        output: temp("results/annotation/phyloP/{sample}.{spim}.sorted.{chrom}.phyloP.bed")
        shell:
            """
            mkdir -p results/annotation/ ; mkdir -p results/annotation/phyloP
            bedtools intersect -a {input.bed} -b {input.phyloP} -loj -sorted |
                awk -F '\t' '{{OFS=FS}}{{if ($8 ==".") {{$9="NaN"}}}}1' |
                cut -f 1-4,9 > {output}
            """
    
    rule join_phyloP:
        input: 
            phylop = expand(["results/annotation/phyloP/{{sample}}.{{spim}}.sorted.{chrom}.phyloP.bed"], chrom = chromosomes),
            bed = "results/input_beds/{sample}.{spim}.sorted.bed"
        output: "results/annotation/phyloP/{sample}.{spim}.phyloP.bed"
        wildcard_constraints:                                                                                                                                                            
            spim = "s.m"
        shell:
            """
            workflow/scripts/join_phylop_scores.py --phylop '{input.phylop}' --bed {input.bed} --out {output}
            """