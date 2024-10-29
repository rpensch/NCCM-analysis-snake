# Convert the vcf files to bed format
rule vcf2bed:
    input:
        expand("resources/input/{{sample}}.{spim}.vcf.gz", spim = spims)
    output:
        expand("results/input_beds/{{sample}}.{spim}.sorted.bed", spim = spims)
    shell:
        "workflow/scripts/vcf2bed.sh {wildcards.sample} {input} {output}"

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