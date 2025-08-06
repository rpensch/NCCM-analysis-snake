# Non-coding (non)constraint positions per gene
rule subtract_coding:
    input:
        coding="resources/coding.{chrom}.bed",
        flanks="resources/gene_set.{chrom}.bed"
    output:
        temp("resources/{chrom}.noncoding_flanks.bed")
    shell:
        """
        bedtools merge -i {input.coding} |
        bedtools subtract -a {input.flanks} -b stdin |
        sort -k1,1 -k2,2n > {output}
        """

rule constraint_pos:
    input:
        phyloP = config['phyloP'] + "/" + "{chrom}.bed.gz"
    output:
        temp("resources/constraint_scores/{chrom}.constraint_positions.bed.gz")
    params:
        phyloP_threshold = config["phyloP_threshold"]
    shell:
        """
        zcat {input.phyloP} | awk -v phylop="{params.phyloP_threshold}" '$5>=phylop' | gzip -c \
        > {output}
        """

rule nonconstraint_pos:
    input:
        phyloP = config['phyloP'] + "/" + "{chrom}.bed.gz"
    output:
        temp("resources/nonconstraint_scores/{chrom}.nonconstraint_positions.bed.gz")
    params:
        phyloP_threshold = config["phyloP_threshold"]
    shell:
        """
        zcat {input.phyloP} | awk -v phylop="{params.phyloP_threshold}" '$5<phylop' | gzip -c \
        > {output}
        """

rule count_constraint_pos:
    input:
        nc_flanks="resources/{chrom}.noncoding_flanks.bed",
        constraint="resources/constraint_scores/{chrom}.constraint_positions.bed.gz"
    output:
        temp("resources/{chrom}.noncoding_flanks.constraint_counts.bed")
    shell:
        """
        bedtools intersect -sorted -wa -c -a {input.nc_flanks} -b {input.constraint} > {output}
        """
    
rule count_nonconstraint_pos:
    input:
        nc_flanks="resources/{chrom}.noncoding_flanks.bed",
        nonconstraint="resources/nonconstraint_scores/{chrom}.nonconstraint_positions.bed.gz"
    output:
        temp("resources/{chrom}.noncoding_flanks.nonconstraint_counts.bed")
    shell:
        """
        bedtools intersect -sorted -wa -c -a {input.nc_flanks} -b {input.nonconstraint} > {output}
        """

rule merge_phylop_counts:
    input:
        gene_set = config["gene_set"], 
        nccps=expand(["resources/{chrom}.noncoding_flanks.constraint_counts.bed"], chrom = chromosomes),
        ncncps=expand(["resources/{chrom}.noncoding_flanks.nonconstraint_counts.bed"], chrom = chromosomes)
    output:
        'resources/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+\
        '.phylop' + config['phyloP_threshold'] + '.phylop_counts.bed'
    shell:
        """
        workflow/scripts/merge_nc_phylop_positions.py \
        --gene_set {input.gene_set} --nccp '{input.nccps}' \
        --ncncp '{input.ncncps}' --output {output}
        """

rule phylop_counts:
    input: 
        'resources/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+\
        '.phylop' + config['phyloP_threshold'] + '.phylop_counts.bed'