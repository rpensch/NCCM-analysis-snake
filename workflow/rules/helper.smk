rule split_coding_bed:
    input:
        config['coding_bed']
    output:
        temp("results/trinuc/coding.{chrom}.bed")
    shell:
        """
        awk -v chr='{wildcards.chrom}' -v OFS='\t' '$1==chr' {input} > {output}
        """

rule split_gene_flanks:
    input:
        gene_set = config["gene_set"]
    output:
        temp("results/trinuc/gene_set.{chrom}.bed")
    shell:
        """
        awk -v chr='{wildcards.chrom}' -v OFS='\t' '$1==chr' {input} > {output}
        """

