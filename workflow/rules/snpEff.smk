# Annotate the vcf files wih snpEff
rule snpEff:
    input:
        vcf = "resources/input/{sample}.{spim}.vcf.gz"
    output: 
        "results/annotation/snpEff/{sample}.{spim}.ann.vcf.gz"
    params:
        config["genome"]
    shell:
        """
        snpEff ann -noStats {params} {input.vcf} |
        gzip -c > {output}
        """