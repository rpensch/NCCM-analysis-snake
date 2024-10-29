# Check the vcf files, separate spims if needed
if config['spim'].lower() == 'combined':

    rule extract_spms:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['spim']
        output:
            temp("resources/input/{sample}.spm.vcf.gz"),
        shadow: "full"
        group: "separate_spims"
        shell:
            """
            vcftools --gzvcf {input} --remove-indels \
            --recode --recode-INFO-all --stdout |
            gzip -c > {output}
            """

    rule extract_sims:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['spim']
        output:
            temp("resources/input/{sample}.sim.vcf.gz")
        shadow: "full"
        group: "separate_spims"
        shell:
            """
            vcftools --gzvcf {input} --keep-only-indels \
            --recode --recode-INFO-all --stdout |
            gzip -c > {output}
            """

elif config['spim'].lower() == 'separated':

    rule link_vcfs:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample][wildcards.spim]
        output:
            "resources/input/{sample}.{spim}.vcf.gz"
        shell:
            "mkdir -p resources/input ; ln -s {input} {output}"