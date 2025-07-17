# Check the vcf files, separate spims if needed
if config['spim'].lower() == 'combined':

    rule extract_spms:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['spim']
        output:
            temp("resources/input/{sample}.spm.vcf.gz")
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

    rule bgzip_spm:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['spm'],
        output:
            spm=temp("resources/input/{sample}.spm.tmp.vcf.bgz"),
            spm_index=temp("resources/input/{sample}.spm.tmp.vcf.bgz.csi")
        shell:
            """
            bcftools sort {input}| bgzip -c > {output.spm}
            bcftools index -tbi {output.spm}
            """
    
    rule bgzip_sim:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['sim']
        output:
            sim=temp("resources/input/{sample}.sim.tmp.vcf.bgz"),
            sim_index=temp("resources/input/{sample}.sim.tmp.vcf.bgz.csi")
        shell:
            """
            bcftools sort {input}| bgzip -c > {output.sim}
            bcftools index -tbi {output.sim}
            """

    rule merge_spim:
        input:
            spm="resources/input/{sample}.spm.tmp.vcf.bgz",
            spm_index="resources/input/{sample}.spm.tmp.vcf.bgz.csi",
            sim="resources/input/{sample}.sim.tmp.vcf.bgz",
            sim_index="resources/input/{sample}.sim.tmp.vcf.bgz.csi"
        output:
            spim=temp("resources/input/{sample}.spim.merged.tmp.vcf.gz")
        shell:
            """
            bcftools concat {input.spm} {input.sim} -a -Oz -o {output.spim}
            """

    rule extract_spms:
        input:
            "resources/input/{sample}.spim.merged.tmp.vcf.gz"
        output:
            temp("resources/input/{sample}.spm.vcf.gz")
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
            "resources/input/{sample}.spim.merged.tmp.vcf.gz"
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

rule preprocess:
    input:
        expand(["resources/input/{sample}.spm.vcf.gz"], sample = samples),
        expand(["resources/input/{sample}.sim.vcf.gz"], sample = samples)