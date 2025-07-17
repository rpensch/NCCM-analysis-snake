# Check the vcf files, separate spims if needed
if config['spim'].lower() == 'combined':

    rule extract_spms:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['spim']
        output:
            temp("resources/input/{sample}.spm.vcf")
        group: "separate_spims"
        shell:
            """
            zcat {input} | grep "^#" > {output}
            zcat {input} | grep -v "^#" | awk '(length($4)==1) && (length($5)==1)' >> {output}
            """

    rule extract_sims:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['spim']
        output:
            temp("resources/input/{sample}.sim.vcf")
        group: "separate_spims"
        shell:
            """
            zcat {input} | grep "^#" > {output}
            zcat {input} | grep -v "^#" | awk '(length($4)>1) || (length($5)>1)' >> {output}
            """

elif config['spim'].lower() == 'separated':

    rule merge_spim:
        input:
            spm=lambda wildcards: vcfs.loc[wildcards.sample]['spm'],
            sim=lambda wildcards: vcfs.loc[wildcards.sample]['sim']
        output:
            temp("resources/input/{sample}.spim.merged.tmp.vcf")
        shell:
            """
            zcat {input.spm} > {output}
            zcat {input.sim} | grep -v "^#" >> {output}
            """

    rule extract_spms:
        input:
            "resources/input/{sample}.spim.merged.tmp.vcf"
        output:
            temp("resources/input/{sample}.spm.vcf")
        shadow: "full"
        group: "separate_spims"
        shell:
            """
            cat {input} | grep "^#" > {output}
            cat {input} | grep -v "^#" | awk '(length($4)==1) && (length($5)==1)' >> {output}
            """

    rule extract_sims:
        input:
            "resources/input/{sample}.spim.merged.tmp.vcf"
        output:
            temp("resources/input/{sample}.sim.vcf")
        shadow: "full"
        group: "separate_spims"
        shell:
            """
            cat {input} | grep "^#" > {output}
            cat {input} | grep -v "^#" | awk '(length($4)>1) || (length($5)>1)' >> {output}
            """

rule compress:
    input:
        "resources/input/{sample}.{spim}.vcf"
    output:
        "resources/input/{sample}.{spim}.vcf.gz"
    shell:
        """
        gzip {input}
        """

rule preprocess:
    input:
        expand(["resources/input/{sample}.spm.vcf.gz"], sample = samples),
        expand(["resources/input/{sample}.sim.vcf.gz"], sample = samples)