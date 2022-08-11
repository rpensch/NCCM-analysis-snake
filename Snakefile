import pandas as pd

configfile: "config/config.yaml"

vcfs = pd.read_table(config["vcfs"]).set_index("sample", drop=False)
samples = list(vcfs["sample"])
spims = ['spm','sim']

if 'spim' in vcfs.columns.tolist():

    spim_type = 'combined'

    rule extract_spms:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample]['spim']
        output:
            temp("input/{sample}.spm.vcf.gz"),
        shadow: "full"
        envmodules:
            "bioinfo-tools",
            "vcftools/0.1.12"
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
            temp("input/{sample}.sim.vcf.gz")
        shadow: "full"
        envmodules:
            "bioinfo-tools",
            "vcftools/0.1.12"
        group: "separate_spims"
        shell:
            """
            vcftools --gzvcf {input} --keep-only-indels \
            --recode --recode-INFO-all --stdout |
            gzip -c > {output}
            """

elif 'spm' in vcfs.columns.tolist() and 'sim' in vcfs.columns.tolist():

    spim_type = 'separated'

    rule link_vcfs:
        input:
            lambda wildcards: vcfs.loc[wildcards.sample][wildcards.spim]
        output:
            temp("input/{sample}.{spim}.vcf.gz")
        shell:
            "ln -s {input} {output}"

rule snpEff:
    input:
        vcf = "input/{sample}.{spim}.vcf.gz"
    output: 
        "snpEff/{sample}.{spim}.ann.vcf.gz"
    params:
        config["genome"]
    envmodules:
        "bioinfo-tools",
        "snpEff/4.3t"
    shell:
        """
        snpEff ann -noStats {params} {input.vcf} |
        gzip -c > {output}
        """

rule vcf2bed:
    input:
        expand("input/{{sample}}.{spim}.vcf.gz", spim = spims)
    output:
        expand("input_beds/{{sample}}.{spim}.sorted.bed", spim = spims)
    envmodules:
        "bioinfo-tools",
        "BEDOPS/2.4.39"
    shell:
        "scripts/vcf2bed.sh {wildcards.sample} {input} {output}"

rule phyloP:
    input:
        bed = "input_beds/{sample}.{spim}.sorted.bed",
        phyloP = config["phyloP"]
    output:
        "phyloP/{sample}.{spim}.phyloP.bed"
    shadow: "full"
    shell:
        """
        mkdir -p phyloP
        scripts/phylop_annotation_bw.sh {input.bed} {output} {input.phyloP}
        """

if spim_type == 'separated':

    rule annotation_summary:
        input:
            spm_vcf = vcfs['spm'].tolist(),
            sim_vcf = vcfs['sim'].tolist(),
            spm_snpEff = expand(["snpEff/{sample}.spm.ann.vcf.gz"], sample = samples),
            sim_snpEff = expand(["snpEff/{sample}.sim.ann.vcf.gz"], sample = samples),
            spm_phyloP = expand(["phyloP/{sample}.spm.phyloP.bed"], sample = samples),
            sim_phyloP = expand(["phyloP/{sample}.sim.phyloP.bed"], sample = samples)
        output:
            config["run_name"]+".annotation_summary.tsv"
        envmodules:
            "python3/3.9.5"
        shell:
            """
            scripts/annotation_summary.py --samples '{samples}' \
            --spm_vcfs '{input.spm_vcf}' --sim_vcfs '{input.sim_vcf}' \
            --spm_anns '{input.spm_snpEff}' --sim_anns '{input.sim_snpEff}' \
            --spm_phylops '{input.spm_phyloP}' --sim_phylops '{input.sim_phyloP}' \
            --out {output}
            """


elif spim_type == 'combined':

    rule annotation_summary:
        input:
            spim_vcf = vcfs['spim'].tolist(),
            spm_snpEff = expand("snpEff/{sample}.spm.ann.vcf.gz", sample = samples),
            sim_snpEff = expand("snpEff/{sample}.sim.ann.vcf.gz", sample = samples),
            spm_phyloP = expand("phyloP/{sample}.spm.phyloP.bed", sample = samples),
            sim_phyloP = expand("phyloP/{sample}.sim.phyloP.bed", sample = samples)
        output:
            config["run_name"]+".annotation_summary.tsv"
        envmodules:
            "python3/3.9.5"
        shell:
            """
            scripts/annotation_summary.py --samples '{samples}'' \
            --spim_vcfs '{input.spim_vcf}' \
            --spm_anns '{input.spm_snpEff}' --sim_anns '{input.sim_snpEff}' \
            --spm_phylops '{input.spm_phyloP}' --sim_phylops '{input.sim_phyloP}' \
            --out {output}
            """

rule annotate:
    input:
        config['run_name']+".annotation_summary.tsv"

rule combine_data: 
    input:
        snpEff = "snpEff/{sample}.{spim}.ann.vcf.gz",
        phyloP = "phyloP/{sample}.{spim}.phyloP.bed"
    output:
        temp("composite_matrix/{sample}.{spim}.composite_matrix.tsv.gz")
    envmodules:
            "python3/3.9.5"
    shell:
        """
        scripts/composite_matrix.py --sample {wildcards.sample} \
        --ann {input.snpEff} --phylop {input.phyloP} \
        --out {output}
        """

rule combine_all_samples:
    input:
        expand("composite_matrix/{sample}.{spim}.composite_matrix.tsv.gz", sample = samples, spim = spims)
    output:
        "composite_matrix/"+config["run_name"]+".composite_matrix.tsv.gz"
    shell:
        """
        for matrix in {input} ; do 
            if [ -f {output} ] ; then
                zcat "$matrix" | sed 1d | gzip -c >> {output}
            else
                zcat "$matrix" | gzip -c > {output}
            fi
        done
        """

rule composite_matrix:
    input:
        annotaiton = config['run_name']+".annotation_summary.tsv",
        composite_matrix = "composite_matrix/"+config["run_name"]+".composite_matrix.tsv.gz"

rule all:
    input:
        "phyloP/T-BlueSkye_vs_N-BlueSkye.spm.phyloP.bed",
        "phyloP/T-BlueSkye_vs_N-BlueSkye.sim.phyloP.bed"