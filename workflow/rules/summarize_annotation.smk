# Summarize the annotation results
if config['spim'].lower() == 'separated':

    rule annotation_summary:
        input:
            spm_vcf = vcfs['spm'].tolist(),
            sim_vcf = vcfs['sim'].tolist(),
            spm_snpEff = expand(["results/annotation/snpEff/{sample}.spm.ann.vcf.gz"], sample = samples),
            sim_snpEff = expand(["results/annotation/snpEff/{sample}.sim.ann.vcf.gz"], sample = samples),
            spm_phyloP = expand(["results/annotation/phyloP/{sample}.spm.phyloP.bed"], sample = samples),
            sim_phyloP = expand(["results/annotation/phyloP/{sample}.sim.phyloP.bed"], sample = samples)
        output:
            "results/annotation/" + config["run_name"] + ".annotation_summary.tsv"
        params:
            samples_list=temp("list_samples.txt"),
            spm_vcf_list=temp("list_spm_vcf_files.txt"),
            sim_vcf_list=temp("list_sim_vcf_files.txt"),
            spm_snpeff_list=temp("list_spm_snpeff_files.txt"),
            sim_snpeff_list=temp("list_sim_snpeff_files.txt"),
            spm_phylop_list=temp("list_spm_phylop_files.txt"),
            sim_phylop_list=temp("list_sim_phylop_files.txt")
        shell:
            """
            echo {samples} | tr " " "\\n" > {params.samples_list}
            echo {input.spm_vcf} | tr " " "\\n" > {params.spm_vcf_list}
            echo {input.sim_vcf} | tr " " "\\n" > {params.sim_vcf_list}
            echo {input.spm_snpEff} | tr " " "\\n" > {params.spm_snpeff_list}
            echo {input.sim_snpEff} | tr " " "\\n" > {params.sim_snpeff_list}
            echo {input.spm_phyloP} | tr " " "\\n" > {params.spm_phylop_list}
            echo {input.sim_phyloP} | tr " " "\\n" > {params.sim_phylop_list}

            workflow/scripts/annotation_summary.py --samples {params.samples_list} \
            --spm_vcfs {params.spm_vcf_list} --sim_vcfs {params.sim_vcf_list} \
            --spm_anns {params.spm_snpeff_list} --sim_anns {params.sim_snpeff_list} \
            --spm_phylops {params.spm_phylop_list} --sim_phylops {params.sim_phylop_list} \
            --out {output}
            """

elif config['spim'].lower() == 'combined':

    rule annotation_summary:
        input:
            spim_vcf = vcfs['spim'].tolist(),
            spm_snpEff = expand(["results/annotation/snpEff/{sample}.spm.ann.vcf.gz"], sample = samples),
            sim_snpEff = expand(["results/annotation/snpEff/{sample}.sim.ann.vcf.gz"], sample = samples),
            spm_phyloP = expand(["results/annotation/phyloP/{sample}.spm.phyloP.bed"], sample = samples),
            sim_phyloP = expand(["results/annotation/phyloP/{sample}.sim.phyloP.bed"], sample = samples)
        output:
            "results/annotation/" + config["run_name"] + ".annotation_summary.tsv"
        params:
            samples_list=temp("list_samples.txt"),
            vcf_list=temp("list_spm_vcf_files.txt"),
            spm_snpeff_list=temp("list_spm_snpeff_files.txt"),
            sim_snpeff_list=temp("list_sim_snpeff_files.txt"),
            spm_phylop_list=temp("list_spm_phylop_files.txt"),
            sim_phylop_list=temp("list_sim_phylop_files.txt")
        shell:
            """
            echo {samples} | tr " " "\\n" > {params.samples_list}
            echo {input.spim_vcf} | tr " " "\\n" > {params.vcf_list}
            echo {input.spm_snpEff} | tr " " "\\n" > {params.spm_snpeff_list}
            echo {input.sim_snpEff} | tr " " "\\n" > {params.sim_snpeff_list}
            echo {input.spm_phyloP} | tr " " "\\n" > {params.spm_phylop_list}
            echo {input.sim_phyloP} | tr " " "\\n" > {params.sim_phylop_list}

            workflow/scripts/annotation_summary.py --samples {params.samples_list} \
            --spim_vcfs {params.vcf_list} \
            --spm_anns {params.spm_snpeff_list} --sim_anns {params.sim_snpeff_list} \
            --spm_phylops {params.spm_phylop_list} --sim_phylops {params.sim_phylop_list} \
            --out {output}
            """

rule annotate:
    input:
        "results/annotation/" + config['run_name'] + ".annotation_summary.tsv"