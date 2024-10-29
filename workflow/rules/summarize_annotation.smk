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
        shell:
            """
            workflow/scripts/annotation_summary.py --samples '{samples}' \
            --spm_vcfs '{input.spm_vcf}' --sim_vcfs '{input.sim_vcf}' \
            --spm_anns '{input.spm_snpEff}' --sim_anns '{input.sim_snpEff}' \
            --spm_phylops '{input.spm_phyloP}' --sim_phylops '{input.sim_phyloP}' \
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
        shell:
            """
            workflow/scripts/annotation_summary.py --samples '{samples}' \
            --spim_vcfs '{input.spim_vcf}' \
            --spm_anns '{input.spm_snpEff}' --sim_anns '{input.sim_snpEff}' \
            --spm_phylops '{input.spm_phyloP}' --sim_phylops '{input.sim_phyloP}' \
            --out {output}
            """

rule annotate:
    input:
        "results/annotation/" + config['run_name'] + ".annotation_summary.tsv"