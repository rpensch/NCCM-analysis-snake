# Combine all annotation and phyloP scores in a matrix
rule combine_data: 
    input:
        snpEff = "results/annotation/snpEff/{sample}.{spim}.ann.vcf.gz",
        phyloP = "results/annotation/phyloP/{sample}.{spim}.phyloP.bed"
    output:
        temp("results/composite_matrix/{sample}.{spim}.composite_matrix.tsv.gz")
    shell:
        """
        workflow/scripts/composite_matrix.py --sample {wildcards.sample} \
        --ann {input.snpEff} --phylop {input.phyloP} \
        --out {output}
        """

rule combine_all_samples:
    input:
        expand("results/composite_matrix/{sample}.{spim}.composite_matrix.tsv.gz", sample = samples, spim = spims)
    output:
        "results/composite_matrix/" + config["run_name"] + ".composite_matrix.tsv.gz"
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
        annotation = "results/annotation/" + config['run_name'] + ".annotation_summary.tsv",
        composite_matrix = "results/composite_matrix/" + config["run_name"] + ".composite_matrix.tsv.gz"
