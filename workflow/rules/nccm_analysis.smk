# Run the NCCM analysis
rule nccm_analysis:
    input:
        matrix = "results/composite_matrix/" + config["run_name"] + ".composite_matrix.tsv.gz",
        gene_set = config["gene_set"]
    output:
        "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"] + ".scan.tsv",
        "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"] + ".top_nccm_genes_spim.tsv",
        "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"] + ".top_nccm_genes.list"
    params:
        phyloP_threshold = config["phyloP_threshold"],
        prefix = "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"]
    threads: config["threads"]
    shell:
        """
        mkdir -p results/nccms ; \
        workflow/scripts/1.1.NCCM_analysis.sh {input.matrix} {input.gene_set} \
        {params.prefix} {params.phyloP_threshold} 10 7 {threads}
        """
