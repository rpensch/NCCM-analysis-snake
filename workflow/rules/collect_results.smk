rule all:
    input:
        "results/annotation/" + config['run_name'] + ".annotation_summary.tsv",
        "results/composite_matrix/" + config["run_name"] + ".composite_matrix.tsv.gz",
        "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"] + ".scan.tsv",
        "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"] + ".top_nccm_genes_spim.tsv",
        "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"] + ".top_nccm_genes.list"
    params:
        prefix = "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"]