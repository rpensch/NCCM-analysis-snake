rule all:
    input:
        "results/annotation/" + config['run_name'] + ".annotation_summary.tsv",
        "results/composite_matrix/" + config["run_name"] + ".composite_matrix.tsv.gz",
        "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"] + ".scan.tsv",
    params:
        prefix = "results/nccms/" + config["run_name"] + ".phylop-" + config["phyloP_threshold"]