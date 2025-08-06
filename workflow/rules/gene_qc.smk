rule gene_qc_overlap:
    input:
        gene_set = config["gene_set"],
        qc_data = config["qc_bed"]
    output:
        'resources/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+ '.' + \
        '.'.join(os.path.basename(config["qc_bed"]).split('.')[:-1])+ '.perc_overlap.bed'
    shell:
        """
        workflow/scripts/gene_qc.sh {input.gene_set} {input.qc_data} {output}
        """

rule gene_qc:
    input:
        'resources/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+ '.' + \
        '.'.join(os.path.basename(config["qc_bed"]).split('.')[:-1])+ '.perc_overlap.bed'
