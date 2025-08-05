import os

# Get the trinucleotide context of potential NCCM sites as NCC positions per gene
rule split_coding_bed:
    input:
        config['coding_bed']
    output:
        temp("results/trinuc/coding.{chrom}.bed")
    shell:
        """
        awk -v chr='{wildcards.chrom}' -v OFS='\t' '$1==chr' {input} > {output}
        """

rule split_gene_flanks:
    input:
        gene_set = config["gene_set"]
    output:
        temp("results/trinuc/gene_set.{chrom}.bed")
    shell:
        """
        awk -v chr='{wildcards.chrom}' -v OFS='\t' '$1==chr' {input} > {output}
        """

rule count_gene_ncc_trinucs:
    input: 
        coding_bed = "results/trinuc/coding.{chrom}.bed",
        gene_flanks = "results/trinuc/gene_set.{chrom}.bed",
        phyloP = config['phyloP'] + "/" + "{chrom}.bed.gz",
        ref_fasta = config['ref_fasta']
    output:
        temp("results/trinuc/gene_flanks.ncc_positions_trinuc_counts.{chrom}.tsv")
    params:
        phyloP_threshold = config["phyloP_threshold"]
    shell:
        """
        zcat {input.phyloP} | awk -v phylop="{params.phyloP_threshold}" '$5>=phylop' | 
        bedtools intersect -sorted -wa -wb -a stdin -b {input.gene_flanks} | 
        bedtools subtract -a stdin -b {input.coding_bed} | 
        awk -v OFS='\t' '{{print $1,$2-1,$3+1,$9}}' |
        bedtools getfasta -fi {input.ref_fasta} -bed stdin -bedOut |
        cut -f 4,5 | sort | uniq -c | awk -v OFS='\t' '{{print $2,$3,$1}}' \
        > {output}
        """

rule total_gene_ncc_trinucs:
    input:
        gene_set = config["gene_set"],
        counts = expand(["results/trinuc/gene_flanks.ncc_positions_trinuc_counts.{chrom}.tsv"], chrom = chromosomes)
    output:
        'results/trinuc/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+".ncc_positions_trinuc_counts.tsv"
    shell:
        """
        workflow/scripts/combine_gene_ncc_trinucs.py --counts '{input.counts}' --gene_flanks {input.gene_set} --output {output}
        """

rule gene_ncc_trinucs:
    input: 'results/trinuc/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+".ncc_positions_trinuc_counts.tsv"
