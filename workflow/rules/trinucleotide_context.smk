# Get the trinucleotide context of potential NCCM sites as NCC positions per gene
rule count_gene_ncc_trinucs_chr:
    input: 
        coding_bed = "results/resources/coding.{chrom}.bed",
        gene_flanks = "results/resources/gene_set.{chrom}.bed",
        phyloP = config['phyloP'] + "/" + "{chrom}.bed.gz",
        ref_fasta = config['ref_fasta']
    output:
        temp("results/resources/gene_flanks.ncc_positions_trinuc_counts.{chrom}.tsv")
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

rule combine_gene_ncc_trinucs:
    input:
        gene_set = config["gene_set"],
        counts = expand(["results/resources/gene_flanks.ncc_positions_trinuc_counts.{chrom}.tsv"], chrom = chromosomes)
    output:
        'results/resources/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+".ncc_positions_trinuc_counts.tsv"
    shell:
        """
        workflow/scripts/combine_gene_ncc_trinucs.py --counts '{input.counts}' --gene_flanks {input.gene_set} --output {output}
        """

rule gene_ncc_trinucs:
    input: 'results/resources/'+'.'.join(os.path.basename(config["gene_set"]).split('.')[:-1])+".ncc_positions_trinuc_counts.tsv"

