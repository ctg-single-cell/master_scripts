import os
configfile: "config_fuma_ct.json"

scrnas = ["90_Siletti_Pons-PnRF_Human_2022", "91_Siletti_Pons-PnEN_Human_2022", "92_Siletti_Pons-PN_Human_2022"]

rule all:
    input:
        expand(os.path.join(config["out_dir"], "{trait}/{scrna}_magma.gsa.out.rmComment"), trait=config["traits"], scrna=scrnas),
        expand(os.path.join(config["out_dir"], "{trait}/magma_celltype_step1_significant_within.txt"), trait=config["traits"])


rule magma_step1:
    input:
        gene_results = os.path.join(config["gwas_dir"], config["gene_region"] + "_" + config["gene_name"],  "{trait}_" + config["gene_region"] + "_" + config["gene_name"] + "_magma.genes.raw"),
        gene_covar = os.path.join(config["scrna_dir"], "{scrna}/means_cell_log_counts_pM.tsv")
    output:
        os.path.join(config["out_dir"], "{trait}/{scrna}_magma.gsa.out.rmComment")
    params:
        magma_executable = config["magma_executable"],
        out_prefix = config["out_dir"] + "/{trait}/{scrna}_magma"
    shell:
        """
        {params.magma_executable} --gene-results {input.gene_results} --gene-covar {input.gene_covar} --model condition-hide=Average direction=greater --out {params.out_prefix};
        grep -v "#" {params.out_prefix}.gsa.out > {params.out_prefix}.gsa.out.rmComment
        """


rule correct_pval_magma_step1:
    input:
        expand(os.path.join(config["out_dir"], "{trait}/{scrna}_magma.gsa.out.rmComment"), scrna=scrnas, trait=config["traits"])
    output:
        os.path.join(config["out_dir"], "{trait}/magma_celltype_step1_significant_within.txt")
    params:
        script = confgig["calc_adj_p_script"],
        dir = os.path.join(config["out_dir"], "{trait}/")
    shell:
        """
        Rscript {params.script} --magma_outdir {params.dir}
        """