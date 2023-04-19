import os
configfile: "config.json"


rule all:
    input:
        expand(os.path.join(config["out_dir"], "{trait}/{scrna}_magma.gsa.out.rmComment"), trait=config["traits"], scrna=config["scrnas"]),
        expand(os.path.join(config["out_dir"], "{trait}/magma_celltype_step1_significant_within.txt"), trait=config["traits"])


rule magma_step1:
    input:
        gene_results = os.path.join(config["data_dir"], "gene_results", config["gene_region"] + "_" + config["gene_name"],  "{trait}_" + config["gene_region"] + "_" + config["gene_name"] + "_magma.genes.raw"),
        gene_covar = os.path.join(config["data_dir"], "gene_covar", "{scrna}/means_cell_log_counts_pM_convert.tsv")
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
        expand(os.path.join(config["out_dir"], "{trait}/{scrna}_magma.gsa.out.rmComment"), scrna=config["scrnas"], trait=config["traits"])
    output:
        os.path.join(config["out_dir"], "{trait}/magma_celltype_step1_significant_within.txt")
    params:
        script = "calc_adj_pval_magma_step1.R",
        dir = os.path.join(config["out_dir"], "{trait}/")
    shell:
        """
        Rscript {params.script} --magma_outdir {params.dir}
        """