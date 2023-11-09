import os

rule all:
    input:
        expand(os.path.join(config["magma_dir"], "mean_{mean_method}_linear/{scrna}/out/{trait}_mean_{mean_method}_linear_magma.gsa.out.rmComment"), trait=config["traits"], scrna=config["scrnas"], mean_method=config["mean_methods"])

rule magma_mean_linear:
    input:
        gene_results = os.path.join(config["gwas_dir"], "{trait}/{trait}_magma.genes.raw"),
        gene_covar = os.path.join(config["magma_dir"], "mean_{mean_method}_linear/{scrna}/gene_covar/mean_{mean_method}_linear.tsv")
    output:
        os.path.join(config["magma_dir"], "mean_{mean_method}_linear/{scrna}/out/{trait}_mean_{mean_method}_linear_magma.gsa.out.rmComment")
    params:
        magma_executable = config["magma_executable"],
        out_prefix = config["magma_dir"] + "/mean_{mean_method}_linear/{scrna}/out/{trait}_mean_{mean_method}_linear_magma"
    shell:
        """
        {params.magma_executable} --gene-results {input.gene_results} --gene-covar {input.gene_covar} --model condition-hide=Average direction=greater --out {params.out_prefix};
        grep -v "#" {params.out_prefix}.gsa.out > {params.out_prefix}.gsa.out.rmComment
        """