import os

rule all:
    input:
        expand(os.path.join(config["magma_dir"], "spec_{spec_method}_linear/{scrna}/{trait}_spec_{spec_method}_linear_magma.gsa.out.rmComment"), trait=config["traits"], scrna=config["scrnas"], spec_method=config["spec_methods"]),
        expand(os.path.join(config["magma_dir"], "spec_{spec_method}_top10/{scrna}/{trait}_spec_{spec_method}_top10_magma.gsa.out.rmComment"), trait=config["traits"], scrna=config["scrnas"], spec_method=config["spec_methods"])

rule magma_spec_linear:
    input:
        gene_results = os.path.join(config["gwas_dir"], "{trait}/{trait}_magma.genes.raw"),
        gene_covar = os.path.join(config["magma_dir"], "spec_{spec_method}_linear/{scrna}/gene_covar/spec_{spec_method}_linear.txt")
    output:
        os.path.join(config["magma_dir"], "spec_{spec_method}_linear/{scrna}/out/{trait}_spec_{spec_method}_linear_magma.gsa.out.rmComment")
    params:
        magma_executable = config["magma_executable"],
        out_prefix = config["magma_dir"] + "spec_{spec_method}_linear/{scrna}/out/{trait}_spec_{spec_method}_linear_magma"
    shell:
        """
        {params.magma_executable} --gene-results {input.gene_results} --gene-covar {input.gene_covar} --model direction=greater --out {params.out_prefix};
        grep -v "#" {params.out_prefix}.gsa.out > {params.out_prefix}.gsa.out.rmComment
        """

rule magma_spec_top10:
    input:
        gene_results = os.path.join(config["gwas_dir"], "{trait}/{trait}_magma.genes.raw"),
        gene_covar = os.path.join(config["magma_dir"], "spec_{spec_method}_top10/{scrna}/gene_covar/spec_{spec_method}_top10.txt")
    output:
        os.path.join(config["magma_dir"], "spec_{spec_method}_top10/{scrna}/out/{trait}_spec_{spec_method}_top10_magma.gsa.out.rmComment")
    params:
        magma_executable = config["magma_executable"],
        out_prefix = config["magma_dir"] + "spec_{spec_method}_top10/{scrna}/out/{trait}_spec_{spec_method}_top10_magma"
    shell:
        """
        {params.magma_executable} --gene-results {input.gene_results} --set-annot {input.gene_covar} --out {params.out_prefix};
        grep -v "#" {params.out_prefix}.gsa.out > {params.out_prefix}.gsa.out.rmComment
        """