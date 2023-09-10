import os

scrnas = []

with open("/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/scrnaseq_data_master.csv", "r") as f:
    for line in f:
        if not line.startswith("Internal_ID"):
            internal_id = line.rstrip("\n").split(",")[0]
            if "Siletti" in internal_id:
                scrnas.append(internal_id)

# # run on ABA level 1
# scrnas = ["1_Allen_MCA_Human_2019", "2_Allen_M1_Human_2020", "3_Allen_MTG_Human_2022", "4_Allen_Cortex.Hippocampus_Mouse_2019"]
# types = ["linear", "top10"]

rule all:
    input:
        expand(os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment"), trait=config["traits"], h5ad_id=scrnas),
        expand(os.path.join(config["output_dir"], "{trait}/top10/magma_celltype_step1_significant_within.txt"), trait=config["traits"])
        # expand(os.path.join(config["output_dir"], "{trait}", "{type}/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment"), trait=config["traits"], type=types, h5ad_id=scrnas),
        # expand(os.path.join(config["output_dir"], "{trait}/linear/magma_celltype_step1_significant_within.txt"), trait=config["traits"])
        # expand(os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out"), trait=config["traits"], h5ad_id=scrnas),
        # expand(os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment"), trait=config["traits"], h5ad_id=scrnas),
        # expand(os.path.join(config["output_dir"], "{trait}/top10/magma_celltype_step1_significant_within.txt"), trait=config["traits"])


## The following 2 rules `process_scrna` and `run_ewce` only need to be run once per scRNAseq dataset
# rule process_scrna:
#     input:
#         h5ad = os.path.join(config["preprocessed_h5ad_dir"], "{h5ad_id}/{h5ad_id}.h5ad")
#     output:
#         os.path.join(config["processed_h5ad_dir"], "{h5ad_id}_processed_magma_celltyping.h5ad")
#     shell:
#         """
#         python process_scrna_magma_celltyping.py --h5ad {input.h5ad} --out_h5ad {output}
#         """

# rule run_ewce:
#     input:
#         h5ad = os.path.join(config["processed_h5ad_dir"], "{h5ad_id}_processed_magma_celltyping.h5ad")
#     output:
#         os.path.join(config["processed_h5ad_dir"], "{h5ad_id}_" + "siletti.rda")
#     params:
#         outdir = config["processed_h5ad_dir"],
#         file_prefix = "{h5ad_id}",
#         groupName = "siletti"
#     shell:
#         """
#         Rscript process_scrna_ewce.R --h5ad {input.h5ad} --outdir {params.outdir} --file_prefix {params.file_prefix} --groupName {params.groupName}
#         """
## End 2 rules that only need to be run once per scRNAseq dataset

rule generate_genecovar_top10:
    input:
        ctd = os.path.join(config["processed_h5ad_dir"], "{h5ad_id}_" + "protein_coding_entrez.rda"),
        magma_gene_out = os.path.join(config["magma_gene_dir"], config["groupName"],  "{trait}_" + config["groupName"] + "_magma.genes.out")
    output:
        level_1 = os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1")
    params:
        out_basename = os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}")
    shell:
        """
        Rscript create_geneCovar_top10.R --ctd {input.ctd} --magma_gene_out {input.magma_gene_out} --out_basename {params.out_basename}
        """

rule run_magma_top10:
    input:
        magma_gene_raw = os.path.join(config["magma_gene_dir"], config["groupName"], "{trait}_" + config["groupName"] + "_magma.genes.raw"),
        level_1 = os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1")
    output:
        os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out")
    params:
        magma_path = config["magma_path"],
        out_basename_level_1 = os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping")
    shell:
        """
        {params.magma_path} --gene-results {input.magma_gene_raw} --set-annot {input.level_1} --out {params.out_basename_level_1}
        """

rule remove_comments:
    input:
        os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out")
    output:
        os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment")
    shell:
        """
        grep -v "#" {input} > {output}
        """

rule correct_pval:
    input:
        top10 = expand(os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment"), h5ad_id=scrnas, trait=config["traits"])
    output:
        top10 = os.path.join(config["output_dir"], "{trait}/top10/magma_celltype_step1_significant_within.txt")
    params:
        script = "/gpfs/home6/tphung/software/master_scripts/FUMA_Celltype_cmd/calc_adj_pval_magma_step1.R",
        top10 = os.path.join(config["output_dir"], "{trait}/top10/")
    shell:
        """
        Rscript {params.script} --magma_outdir {params.top10}
        """

# rule generate_genecovar:
#     input:
#         ctd = os.path.join(config["processed_h5ad_dir"], "{h5ad_id}_" + "protein_coding_entrez.rda"),
#         magma_gene_out = os.path.join(config["magma_gene_dir"], config["groupName"], "{trait}_" + config["groupName"] + "_magma.genes.out")
#     output:
#         level_1 = os.path.join(config["output_dir"], "{trait}", "linear/tmp_{trait}_{h5ad_id}_level_1")
#     params:
#         out_basename = os.path.join(config["output_dir"], "{trait}", "linear/tmp_{trait}_{h5ad_id}")
#     shell:
#         """
#         Rscript create_geneCovar.R --ctd {input.ctd} --magma_gene_out {input.magma_gene_out} --out_basename {params.out_basename}
#         """
    
# rule fix_geneCovar:
#     input:
#         level_1 = os.path.join(config["output_dir"], "{trait}", "linear/tmp_{trait}_{h5ad_id}_level_1")
#     output:
#         level_1 = os.path.join(config["output_dir"], "{trait}", "linear/{trait}_{h5ad_id}_level_1")
#     params:
#     shell:
#         """
#         python fix_geneCovar.py --input {input.level_1} --output {output.level_1}
#         """

# rule run_magma:
#     input:
#         magma_gene_raw = os.path.join(config["magma_gene_dir"], config["groupName"], "{trait}_" + config["groupName"] + "_magma.genes.raw"),
#         level_1 = os.path.join(config["output_dir"], "{trait}", "linear/{trait}_{h5ad_id}_level_1")
#     output:
#         os.path.join(config["output_dir"], "{trait}", "linear/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out")
#     params:
#         magma_path = config["magma_path"],
#         out_basename_level_1 = os.path.join(config["output_dir"], "{trait}", "linear/{trait}_{h5ad_id}_level_1_magma_celltyping")
#     shell:
#         """
#         {params.magma_path} --gene-results {input.magma_gene_raw} --gene-covar {input.level_1} --model direction=pos --out {params.out_basename_level_1}
#         """

# rule remove_comments:
#     input:
#         os.path.join(config["output_dir"], "{trait}", "{type}/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out")
#     output:
#         os.path.join(config["output_dir"], "{trait}", "{type}/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment")
#     shell:
#         """
#         grep -v "#" {input} > {output}
#         """

# rule correct_pval:
#     input:
#         linear = expand(os.path.join(config["output_dir"], "{trait}", "linear/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment"), h5ad_id=scrnas, trait=config["traits"]),
#         top10 = expand(os.path.join(config["output_dir"], "{trait}", "top10/{trait}_{h5ad_id}_level_1_magma_celltyping.gsa.out.rmComment"), h5ad_id=scrnas, trait=config["traits"])
#     output:
#         linear = os.path.join(config["output_dir"], "{trait}/linear/magma_celltype_step1_significant_within.txt"),
#         top10 = os.path.join(config["output_dir"], "{trait}/top10/magma_celltype_step1_significant_within.txt")
#     params:
#         script = "/gpfs/home6/tphung/software/master_scripts/FUMA_Celltype_cmd/calc_adj_pval_magma_step1.R",
#         linear = os.path.join(config["output_dir"], "{trait}/linear/"),
#         top10 = os.path.join(config["output_dir"], "{trait}/top10/")
#     shell:
#         """
#         Rscript {params.script} --magma_outdir {params.linear};
#         Rscript {params.script} --magma_outdir {params.top10}
#         """
