# This snakemake is to run the Rscript `postprocess_cellect_nonABA.R`
# Create on 2023-07-24
# Author: Tanya Phung (t.n.phung@vu.nl)
# Change log
# Date: 2023-08-03: editted for cellect-ldsc

import os

traits = ["sumstats_f1_N", "sumstats_f2_N", "sumstats_f3_N", "sumstats_f4_N"]
scrnas = []

scrnas_id_rm = ["31_Siletti_CerebralCortex-FuGt-TF_Human_2022", "54_Siletti_CerebralNuclei-CaB_Human_2022"] #added on 2023-08-03 because these 2 didn't work with cellect-ldsc
with open("/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/scrnaseq_data_master.csv", "r") as f:
    for line in f:
        if not line.startswith("Internal_ID"):
            internal_id = line.rstrip("\n").split(",")[0]
            if "Siletti" in internal_id:
                new_id = internal_id.replace(".", "-")
                if new_id not in scrnas_id_rm:
                    scrnas.append(new_id)

rule all:
    input:
        expand("/gpfs/home6/tphung/projects/alc_jeanne/out/ct_enrichment/cellect/CELLECT-LDSC/results/{trait}/{scrnaseq_id}_cellect_magma_significant_within.txt", trait=traits, scrnaseq_id=scrnas),
        expand("/gpfs/home6/tphung/projects/alc_jeanne/out/ct_enrichment/cellect/CELLECT-MAGMA/results/{trait}/{scrnaseq_id}_cellect_magma_significant_within.txt", trait=traits, scrnaseq_id=scrnas)

rule postprocess_cellect_magma:
    input:
        "/gpfs/home6/tphung/projects/alc_jeanne/out/ct_enrichment/cellect/CELLECT-MAGMA/results/prioritization.csv"
    output:
        "/gpfs/home6/tphung/projects/alc_jeanne/out/ct_enrichment/cellect/CELLECT-MAGMA/results/{trait}/{scrnaseq_id}_cellect_magma_significant_within.txt"
    params:
        trait = "{trait}",
        scrnaseq_id = "{scrnaseq_id}"
    shell:
        """
        Rscript postprocess_cellect_nonABA.R --input {input} --trait {params.trait} --scrnaseq_id {params.scrnaseq_id} --outfile {output}
        """

rule postprocess_cellect_ldsc:
    input:
        "/gpfs/home6/tphung/projects/alc_jeanne/out/ct_enrichment/cellect/CELLECT-LDSC/results/prioritization.csv"
    output:
        "/gpfs/home6/tphung/projects/alc_jeanne/out/ct_enrichment/cellect/CELLECT-LDSC/results/{trait}/{scrnaseq_id}_cellect_magma_significant_within.txt"
    params:
        trait = "{trait}",
        scrnaseq_id = "{scrnaseq_id}"
    shell:
        """
        Rscript postprocess_cellect_nonABA.R --input {input} --trait {params.trait} --scrnaseq_id {params.scrnaseq_id} --outfile {output}
        """
