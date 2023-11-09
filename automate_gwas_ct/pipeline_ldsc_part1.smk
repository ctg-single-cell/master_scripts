import os

# IMPORTANT: in this snakemake, the list of cell types from the scRNAseq data is pre-generated. In Rachel's original pipeline, dynamic was used to get around not knowing aprior the cell type annotation per dataset. While this is a clever solution, I tend to avoid this because it is difficult for me to troubleshoot if anything goes wrong
# Therefore, I prefer to create a list of cell type apriori. The disadvantage of this is that the cell type annotation is going to differ between datasets so it is a bit more manual work to ensure that the cell type annotation list is correct with each dataset. However, I don't think that this is going to be a huge problem because we can pass the file that contains the cell types on each line as a config
# Generate 2 lists to list the cell type annotation
# list 1: all_cts: all cell types with All
# list 2: cts_wo_all: all cell types without All
# example python scripts:

# with open(config["celltype_list"], "r") as f: #assuming that this file contains the cell types on each line
#    ct_wo_all = [line.rstrip("\n") for line in f]

# all_cts = ct_wo_all + ["All"]

with open(config["cell_type_list"], "r") as f:
   ct_wo_all = [line.rstrip("\n") for line in f]
all_cts = ct_wo_all + ["All"]

# all_cts = ["Committed_oligodendrocyte_precursor", "Oligodendrocyte", "Oligodendrocyte_precursor", "Splatter", "Upper_rhombic_lip", "Cerebellar_inhibitory", "Miscellaneous", "Astrocyte", "Bergmann_glia", "Ependymal", "Vascular", "Fibroblast", "Microglia", "All"]#testing for 41_Siletti_Cerebellum.CBV_Human_2022 only
# cts_wo_all = ["Committed_oligodendrocyte_precursor", "Oligodendrocyte", "Oligodendrocyte_precursor", "Splatter", "Upper_rhombic_lip", "Cerebellar_inhibitory", "Miscellaneous", "Astrocyte", "Bergmann_glia", "Ependymal", "Vascular", "Fibroblast", "Microglia"]

rule all: 
   input:
      expand(os.path.join(config["ldsc_dir"], "spec_{method}/{scrna}/LD_annotations/{ct}.{chr}.annot.gz"), method=config["spec_methods"], scrna=config["scrnas"], ct=all_cts, chr=config["chrs"]),
      expand(os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/LD_annotations/{ct}.{chr}.l2.ldscore.gz"), method=config["spec_methods"], scrna=config["scrnas"], ct=all_cts, chr=config["chrs"]),
      expand(os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats_{ct}.results"), method=config["spec_methods"], scrna=config["scrnas"], ct=ct_wo_all, trait=config["traits"])


rule genelist_to_annot:
   output:
      os.path.join(config["ldsc_dir"], "spec_{method}/{scrna}/LD_annotations/{ct}.{chr}.annot.gz")
   params:
      script=os.path.join(config["ldsc_path"],"make_annot.py"),
      script_params=" --gene-set-file " + os.path.join(config["ldsc_dir"], "spec_{method}/{scrna}/genelists/{ct}.txt") + " --gene-coord-file " + os.path.join(config["reference_data"],"all_ensmbl_genes_coord.txt") + " --windowsize 100000 --bimfile " + os.path.join(config["reference_data"],"1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim") + " --annot " + os.path.join(config["ldsc_dir"], "spec_{method}/{scrna}/LD_annotations/{ct}.{chr}.annot.gz")
   resources:
      mem_mb=1000
   conda: 
      "/gpfs/home6/tphung/software/ldsc/environment.yml"
   shell:
      """
      python {params.script} {params.script_params}
      """

rule annot_to_ldscore:
   input:
      file=os.path.join(config["ldsc_dir"], "spec_{method}/{scrna}/LD_annotations/{ct}.{chr}.annot.gz")
   output:
      os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/LD_annotations/{ct}.{chr}.l2.ldscore.gz"),
   params:
      script=os.path.join(config["ldsc_path"],"ldsc.py"),
      script_params=" --bfile " + os.path.join(config["reference_data"],"1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}") + " --l2 --ld-wind-cm 1 --yes-really --thin-annot --annot " + os.path.join(config["ldsc_dir"], "spec_{method}/{scrna}/LD_annotations/{ct}.{chr}.annot.gz") + " --print-snps " + os.path.join(config["reference_data"],"baseline_v2.1_snplists/snplist_chrom{chr}") + " --out " + os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/LD_annotations/{ct}.{chr}")
   conda:
      "/gpfs/home6/tphung/software/ldsc/environment.yml"
   resources:
      mem_mb=3000
   shell:
      """
      python {params.script} {params.script_params}
      """

rule ldsc_regression:
   output: 
      os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats_{ct_wo_all}.results")
   params:
      script=os.path.join(config["ldsc_path"],"ldsc.py"),
      script_params=" --h2 " + os.path.join(config["gwas_dir"],"{trait}/{trait}.sumstats.gz") + " --frqfile-chr " + os.path.join(config["reference_data"],"1000G_Phase3_frq/1000G.EUR.QC.") + " --ref-ld-chr " + os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/LD_annotations/{ct_wo_all}.,") + os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/LD_annotations/All.,") + os.path.join(config["reference_data"],"baselineLD_v2.1_annots/baselineLD.")  + " --w-ld-chr " + os.path.join(config["reference_data"],"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.") + " --overlap-annot --print-coefficients --print-delete-vals --out " + os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats_{ct_wo_all}")
   resources:
      mem_mb=3000
   conda: 
      "/gpfs/home6/tphung/software/ldsc/environment.yml"
   shell:
      """
      python {params.script} {params.script_params}
      """


