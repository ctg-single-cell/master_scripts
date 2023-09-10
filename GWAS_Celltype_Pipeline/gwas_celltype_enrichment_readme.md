# Overview
- This document details how to run gwas celltype enrichment 
- Author: Tanya Phung (t.n.phung@vu.nl)
- version 1.0 (date: 2023-09-07)
## Change log: 

# Preprocess the scRNAseq file

# fuma_ct
## Step 1: Create `magma.genes.raw` from GWAS summary statistics
- From GWAS summary statistics for each trait, create `magma.genes.raw`
    - Example: https://github.com/ctg-single-cell/master_scripts/blob/main/Preprocessing_GWAS/MAGMA/run_magma_gene.do
    - Please check on Snellius to see if your trait of interest has already been computed. Location on Snellius: `/gpfs/work4/0/vusr0455/sumstats/magma_gene/`
  
## Step 2: Compute mean from scRNAseq
### Quick start:
- Edit the snakemake file (on snellius: `/gpfs/work5/0/vusr0480/Processed_scRNA/code/fuma_ct/snakefile_compute_sumstat_magma.smk`). Specifically:
    - Line 2
    - Line 3 
        - If your scRNAseq file name is in this format: `1_Allen_MCA_Human_2019.h5ad`, then the id is `1_Allen_MCA_Human_2019`
    - Line 24, 25, and 26 to give the correct path to the script `compute_sumstat_magma.py` (this script can also be found on github)
- Dry run the snakemake file with the `-np` to make sure that the snakemake file is coded correctly
- Submit the snakemake file via `sbatch`
    - An example submit script is here: `/gpfs/work5/0/vusr0480/Processed_scRNA/code/fuma_ct/submit_smk.sh` but you need to edit it appropriately
### More detailed information
- To compute mean from the scRNAseq data, use the Python script https://github.com/ctg-single-cell/master_scripts/blob/main/Compute_FUMA_mean/compute_sumstat_magma.py
- A snakemake file is already prepared. Unfortunately this snakemake file contains many hard-coded items and one would need to edit it in order to run with new data (see above for what you need to do exactly)
  
## Step 3: Run MAGMA gene property
- Edit the config file `config_fuma_ct` (find this file in the folder fuma_ct)
    - On line 9, the Rscript is here: https://github.com/ctg-single-cell/master_scripts/blob/main/FUMA_Celltype_cmd/calc_adj_pval_magma_step1.R
- Edit the snakemake file `fuma_celltype_step1.smk` (find this file in the folder fuma_ct)
    - On line 4, edit the list `scrnas`. If your scRNAseq file name is in this format: `1_Allen_MCA_Human_2019.h5ad`, then the item in the `scrna` list is: `1_Allen_MCA_Human_2019`

## Step 4: Overview of outputs
- From the snakemake file `fuma_celltype_step1.smk`
    - Rule `magma_step1` generates the file ending in `_magma.gsa.out`
    - Rule `correct_pval_magma_step1` generates the file ending  in `_magma.gsa.out.rmComment`
        - From this file, you can calculate adjusted p values from the `P` column

# cellect
## Step 1: Download software
- Please follow instructions from https://github.com/ctg-single-cell/prioritize_ct_methods_comp/tree/main/codes/cellect in order to download cellect and set up the conda environment

## Step 2: Process GWAS
- Essentially cellect uses a similar script to LDSC to process GWAS but it's slightly different. 
- Check an example here: `/gpfs/work5/0/vusr0480/Processed_GWAS/for_cellect/code/process_gwas.sh`
    - Please check this directory here `/gpfs/work5/0/vusr0480/Processed_GWAS/for_cellect/out/` in case your trait of interest has already been processed

## Step 3: Run cellex from scRNAseq
- cellect uses "specificity" computed from cellex
- For some scRNAseq data, I have already computed specificity using cellex and the output files are located here: `/gpfs/work5/0/vusr0480/Processed_scRNA/data/cellex/`
- For any new scRNAseq that has not been processed yet, please use the script `process_scrna_cellex_alc.py` as a guide (find this script in the attached folder)

## Step 4: Run cellect
- As detailed in `Step 3: run snakemake` from the readme (https://github.com/ctg-single-cell/prioritize_ct_methods_comp/blob/main/codes/cellect/README.md), you need to edit the config file and then run the snakemake file from cellect.
- To submit cellect snakemake file to snellius, this is an example I use:

```
conda activate cellect
cd CELLECT/ #this is the path to the CELLECT github repo
snakemake --use-conda -j 10 -s /gpfs/home6/tphung/software/CELLECT/cellect-ldsc.snakefile --configfile /gpfs/home6/tphung/projects/alc_jeanne/code/cellect/config_ldsc_alc_jeanne.yml --rerun-incomplete -k
```

- Where is the output?
    - The output would be saved in `{BASE_OUTPUT_DIR}/CELLECT-MAGMA/results/prioritization.csv` or `{BASE_OUTPUT_DIR}/CELLECT-LDSC/results/prioritization.csv` where `{BASE_OUTPUT_DIR}` is specified in the config file for cellect

## Step 5: Postprocess
- Since the results for cellect for all traits and all scRNAseq datasets are all in one file, I wrote an R script (`postprocess_cellect_nonABA.R`) to tabulate the results for each trait across all scRNAseq datasets for each level (p values were adjusted within dataset) (find this file in the folder cellect)
- The snakemake file `postprocess_cellect_nonABA.smk` implements the Rscript (find this file in the folder fuma_ct)
  
## MAGMA CellTyping top 10
- Brief background: MAGMA CellTyping is a method where specificity is the "vanilla" calculation of specificity (for each gene, mean expression per cell type is divided over sum expression across cell types for that gene). Briefly, here are the steps:
    - Use the R package `EWCE` to compute specificity: I was able to implement this on my end
    - Use the R package MAGMA.Celltyping to run: I was not able to get this to work so I extract out the relevant functions in R to run this
### Detailed instruction
#### Step 1: Process GWAS
- For this, you can use the processed GWAS that was done previously for `fuma_ct` (see above)

### Step 2: Process scRNAseq data 
- The rules for this step are included in the snakemake file `magma_celltyping.smk` (find this script in the attached folder)
- Here, we will be using the EWCE package in R to compute specificity
    - Rule: `process_scrna`
        - Please edit the python script `process_scrna_magma_celltying.py` if necessary. This script is used to filter the h5ad file before running EWCE (find this file in the folder celltyping_top10)
    - Rule: `run_ewce`
        - Please edit the R script `process_scrna_ewce.R` if necessary. This script is calling the R package EWCE to compute specificity in a format that is compatible with MAGMA CellTyping (find this file in the folder celltyping_top10)

### Step 3: Run MAGMA
- The rules for this step are included in the snakemake file `magma_celltyping.smk` (find this file in the folder celltyping_top10)
- Rule: `generate_genecovar_top10`
    - This rule run the R script `create_geneCovar_top10.R`. 
        - Please edit this script as necessary.
        - Specifically, on line 27, please edit the path to the MAGMA_Celltyping repo
- Rule: `run_magma_top10`
    - Please edit if necessary

### Step 4: Postprocess
- The riles for this step are included in the snakemake file `magma_celltyping.smk` (find this file in the folder celltyping_top10)
- Rule: `remove_comments`
    - This rule basically just remmove the character # in the magma output file. Please edit if necessary
- Rule: `correct_pval`
    - This rule is the same as the rule in the fuma_ct `fuma_celltype_step1.smk` file