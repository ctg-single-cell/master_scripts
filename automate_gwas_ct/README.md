- For examples please see the shared directory on snellius: `/gpfs/work5/0/vusr0480/Pipelines/gwas_ct_pipeline/`
- Repo structure: 
    - snakemake and config files are found under `master_scripts/automate_gwas_ct`
    - other helper scripts that are mentioned throughout this documentation can be found under `master_scripts/automate_gwas_ct/helper_scripts`. If there is a script that is mentioned in this readme but can't be found, please contact Tanya
## ldsc
### specificity
- The in-house pipeline implements 3 ways of computing specificity called `spec_vanilla`, `spec_ewce`, and `spec_cellex`. For more information about how to compute specificity, please refer to the file `compute_spec.ipynb`

- Step 1: Set up directoy
    - TL;DR:
    
    ```
    mkdir ldsc
    for i in spec_vanilla spec_ewce spec_cellec; do
        mkdir ldsc/${i}
        for j in 41_Siletti_Cerebellum.CBV_Human_2022 #add here additional scRNAseq if desired
            mkdir ldsc/${i}/${j}
            mkdir ldsc/${i}/${j}/genelists
    ``` 
    - Detailed explanation:
        - Create a directory named `ldsc`
        - Subdirectories within `ldsc`: `spec_vanilla`, `spec_ewce`, `spec_cellex`
        - Within each subdirectory, create a directory for each scrnaseq dataset
        - The subdirectories under `spec_vanilla` is:
            - `genelists`: this is storing the genes that are specific for each cell type
            - `LD_annotations`: output from running ldsc (you don't need to create this because it will be automatically created from the snakemake pipeline)
            - `ldsc-seq_out`: final output from running ldsc (you don't need to create this because it will be automatically created from the snakemake pipeline)
        - An example of a directory structure:
        ``` 
        ├─ ldsc
        │   ├── spec_vanilla
        │   │   ├── 41_Siletti_Cerebellum.CBV_Human_2022
        │   │   │   ├── genelists
        │   │   │   ├── LD_annotations
        │   │   │   ├── ldsc-seq_out
        │   ├── spec_ewce
        │   │   ├── 41_Siletti_Cerebellum.CBV_Human_2022
        │   │   │   ├── genelists
        │   │   │   ├── LD_annotations
        │   │   │   ├── ldsc-seq_out
        │   ├── spec_cellex
        │   │   ├── 41_Siletti_Cerebellum.CBV_Human_2022
        │   │   │   ├── genelists
        │   │   │   ├── LD_annotations
        │   │   │   ├── ldsc-seq_out
        ```

- Step 2: Create genelists 
    - Follow the instruction from the file `compute_spec.ipynb` in order to generate the genelists 
    - Then, copy to the appropriate folder
- Step 3: Edit the confile file `pipeline_config.json`. Specifically: 
    - Line 3: `ldsc_path`
    - Line 5: `gwas_dir`
    - Line 7: `traits`
    - Line 8: `scrnas`
    - Line 9: `scrna_dir`
    - Line 10: `ldsc_dir` (this is the path to the directory we created in Step 1)
- Step 4: Create a text file that contain a list of cell types for the scRNAseq data that you are working on
    - For example, for the scRNAseq data `41_Siletti_Cerebellum.CBV_Human_2022`:
    
    ```
    cat cell_type_list.txt
    Oligodendrocyte
    Committed_oligodendrocyte_precursor
    Oligodendrocyte_precursor
    Splatter
    Upper_rhombic_lip
    Cerebellar_inhibitory
    Miscellaneous
    Astrocyte
    Bergmann_glia
    Ependymal
    Vascular
    Microglia
    Fibroblast
    ``` 
- Step 5: Edit the snakemake file `pipeline_ldsc_part1.smk`
    - Edit line 38, 53, and 70 to give the correct path
    - Check to make sure that there are no any other spots with hard-coded path
- Step 6: Edit the run file `run_smk.sh`
    - If necessary, change the number of node and time
    - If necessary, change the name of the conda environment
    - Update `--config cell_type_list` to the correct path
- Step 7: Postprocess results from running ldsc by editting and running the snakemake file `pipeline_ldsc_part2.smk`
    - Edit line 12 to give the correct path
    - In this snakemake pipeline, in the rule `aggregate_cts`, we are aggregating results across all cell types. Then, in the rule `compute_p`, we are computing P values from coefficient z score and compute adjust p values using bonferroni correction


## magma
- For magma, we will do the following: 
    - Using mean gene expression per cell type, run magma gene property (my naming convention is to call this `linear`)
    - Using specificity:
        - Run magma gene property (`linear`)
        - Obtain the top 10% (or whatever percentage) run magma gene set (`top10`)
    - In brief, there are 2 ways to compute mean gene expression per cell type (`vanilla` and `ewce` - check below for more information) and 3 ways to compute specificty (`vanilla`, `ewce`, and `cellex` - check the section above for more information). Then these are the directories we are going to have under `magma`:
        -  `mean_vanilla_linear`
        -  `mean_ewce_linear`
        -  `spec_vanilla_top10`
        -  `spec_ewce_top10`
        -  `spec_cellex_top10`
        -  `spec_vanilla_linear`
        -  `spec_ewce_linear`
        -  `spec_cellex_linear`

- Set up directoy
```
mkdir magma
for i in mean_vanilla_linear mean_ewce_linear spec_vanilla_top10 spec_ewce_top10 spec_cellex_top10 spec_vanilla_linear spec_ewce_linear spec_cellex_linear; do
    mkdir magma/${i} 
    for j in 41_Siletti_Cerebellum.CBV_Human_2022 #add here additional scRNAseq if desired
        mkdir magma/${i}/${j}
        mkdir magma/${i}/${j}/gene_covar
        mkdir magma/${i}/${j}/out
``` 
**NOTE**: The above code is just for demonstration. One should put this in a bash script to automate it. For simplicity I am not providing that script here. 
- Detailed explanation:
    - From your main directory, create a directory named `magma`
    - Subdirectories within `magma` (for mean): `mean_vanilla_linear`, `mean_ewce_linear`, `spec_vanilla_top10`, `spec_ewce_top10`, `spec_cellex_top10`, `spec_vanilla_linear`, `spec_ewce_linear`, `spec_cellex_linear`
    - Within each subdirectory, create a directory for each scrnaseq dataset
    - The subdirectories under `magma_vanilla_linear` is:
        - `gene_covar`: this is storing the covariate file used in magma
        - `out`: this is where the output is going to be
    ``` 
    ├─ magma
    │   ├── mean_vanilla_linear
    │   │   ├── 41_Siletti_Cerebellum.CBV_Human_2022
    │   │   │   ├── gene_covar
    │   │   │   ├── out
    │   ├── mean_ewce_linear
    │   │   ├── 41_Siletti_Cerebellum.CBV_Human_2022
    │   │   │   ├── gene_covar
    │   │   │   ├── out
    ```

### mean
- The in-house pipeline implements 2 ways of computing mean called `mean_vanilla`, and `mean_ewce`. For more information about how to compute mean, please refer to the file `compute_mean.ipynb`

- Step 1: Create gene_covar
    - Follow the instruction from the file `compute_mean.ipynb` in order to generate the gene covariate file
    - Then, copy to the appropriate folder. Ideally you should set it up so that the output after computing mean should be put directly into the folder.
- Step 2: Edit the confile file `pipeline_config.json`. Specifically: 
    - Line 2: `magma_executable`
    - Line 6: `magma_dir`
    - Line 7: `traits`
    - Line 8: `scrnas`
    - Line 9: `scrna_dir`
- Step 3: Edit the snakemake file `pipeline_magma_mean.smk`
    - Read over the snakemakefile `pipeline_magma_mean.smk`. In theory there should be nothing to be changed but double-check
- Step 4: Edit the run file `run_smk.sh`
    - If necessary, change the number of node and time
    - If necessary, change the name of the conda environment
    - Command to call snakemake: `snakemake -s pipeline_magma_mean.smk --configfile pipeline_config.json -j --rerun-incomplete -k`
- Step 5: Where to find outputs: 
    - Outputs should be found in for example: `magma/mean_vanilla_linear/41_Siletti_Cerebellum.CBV_Human_2022/out`:
    
    ```
    head magma/mean_vanilla_linear/41_Siletti_Cerebellum.CBV_Human_2022/out/SZ3_mean_vanilla_linear_magma.gsa.out.rmComment
    VARIABLE                             TYPE  NGENES         BETA     BETA_STD           SE            P FULL_NAME
    Committed_oligodendrocyte_pr...     COVAR   16741   -0.0021193   -0.0033898     0.010787      0.57787 Committed_oligodendrocyte_precursor
    Oligodendrocyte                     COVAR   16741     0.020833     0.031745    0.0099088     0.017765 Oligodendrocyte
    Oligodendrocyte_precursor           COVAR   16741     0.035745     0.056809     0.012469     0.002077 Oligodendrocyte_precursor
    Splatter                            COVAR   16741    0.0053039     0.010161     0.006445      0.20527 Splatter
    Upper_rhombic_lip                   COVAR   16741     0.052883     0.078204     0.011051   8.6037e-07 Upper_rhombic_lip
    Cerebellar_inhibitory               COVAR   16741     0.035635     0.071215     0.008593   1.6929e-05 Cerebellar_inhibitory
    Miscellaneous                       COVAR   16741    0.0053327    0.0061336     0.012038      0.32889 Miscellaneous
    Astrocyte                           COVAR   16741   -0.0011997   -0.0016138      0.01218      0.53923 Astrocyte
    Bergmann_glia                       COVAR   16741     0.026267     0.040584     0.011471      0.01102 Bergmann_glia
    ```


### specificity
- Similar to the ldsc pipeline, the in-house pipeline implements 3 ways of computing specificity. For more information about how to compute specificity, please refer to the file `compute_spec.ipynb`
    - Note that the file formats for `ldsc` and `magma` are different. Therefore, please refer to the file `compute_spec.ipynb` so that this can be set up properly

- Step 1: Create gene_covar
    - Follow the instruction from the file `compute_spec.ipynb` in order to generate the gene covariate file
    - Then, copy to the appropriate folder. Ideally you should set it up so that the output after computing mean should be put directly into the folder.
- Step 2: check and edit the confile file `pipeline_config.json` if necessary. Specifically: 
    - Line 2: `magma_executable`
    - Line 6: `magma_dir`
    - Line 7: `traits`
    - Line 8: `scrnas`
    - Line 9: `scrna_dir`
- Step 3: Edit the snakemake file `pipeline_magma_spec.smk`
    - Read over the snakemakefile `pipeline_magma_spec.smk`. In theory there should be nothing to be changed but double-check
- Step 4: Edit the run file `run_smk.sh`
    - If necessary, change the number of node and time
    - If necessary, change the name of the conda environment
    - Command to call snakemake: `snakemake -s pipeline_magma_spec.smk --configfile pipeline_config.json -j --rerun-incomplete -k`
- Step 5: Where to find outputs: 
    - Outputs should be found in for example: `magma/spec_vanilla_linear/41_Siletti_Cerebellum.CBV_Human_2022/out`:
    
    ```
    head magma/spec_vanilla_top10/41_Siletti_Cerebellum.CBV_Human_2022/out/SZ3_spec_vanilla_top10_magma.gsa.out.rmComment
    VARIABLE                             TYPE  NGENES         BETA     BETA_STD           SE            P FULL_NAME
    Committed_oligodendrocyte_pr...       SET     898    -0.042591   -0.0090694     0.037886      0.86953 Committed_oligodendrocyte_precursor
    Oligodendrocyte                       SET    1016     0.071267     0.016089     0.037068     0.027273 Oligodendrocyte
    Oligodendrocyte_precursor             SET     970     0.078136     0.017258     0.040024     0.025462 Oligodendrocyte_precursor
    Splatter                              SET     353   -0.0073949   -0.0010021     0.059387      0.54955 Splatter
    Upper_rhombic_lip                     SET    1052      0.10493     0.024081     0.036453    0.0019999 Upper_rhombic_lip
    Cerebellar_inhibitory                 SET    1043     0.093274     0.021319     0.036935    0.0057828 Cerebellar_inhibitory
    Miscellaneous                         SET     971    -0.041474   -0.0091648     0.037633      0.86478 Miscellaneous
    Astrocyte                             SET     976    -0.054309     -0.01203     0.036983      0.92901 Astrocyte
    Bergmann_glia                         SET     975     0.067776     0.015006     0.038206     0.038043 Bergmann_glia
    ```