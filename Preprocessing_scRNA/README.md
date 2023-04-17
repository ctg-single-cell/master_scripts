# Preprocess publicly available scRNAseq data
## First step
- Preprocess definition: 
  - Document relevant information in a master spreadsheet `scrnaseq_data_master.csv`
    - Each scRNAseq dataset has an `Internal_ID` such as `idX` where X is the next number in the spreadsheet
    - Naming convention for each scRNAseq dataset:
      - {Author/Group}_{Tissue/Area}_{Year}_{Extra}
      - Example: `Allen_MCA_2019`, `Siletti_Cerebellum_2022`
  - Basic filtering:
    - Filter #1: keep cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    - Filter #2: keep cells where mt percentage is less than 10% for human and 5% for mouse
  - Convert gene name convention to ensembl
    - Here, genes that do not have an ensembl id will be annotated as NaN. These genes are kept (instead of removed) because it could be that we might need to convert the original gene name to something else in the future, so it's best to keep as NaN.
- Initial setup:
  - For documenting, the working directory is: `/project/prjstphung/Preprocessing_scRNA/`.
  - This is the file structure in this directory:
  ```commandline
  |-code
  | |-initial_setup.py
  | |-qc_scrna_id1.py
  | |-single_cell_helper_functions_v3.py
  | |-conversion_files
  | | |-gene_names_human_incl_rat.txt
  | | |-gene_names_human.txt
  | | |-gene_names_mouse.txt
  |-data
  |-scrnaseq_data_master.csv
  ```
  - For more details: 
    - the scripts `initial_setup.py` and `qc_scrna_id1.py` can be found in this repo
    - copy the latest version of `single_cell_helper_functions_v3.py` from https://github.com/ctg-single-cell/master_scripts/blob/main/Conversion_gene_names/single_cell_helper_functions_v3.py )#TODO: to make this more automatic)
    - copy over the converion files as well

## Step by step instruction for each scRNAseq
1. Setup
- Create a new entry in the excel spread sheet `scrnaseq_data_master.csv` with `idX` where X is the next number.
- Run Python script:
```
pwd
/project/prjstphung/Preprocessing_scRNA
python code/initial_setup.py --id 1 --work_dir data/
```
- What does the above script do? 
  + create a folder: `id1` in `{working_dir}/data`. Here, `{working_dir}` is `/project/prjstphung/Preprocessing_scRNA/`
  + create a `readme.md` in `{working_dir}/data/id1`
  + updated structure:
  ```
   |-code
   | |-initial_setup.py
   | |-qc_scrna_id1.py
   | |-single_cell_helper_functions_v3.py
   | |-conversion_files
   | | |-gene_names_human_incl_rat.txt
   | | |-gene_names_human.txt
   | | |-gene_names_mouse.txt
   |-data
   | |-id1
   | | |-readme.md
   |-scrnaseq_data_master.csv
  ```

3. Download scrnaseq data to `{working_dir}/data/id1`
- Document downloading of data in `{working_dir}/data/id1/readme.md`
- View the data, etc... in order to understand the data structure so that you can update the script `qc_scrna_id1.py` (here, this script has already been modified for this particular data)
- A few key things to know about the scRNAseq data:
    + what is the format of the data? Is this `csv`, `h5ad`, `loom`? 
    + make sure that `adata.var` contains information for genes
    + make sure that `adata.obs` contains information for cells 
    + is the cell type annotation in the metadata? if so, keep track of what the column name for this
    + what is the gene name convention (symbol, ensembl, or other?)
- Document as much information as possible in `{working_dir}/data/id1/readme.md`
  
4. Customize the python script `qc_scrna_id1.py` 
- #TODO: In the future, we will have a more general script `qc_scrna_general.py` where one would need to make a copy of and edit information as appropriate.

5. Run
```
python code/qc_scrna_id1.py > code/log_qc_scrna_id1.txt
```
- The standard out here is the log from the conversion script (i.e. when we call `adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add)`). #TODO: need to fold this in to the job log.

6. Examination of outputs
- `plot1.jpeg`: highest expressed gene plot
- `plot2.jpeg`: violin plots of n_genes_by_counts, total_counts, pct_counts_mt, stratified by cell types
- `table2.csv`: a table recording the number of genes and cells originally and after each filter
- `table3.csv`: a table tabulating the number of cells per cell type
- `log.txt`: log file
- `{id}X.h5ad`: cleaned h5ad file
  - post applying the 2 above filters