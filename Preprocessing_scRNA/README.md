# Preprocess publicly available scRNAseq data
## General information
- Document relevant information in a master spreadsheet `scrnaseq_data_master.csv`. 
  - This document is stored & updated locally on Tanya's computer. The latest version could also be found on snellius: `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/scrnaseq_data_master.csv`. **NOTES**: it is very important to make sure to sync the local version to the one on snellius.
  - Each scRNAseq dataset has an `Internal_ID` such as `X_Allen_MCA_Human_2019` where:
    - X is the next number in the spreadsheet
    - {FirstAuthor/Consortium/Group}_{Tissue/Region/Area}_{Species}_{Year}
- General workflow:
  - If the downloaded scRNAseq data contain cells for all tissues, subset per tissue first
  - Basic filtering per tissue/region:
    - Filter #1: keep cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    - Filter #2: keep cells where mt percentage is less than 10% for human and 5% for mouse
  - Convert gene name convention to ensembl
    - Here, genes that do not have an ensembl id will be annotated as NaN. These genes are kept (instead of removed) because it could be that we might need to convert the original gene name to something else in the future, so it's best to keep as NaN.
    - For scRNAseq datasets that include both the ensembl id and gene symbol in the `obs` layer, I decided to convert to ensemble id using the gene symbol because I find that our converted ensembl id gene list is typically fewer than the one that is included with the paper. This is probably due to the fact that the conversion file we are using have to have the gene symbols be approved. 
  - Save the "preprocessed" scRNAseq data:
    - This is the file that one should use to calculate different statistics prior to going to the different pipeline (such as mean or specificity)
    - Naming convention: {Internal_ID}.h5ad where `Internal_ID` is in the spreadsheet `scrnaseq_data_master.csv`.
    - An example location of the "preprocesed" scRNAseq data: `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/1_Allen_MCA_Human_2019/1_Allen_MCA_Human_2019.h5ad`
    - Brief description of the layers in this `h5ad` file:
      - Example:
      ```
      # change directory to: `gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/1_Allen_MCA_Human_2019/` on snellius
      adata = anndata.read("1_Allen_MCA_Human_2019.h5ad")
      ```
      - `adata.X`: this is the raw count. In the preprocessing steps, if the raw counts were stored in the `adata.raw.X` layer, we update the `adata.X` layer to be the raw count.
      ```
      adata.X
      array([[  0.,   0.,   0., ...,   0.,   0.,   1.],
       [  0.,   0.,   0., ..., 748.,   0.,   0.],
       [  0.,   0.,   0., ...,   0.,   0.,   0.],
       ...,
       [  0.,   0.,   0., ..., 129.,   0.,   0.],
       [  0.,   0.,   0., ..., 376.,   0.,   0.],
       [  0.,   0.,   0., ...,   0.,   0.,   0.]], dtype=float32)
      ```
        + Here, we assume that since the values are whole integer this is probably the raw count (we could never be 100% sure though).
      - `adata.obs`: this is the observation dataframe where each row contains the metadata for each cell. 
        + This dataframe should be the same as the original with the addition of 5 columns created from the command `sc.pp.calculate_qc_metrics`"
          + `n_genes_by_counts`: the number of genes with at least 1 count in a cell
          + `total_counts`: total number of counts for a cell 
          + `total_counts_mt`
          + `pct_counts_mt`
          + `n_genes`
        + The column for cell type annotation is not changed here because there could be different levels for annotation. Rather, please refer to the readme for each specific dataset to know what the column name is for the cell type. 
      - `adata.var`: this is the dataframe storing the information for the genes:
        + row index is the gene symbol with the name of the row index being `symbol`. 
        + the converted ensembl id is stored in the column with column name `ensembl`. This is the one we should use for downstream analyses
        ```
        adata.var
                        mt  n_cells_by_counts  mean_counts  pct_dropout_by_counts  total_counts  n_cells          ensembl
        symbol
        3.8-1.2      False                 23     0.012344              99.953457         610.0       23              NaN
        3.8-1.3      False                 48     0.065241              99.902867        3224.0       48              NaN
        3.8-1.4      False                 27     0.016796              99.945363         830.0       27              NaN
        3.8-1.5      False                 28     0.004938              99.943339         244.0       28              NaN
        5-HT3C2      False               1093     0.238481              97.788211       11785.0     1093  ENSG00000233576
        ...            ...                ...          ...                    ...           ...      ...              ...
        ZYX          False              24628    33.814396              50.162899     1671006.0    24628  ENSG00000159840
        ZZEF1        False              37165    65.633690              24.793087     3243420.0    37165  ENSG00000074755
        ZZZ3         False              40929   166.196274              17.176275     8212921.0    40929  ENSG00000036549
        bA255A11.4   False                 85     0.019002              99.827994         939.0       85              NaN
        bA395L14.12  False               6087     0.228363              87.682377       11285.0     6087              NaN
        ```

- Initial setup:
  - For documenting, the working directory is: `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/` (on snellius).


## Step by step instruction for each scRNAseq
1. Setup
- Create a new entry in the excel spread sheet `scrnaseq_data_master.csv` with `X_{FirstAuthor/Consortium/Group}_{Tissue/Region/Area}_{Species}_{Year}` where X is the next number.
- Run snakemake script:
```
pwd
/project/prjstphung/Preprocessing_scRNA
snakemake -s snakefile_initial_setup.smk -c1
```
- The snakemake script calls the python script `initial_setup.py` to do the following:
  + create a folder: `X_{FirstAuthor/Consortium/Group}_{Tissue/Region/Area}_{Species}_{Year}` in `{working_dir}/data`. Here, `{working_dir}` is `/project/prjstphung/Preprocessing_scRNA/`
  + create a `readme.md` in `{working_dir}/data/7_Siletti_Cerebellum_Human_2022`
  + the beauty of snakemake here is that it (hopefully) will create a new folder for new rows that we are adding. But to be safe, always run a dry-run first to make sure that it will only create new directories for the new rows. 

3. Download scrnaseq data to `{working_dir}/data/7_Siletti_Cerebellum_Human_2022`
- Document downloading of data in `{working_dir}/data/7_Siletti_Cerebellum_Human_2022/readme.md`
- Exception: sometimes the scRNAseq data contains information for all of the cells from multiple areas. For our purpose, we would like to separate them out by tissue/region first. This is an example:
  + Create folder `Caoetal_2020_Fetal` in `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/`
  + Download the `h5ad` file to `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/Caoetal_2020_Fetal`
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

5. Use snakemake 
- Detailed instruction to be added later but for an example, see: `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/code/snakefile_qc_scrna_caoetal2020.smk`

6. Examination of outputs
- `{Internal_ID}.h5ad`: processed scRNAseq data
- `plot1.jpeg`: highest expressed gene plot
- `{Internal_ID}_obs_for_plot2.csv`: data for plotting in R violin plots of n_genes_by_counts, total_counts, pct_counts_mt, stratified by cell types
- `table2.csv`: a table recording the number of genes and cells originally and after each filter
- `table3.csv`: a table tabulating the number of cells per cell type
- `log.txt`: log file
