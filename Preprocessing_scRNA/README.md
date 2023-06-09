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
    - An example location of the "preprocesed" scRNAseq data: `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/219_Caoetal2020_MuscleOrgan_Human_2020_Level2/219_Caoetal2020_MuscleOrgan_Human_2020_Level2.h5ad`
    - Brief description of the layers in this `h5ad` file:
      - Example:
      ```
      # change directory to: `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/210_Caoetal2020_Eye_Human_2020_Level2/` on snellius
      adata = anndata.read("210_Caoetal2020_Eye_Human_2020_Level2.h5ad")
      ```
      - `adata.X`: this is the raw count. In the preprocessing steps, if the raw counts were stored in the `adata.raw.X` layer, we update the `adata.X` layer to be the raw count.
      ```
      adata.X.A[1,50:300]
      array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
       0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 4., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], dtype=float32)
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
                  feature_id  feature_is_filtered feature_reference feature_biotype     mt  n_cells_by_counts  mean_counts  pct_dropout_by_counts  total_counts  n_cells          ensembl
        symbol
        SLC16A13     ENSG00000174327                False    NCBITaxon:9606            gene  False                  6     0.000468              99.953242           6.0        5  ENSG00000174327
        MMP25        ENSG00000008516                False    NCBITaxon:9606            gene  False                  3     0.000234              99.976621           3.0        3  ENSG00000008516
        MAPK8IP1     ENSG00000121653                False    NCBITaxon:9606            gene  False                 97     0.008339              99.244077         107.0       95  ENSG00000121653
        FAF1         ENSG00000185104                False    NCBITaxon:9606            gene  False               2761     0.282185              78.483479        3621.0     2685  ENSG00000185104
        RP1-30E17.2  ENSG00000225689                False    NCBITaxon:9606            gene  False                 63     0.005455              99.509040          70.0       63              NaN
        ...                      ...                  ...               ...             ...    ...                ...          ...                    ...           ...      ...              ...
        NRG3         ENSG00000185737                False    NCBITaxon:9606            gene  False               4877     1.097413              61.993454       14082.0     4687  ENSG00000185737
        CSTA         ENSG00000121552                False    NCBITaxon:9606            gene  False                  4     0.000312              99.968828           4.0        3  ENSG00000121552
        ERMN         ENSG00000136541                False    NCBITaxon:9606            gene  False                  4     0.000312              99.968828           4.0        4  ENSG00000136541
        CACNG8       ENSG00000142408                False    NCBITaxon:9606            gene  False                113     0.009430              99.119389         121.0      110  ENSG00000142408
        N4BP1        ENSG00000102921                False    NCBITaxon:9606            gene  False                343     0.030938              97.326995         397.0      333  ENSG00000102921
        ```
        + Above, the column `feature_id` was the ensembl column that came originally with the scRNAseq data. The `ensembl` column is the one we added using the conversion script (see lines 15 above for explanation)

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
