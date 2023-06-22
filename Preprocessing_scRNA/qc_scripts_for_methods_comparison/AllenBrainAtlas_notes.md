---
output:
  pdf_document: default
  html_document: default
---
# Documentation on preprocessing Allen Brain Atlas
- The Allen Brain Atlas (ABA) data will be used for the method comparison
## General information
- Webpage: http://portal.brain-map.org/atlases-and-data/rnaseq
- Summary for the manuscript (Supplementary Table 1): https://vunl-my.sharepoint.com/:x:/r/personal/t_n_phung_vu_nl/Documents/Projects/Prioritize_CellType_Methods_Comparison/manuscript/si_materials.xlsx?d=w6131b9cbf14442658d4159e19a419747&csf=1&web=1&e=Bz0SDf
- In summary, there are 3 datasets for human and 2 datasets for mouse
    + multiple cortical areas (smart-seq 2019):
        + This is the dataset where they combine from previous datasets (MTG smart-seq 2018 and V1, ACC smart-seq 2018). Therefore, we will not be using the ones from 2018
        + Internal ID in the spreadsheet `scrnaseq_data_master.csv`: `1_Allen_MCA_Human_2019`
    + m1 (10x 2020)
        + Internal ID: `2_Allen_M1_Human_2020`
    + mtg (10x 2022)
        + part of the SEA-AD studies
        + Internal ID: `3_Allen_MTG_Human_2022`
    + whole cortex and hippocampus (smart-seq 2019)
        + Internal ID: `4_Allen_Cortex.Hippocampus_Mouse_2019`
    + whole cortex and hippocampus (10x 2020)
        + Internal ID: `5_Allen_Cortex.HippoCampus_Mouse_2020`
        + **NOTES**: the preprocessing for this is not yet available because of the large file size
- Location for preprocessed h5ad files (you can use these files as inputs to calculate mean, specifity, etc...)
    + `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/{Interal_ID}.h5ad`

## Description of the preprocessed h5ad
- This section is to describe the preprocessed h5ad file (i.e. after running the qc scripts below)
- Use example from: `1_Allen_MCA_Human_2019`
    - Preprocessed h5ad location on snellius: `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/1_Allen_MCA_Human_2019/1_Allen_MCA_Human_2019.h5ad`
- Matrix layer:

```
adata.X
<49417x48366 sparse matrix of type '<class 'numpy.float32'>'
        with 404815336 stored elements in Compressed Sparse Row format>

# This is the raw count because:
adata.X.A
array([[  0.,   0.,   0., ...,   0.,   0.,   1.],
       [  0.,   0.,   0., ..., 748.,   0.,   0.],
       [  0.,   0.,   0., ...,   0.,   0.,   0.],
       ...,
       [  0.,   0.,   0., ..., 129.,   0.,   0.],
       [  0.,   0.,   0., ..., 376.,   0.,   0.],
       [  0.,   0.,   0., ...,   0.,   0.,   0.]], dtype=float32)
```
- `obs` layer:

```
adata.obs
               sample_name   exp_component_name specimen_type cluster_color  ...  total_counts total_counts_mt pct_counts_mt  n_genes
0      F2S4_160113_027_A01  LS-15005h_S01_E1-50       nucleus           NaN  ...      975405.0             0.0           0.0     6778
1      F2S4_160113_027_B01  LS-15005h_S02_E1-50       nucleus       #E170FE  ...     1679385.0             0.0           0.0     7790
2      F2S4_160113_027_C01  LS-15005h_S03_E1-50       nucleus       #8E5864  ...     1600902.0             0.0           0.0     8696
3      F2S4_160113_027_D01  LS-15005h_S04_E1-50       nucleus       #8B5862  ...     1372513.0             0.0           0.0     7966
4      F2S4_160113_027_E01  LS-15005h_S05_E1-50       nucleus       #CF6EC9  ...     1689105.0             0.0           0.0     9377
...                    ...                  ...           ...           ...  ...           ...             ...           ...      ...
49412  F2S4_190227_100_C01  SM-GE4QU_S187_E1-50       nucleus       #312E27  ...      806351.0             0.0           0.0     4867
49413  F2S4_190227_100_E01  SM-GE4QU_S189_E1-50       nucleus       #BFC124  ...     1094163.0             0.0           0.0     7687
49414  F2S4_190227_100_F01  SM-GE4QU_S190_E1-50       nucleus       #8B5862  ...     1079012.0             0.0           0.0     6662
49415  F2S4_190227_100_G01  SM-GE4QU_S191_E1-50       nucleus       #71AF9A  ...     1016524.0             0.0           0.0     4245
49416  F2S4_190227_100_H01  SM-GE4QU_S192_E1-50       nucleus       #71AF9A  ...      924290.0             0.0           0.0     4197

[49417 rows x 46 columns]
```
  - cell id: is stored in the `sample_name` column
  - cell types: are stored in 3 layers:
      - `cell_type_level_1`: 
      
      ```
      adata.obs["cell_type_level_1"].unique()
      [NaN, 'GABAergic', 'Glutamatergic', 'Non-neuronal']
      Categories (3, object): ['GABAergic', 'Glutamatergic', 'Non-neuronal']
      ```
      - `cell_type_level_2`:
      
      ```
      adata.obs["cell_type_level_2"].unique()
      [NaN, 'VIP', 'LAMP5', 'IT', 'PAX6', ..., 'L5 ET', 'Pericyte', 'Endothelial', 'L4 IT', 'VLMC']
      Length: 20
      Categories (19, object): ['Astrocyte', 'Endothelial', 'IT', 'L4 IT', ..., 'Pericyte', 'SST', 'VIP', 'VLMC']
      ```
      
      - `cell_type_level_3`:  
      
      ```
      adata.obs["cell_type_level_3"].unique()
      [NaN, 'Inh L2-5 VIP TOX2', 'Inh L1 LAMP5 GGT8P', 'Inh L1 LAMP5 NDNF', 'Inh L1-3 VIP ZNF322P1', ..., 'Exc L3-5 FEZF2 DCN', 'Exc L4 RORB CCDC168', 'Exc L3 LINC00507 CTXN3', 'Exc L3 THEMIS PLA2G7', 'Exc L5 FEZF2 DYRK2']
      Length: 121
      Categories (120, object): ['Astro L1 FGFR3 FOS', 'Astro L1 FGFR3 MT1G', 'Astro L1-6 FGFR3 ETNPPL',
                              'Endo L2-5 CLDN5', ..., 'Oligo L4-6 MOBP COL18A1', 'Oligo L4-6 OPALIN',
                              'Peri L1-6 MUSTN1', 'VLMC L1-3 CYP1B1']
      ```
    - Additional columns added on 6-22-2023 to record whether a cell should be keep (True) or removed (False) based on the following conditions:
        - a cell is removed if the there is no cell type annotation (NaN)
        - a cell if removed if it's annotated to a cell type that only has 1 cell 
        - These columns are: `keep_level_1`, `keep_level_2`, and `keep_level_3`
        - Sanity check to make sure that this implementation is correct:
        
        ```
        # for 1_Allen_MCA_Human_2019, level_1 cell type, there are NaN
        adata = anndata.read("1_Allen_MCA_Human_2019.h5ad")
        adata.obs["cell_type_level_1"].isna().sum()
        1985
        adata.obs[adata.obs["keep_level_1"]==False].shape[0]
        1985
        ```
        

- `var` layer:
    - the row index is the gene symbol. The name of the row index is `symbol`
        - **IMPORTANT**: For the mouse dataset `4_Allen_Cortex.Hippocampus_Mouse_2019/4_Allen_Cortex.Hippocampus_Mouse_2019.h5ad`, please beware that the row index gene symbol is the **mouse** gene symbol. 
    - the gene name in ensembl format is in the column `ensembl`
    
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

    [48366 rows x 7 columns]
    ```

## Relevant information
- For each dataset, I also outputed a few files that could be useful to look at (I'm using these files for the Rshiny app to view the QC internally but it's not set up to share easily yet)
- `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/{Internal_ID}/table2.csv`: this file recorded the number of cells and genes before and after each filtering
- `/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/{Internal_ID}/table3.csv`: this file recorded the number of cells per cell type (all levels)
    

## 1_Allen_MCA_Human_2019
### Download:
```
pwd; date
/home/tphung/tphung_proj/Preprocessing_scRNA/data/1_Allen_MCA_Human_2019
Fri Jun  9 13:14:34 CEST 2023
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv

# download the metadata
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv

ls -lrth
total 5.1G
-rw-r--r-- 1 tphung prjs0480 5.1G Oct 20  2021 matrix.csv
-rw-r--r-- 1 tphung prjs0480  13M Oct 20  2021 metadata.csv
-rw-r--r-- 1 tphung prjs0480   15 May 17 09:24 readme.md
```

### Examine the data
```
# Since the downloaded data are in `csv` format and the metadata is in a different file, use the following code to read them into an adata object
adata = anndata.read_csv("matrix.csv")
meta = pd.read_csv("metadata.csv")
adata.obs = meta
```
- Examine the matrix
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
    + Since the values are integers, I assume that this is the raw count

- Examine the observation
```
adata.obs
               sample_name   exp_component_name  ... outlier_call              outlier_type
0      F2S4_160113_027_A01  LS-15005h_S01_E1-50  ...         True  Outlier L1-3 SST OR2AD1P
1      F2S4_160113_027_B01  LS-15005h_S02_E1-50  ...        False                       NaN
2      F2S4_160113_027_C01  LS-15005h_S03_E1-50  ...        False                       NaN
3      F2S4_160113_027_D01  LS-15005h_S04_E1-50  ...        False                       NaN
4      F2S4_160113_027_E01  LS-15005h_S05_E1-50  ...        False                       NaN
...                    ...                  ...  ...          ...                       ...
49412  F2S4_190227_100_C01  SM-GE4QU_S187_E1-50  ...        False                       NaN
49413  F2S4_190227_100_E01  SM-GE4QU_S189_E1-50  ...        False                       NaN
49414  F2S4_190227_100_F01  SM-GE4QU_S190_E1-50  ...        False                       NaN
49415  F2S4_190227_100_G01  SM-GE4QU_S191_E1-50  ...        False                       NaN
49416  F2S4_190227_100_H01  SM-GE4QU_S192_E1-50  ...        False                       NaN

[49417 rows x 41 columns]

# find columns
adata.obs.columns
Index(['sample_name', 'exp_component_name', 'specimen_type', 'cluster_color',
       'cluster_order', 'cluster_label', 'class_color', 'class_order',
       'class_label', 'subclass_color', 'subclass_order', 'subclass_label',
       'full_genotype_color', 'full_genotype_order', 'full_genotype_label',
       'donor_sex_color', 'donor_sex_order', 'donor_sex_label', 'region_color',
       'region_order', 'region_label', 'cortical_layer_color',
       'cortical_layer_order', 'cortical_layer_label',
       'cell_type_accession_color', 'cell_type_accession_order',
       'cell_type_accession_label', 'cell_type_alias_color', 'cell_type_order',
       'cell_type_alias_label', 'cell_type_alt_alias_color',
       'cell_type_alt_alias_order', 'cell_type_alt_alias_label',
       'cell_type_designation_color', 'cell_type_designation_order',
       'cell_type_designation_label', 'external_donor_name_color',
       'external_donor_name_order', 'external_donor_name_label',
       'outlier_call', 'outlier_type'],
      dtype='object')
```
    + `class_label`: 4 cell types
    + `subclass_label`: 20 cell types
    + `cluster_label`: 121 cell types 

- Examine the var
```
adata.var
Empty DataFrame
Columns: []
Index: [3.8-1.2, 3.8-1.3, 3.8-1.4, 3.8-1.5, 5-HT3C2, A1BG, A1BG-AS1, A1CF, A2M, A2M-AS1, A2ML1, A2MP1, A3GALT2, A4GALT, A4GNT, AA06, AAAS, AACS, AACSP1, AADAC, AADACL2, AADACL2-AS1, AADACL3, AADACL4, AADACP1, AADAT, AAED1, AAGAB, AAK1, AAMDC, AAMP, AANAT, AAR2, AARD, AARS, AARS2, AARSD1, AASDH, AASDHPPT, AASS, AATBC, AATF, AATK, AATK-AS1, ABAT, ABCA1, ABCA10, ABCA11P, ABCA12, ABCA13, ABCA17P, ABCA2, ABCA3, ABCA4, ABCA5, ABCA6, ABCA7, ABCA8, ABCA9, ABCA9-AS1, ABCB1, ABCB10, ABCB10P1, ABCB10P3, ABCB10P4, ABCB11, ABCB4, ABCB5, ABCB6, ABCB7, ABCB8, ABCB9, ABCC1, ABCC10, ABCC11, ABCC12, ABCC13, ABCC2, ABCC3, ABCC4, ABCC5, ABCC5-AS1, ABCC6, ABCC6P1, ABCC6P2, ABCC8, ABCC9, ABCD1, ABCD1P1, ABCD1P2, ABCD1P3, ABCD1P4, ABCD1P5, ABCD2, ABCD3, ABCD4, ABCE1, ABCF1, ABCF2, ABCF2P1, ...]
```
    + This is an empty dataframe with the row index being the gene symbol
- Script for preprocessing: `qc_scrna_id1.py`
- Table recording the number of cells and genes
  
|Description       | Number of cells      | Number of genes       |
|  ---  |  ---  |  ---  |
| Original      | 49,417 cells      | 50,281 genes      |
| First filter      | 49,417 cells       | 48,366 genes       | 
| Second filter      | 49,417 cells       | 48,366 genes       | 



## 2_Allen_M1_Human_2020
### Download the data
```
pwd; date
/home/tphung/tphung_proj/Preprocessing_scRNA/data/2_Allen_M1_Human_2020
Sun Jun 11 23:27:04 CEST 2023

wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv
```

### Examine the data
```
adata = anndata.read_csv("matrix.csv")
<!-- adata
AnnData object with n_obs × n_vars = 76533 × 50281 -->
meta = pd.read_csv("metadata.csv")
<!-- meta.shape
(76533, 39) -->
adata.obs = meta
```
- Examine the matrix
```
adata.X
array([[ 0.,  0.,  0., ..., 19.,  0.,  0.],
       [ 0.,  0.,  0., ...,  3.,  0.,  0.],
       [ 0.,  0.,  0., ...,  5.,  0.,  0.],
       ...,
       [ 0.,  0.,  0., ...,  0.,  0.,  0.],
       [ 0.,  0.,  0., ...,  5.,  0.,  0.],
       [ 0.,  0.,  0., ...,  0.,  0.,  0.]], dtype=float32)
```
    + Since the values are integers, I assume that this is the raw count

- Examine the observation
```
adata.obs
                               sample_name                      exp_component_name  ... outlier_call outlier_type
0      AAACCCAAGGATTTCC-LKTX_190129_01_A01  AAACCCAAGGATTTCC-21L8TX_180927_001_A01  ...        False          NaN
1      AAACCCAAGTATGGCG-LKTX_190129_01_A01  AAACCCAAGTATGGCG-21L8TX_180927_001_A01  ...        False          NaN
2      AAACCCACAAAGTGTA-LKTX_190129_01_A01  AAACCCACAAAGTGTA-21L8TX_180927_001_A01  ...        False          NaN
3      AAACCCACACTACTTT-LKTX_190129_01_A01  AAACCCACACTACTTT-21L8TX_180927_001_A01  ...        False          NaN
4      AAACCCACAGTGAGCA-LKTX_190129_01_A01  AAACCCACAGTGAGCA-21L8TX_180927_001_A01  ...        False          NaN
...                                    ...                                     ...  ...          ...          ...
76528  TTTGTTGAGATGGCGT-LKTX_190130_01_H01  TTTGTTGAGATGGCGT-35L8TX_181108_001_D01  ...        False          NaN
76529  TTTGTTGCACAGCCAC-LKTX_190130_01_H01  TTTGTTGCACAGCCAC-35L8TX_181108_001_D01  ...        False          NaN
76530  TTTGTTGCAGAGACTG-LKTX_190130_01_H01  TTTGTTGCAGAGACTG-35L8TX_181108_001_D01  ...        False          NaN
76531  TTTGTTGCATAATGAG-LKTX_190130_01_H01  TTTGTTGCATAATGAG-35L8TX_181108_001_D01  ...        False          NaN
76532  TTTGTTGTCTACTCAT-LKTX_190130_01_H01  TTTGTTGTCTACTCAT-35L8TX_181108_001_D01  ...        False          NaN

[76533 rows x 39 columns]

# find columns
adata.obs.columns
Index(['sample_name', 'exp_component_name', 'cluster_label', 'cluster_color',
       'cluster_order', 'class_label', 'class_color', 'class_order',
       'subclass_label', 'subclass_color', 'subclass_order', 'donor_sex_label',
       'donor_sex_color', 'donor_sex_order', 'region_label', 'region_color',
       'region_order', 'cortical_layer_label', 'cortical_layer_color',
       'cortical_layer_order', 'cell_type_accession_label',
       'cell_type_accession_color', 'cell_type_accession_order',
       'cell_type_alias_label', 'cell_type_alias_color',
       'cell_type_alias_order', 'cell_type_alt_alias_label',
       'cell_type_alt_alias_color', 'cell_type_alt_alias_order',
       'cell_type_designation_label', 'cell_type_designation_color',
       'cell_type_designation_order', 'external_donor_name_label',
       'external_donor_name_color', 'external_donor_name_order',
       'specimen_type', 'full_genotype_label', 'outlier_call', 'outlier_type'],
      dtype='object')
```
    + `class_label`: 3 cell types
    + `subclass_label`: 20 cell types
    + `cluster_label`: 127 cell types 

- Examine the var
```
adata.var
Empty DataFrame
Columns: []
Index: [DDX11L1, WASH7P, MIR6859-1, MIR1302-2, FAM138A, LOC105379212, OR4G4P, OR4G11P, OR4F5, LOC105379213, CICP27, LOC729737, LOC100996442, LOC105379214, LOC102725121, LOC102723897, MIR6859-2, RPL23AP21, LOC102723917, LOC105379431, RPL23AP24, LOC105379248, OR4F29, CICP7, LOC100132287, LOC100134822, LOC105378947, LOC101928626, MTND1P23, MTND2P28, MIR6723, OR4F16, LOC105378582, CICP3, LOC100133331, LOC100288069, LOC105378581, LOC100287934, LOC105378580, FAM87B, LINC00115, LINC01128, FAM41C, TUBB8P11, LOC284600, LOC100130417, SAMD11, NOC2L, KLHL17, PLEKHN1, PERM1, LOC105378583, HES4, RPL39P12, ISG15, AGRN, LOC100288175, LOC105369174, LOC105378948, RNF223, C1orf159, LOC105378584, LINC01342, MIR200B, MIR200A, MIR429, TTLL10-AS1, TTLL10, TNFRSF18, TNFRSF4, SDF4, B3GALT6, FAM132A, UBE2J2, LOC101928895, SCNN1D, ACAP3, MIR6726, PUSL1, CPSF3L, MIR6727, CPTP, TAS1R3, DVL1, MIR6808, MXRA8, AURKAIP1, CCNL2, LOC148413, MRPL20, ANKRD65, LOC105378585, TMEM88B, LOC102724312, VWA1, ATAD3C, ATAD3B, ATAD3A, TMEM240, SSU72, ...]
```
    + This is an empty dataframe with the row index being the gene symbol
- Script for preprocessing: `qc_scrna_id2.py`
- Table recording the number of cells and genes
  
|Description       | Number of cells      | Number of genes       |
|  ---  |  ---  |  ---  |
| Original      | 76,533 cells      | 50,281 genes      |
| First filter      | 76,533 cells       | 32,858 genes       | 
| Second filter      | 76,533 cells       | 32,858 genes       | 
| Ngenes with ensembl      | 76,533 cells       | 20,479 genes       | 

## 3_Allen_MTG_Human_2022
### Download the data
```
pwd; date
/home/tphung/tphung_proj/Preprocessing_scRNA/data/3_Allen_MTG_Human_2022
Mon Jun 12 10:00:43 CEST 2023

wget
```

### Examine the data
```
adata = anndata.read("Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad", backed="r")
<!-- adata
AnnData object with n_obs × n_vars = 166868 × 36601 backed at 'Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad'
    obs: 'sample_name', 'donor_sex_label', 'external_donor_name_label', 'species_label', 'age_label', 'region_label', 'cortical_layer_label', 'full_genotype_label', 'QCpass', 'cluster_label', 'cluster_confidence', 'subclass_label', 'subclass_confidence', 'class_label', 'class_confidence', 'GA_QCpass', 'GA_cluster_label', 'GA_subclass_label', 'GA_neighborhood_label', 'CA_QCpass', 'CA_cluster_label', 'CA_subclass_label', 'CA_neighborhood_label', 'cluster_color', 'cluster_order', 'subclass_color', 'subclass_order', 'class_color', 'class_order', 'GA_cluster_color', 'GA_cluster_order', 'GA_subclass_color', 'GA_subclass_order', 'CA_cluster_color', 'CA_cluster_order', 'CA_subclass_color', 'CA_subclass_order', 'cell_type_accession_label' -->
```
- Examine the matrix
```
adata.X[100,1500:2000].A
array([[0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 2., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
        1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 2., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 2., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 2., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 1.,
        0., 1., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 1., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 1., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 1., 0.,
        0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0.]], dtype=float32)
```
    + Since the values are integers, I assume that this is the raw count

- Examine the observation
```
adata.obs
                                                                                   sample_name donor_sex_label  ... CA_subclass_order cell_type_accession_label
specimen_name                                                                                                   ...                                          
AAACCCACAACTCATG-LKTX_191204_01_A01-1156636525  AAACCCACAACTCATG-LKTX_191204_01_A01-1156636525               M  ...              18.0             CS202204130_8
AAACCCACAATAGTAG-LKTX_191204_01_A01-1156636525  AAACCCACAATAGTAG-LKTX_191204_01_A01-1156636525               M  ...               NaN                       NaN
AAACCCACACGGTGTC-LKTX_191204_01_A01-1156636525  AAACCCACACGGTGTC-LKTX_191204_01_A01-1156636525               M  ...               7.0           CS202204130_105
AAACCCACACTCTGCT-LKTX_191204_01_A01-1156636525  AAACCCACACTCTGCT-LKTX_191204_01_A01-1156636525               M  ...              15.0            CS202204130_89
AAACCCACATCAGCAT-LKTX_191204_01_A01-1156636525  AAACCCACATCAGCAT-LKTX_191204_01_A01-1156636525               M  ...               9.0            CS202204130_93
...                                                                                        ...             ...  ...               ...                       ...
TTGGGCGTCCTCTTTC-L8TX_200107_01_A09-1156636575  TTGGGCGTCCTCTTTC-L8TX_200107_01_A09-1156636575               F  ...              15.0                       NaN
TTGGGCGTCGTGGGTC-L8TX_200107_01_A09-1156636575  TTGGGCGTCGTGGGTC-L8TX_200107_01_A09-1156636575               F  ...              15.0                       NaN
TTGGGTAGTTTGGGAG-L8TX_200107_01_A09-1156636575  TTGGGTAGTTTGGGAG-L8TX_200107_01_A09-1156636575               F  ...              14.0                       NaN
TTGTGGAAGTGGACGT-L8TX_200107_01_A09-1156636575  TTGTGGAAGTGGACGT-L8TX_200107_01_A09-1156636575               F  ...              14.0                       NaN
TTTCGATAGGAATTAC-L8TX_200107_01_A09-1156636575  TTTCGATAGGAATTAC-L8TX_200107_01_A09-1156636575               F  ...              23.0                       NaN

[166868 rows x 38 columns]

# find columns
adata.obs.columns
Index(['sample_name', 'donor_sex_label', 'external_donor_name_label',
       'species_label', 'age_label', 'region_label', 'cortical_layer_label',
       'full_genotype_label', 'QCpass', 'cluster_label', 'cluster_confidence',
       'subclass_label', 'subclass_confidence', 'class_label',
       'class_confidence', 'GA_QCpass', 'GA_cluster_label',
       'GA_subclass_label', 'GA_neighborhood_label', 'CA_QCpass',
       'CA_cluster_label', 'CA_subclass_label', 'CA_neighborhood_label',
       'cluster_color', 'cluster_order', 'subclass_color', 'subclass_order',
       'class_color', 'class_order', 'GA_cluster_color', 'GA_cluster_order',
       'GA_subclass_color', 'GA_subclass_order', 'CA_cluster_color',
       'CA_cluster_order', 'CA_subclass_color', 'CA_subclass_order',
       'cell_type_accession_label'],
      dtype='object')
```
    + `class_label`: 3 cell types, exluding nan (['Neuronal: GABAergic', nan, 'Neuronal: Glutamatergic', 'Non-neuronal and Non-neural'])
    + `subclass_label`: 24 cell types (excluding nan)
    + `cluster_label`: 127 cell types (excluding nan)

- Examine the var
```
adata.var
Empty DataFrame
Columns: []
Index: [MIR1302-2HG, FAM138A, OR4F5, AL627309.1, AL627309.3, AL627309.2, AL627309.5, AL627309.4, AP006222.2, AL732372.1, OR4F29, AC114498.1, OR4F16, AL669831.2, LINC01409, FAM87B, LINC01128, LINC00115, FAM41C, AL645608.6, AL645608.2, AL645608.4, LINC02593, SAMD11, NOC2L, KLHL17, PLEKHN1, PERM1, AL645608.7, HES4, ISG15, AL645608.1, AGRN, AL645608.5, AL645608.8, RNF223, C1orf159, AL390719.3, LINC01342, AL390719.2, TTLL10-AS1, TTLL10, TNFRSF18, TNFRSF4, SDF4, B3GALT6, C1QTNF12, AL162741.1, UBE2J2, LINC01786, SCNN1D, ACAP3, PUSL1, INTS11, AL139287.1, CPTP, TAS1R3, DVL1, MXRA8, AURKAIP1, CCNL2, MRPL20-AS1, MRPL20, AL391244.2, ANKRD65, AL391244.1, TMEM88B, LINC01770, VWA1, ATAD3C, ATAD3B, ATAD3A, TMEM240, SSU72, AL645728.1, FNDC10, AL691432.4, AL691432.2, MIB2, MMP23B, CDK11B, FO704657.1, SLC35E2B, CDK11A, SLC35E2A, NADK, GNB1, AL109917.1, CALML6, TMEM52, CFAP74, AL391845.2, GABRD, AL391845.1, PRKCZ, AL590822.2, PRKCZ-AS1, FAAP20, AL590822.1, SKI, ...]

[36601 rows x 0 columns]
```
    + This is an empty dataframe with the row index being the gene symbol
- Script for preprocessing: `qc_scrna_id3.py`
- Table recording the number of cells and genes
  
|Description       | Number of cells      | Number of genes       |
|  ---  |  ---  |  ---  |
| Original      | 166,868 cells      | 36,601 genes      |
| First filter      | 166,845 cells       | 34,557 genes       | 
| Second filter      | 164,053 cells       | 34,557 genes       | 
| Ngenes with ensembl      | 164,053 cells       | 22,654 genes       | 

## 4_Allen_Cortex.Hippocampus_Mouse_2019
### Download the data
```
pwd; date
/home/tphung/tphung_proj/Preprocessing_scRNA/data/4_Allen_Cortex.Hippocampus_Mouse_2019
Fri Jun  9 15:42:19 CEST 2023

wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/matrix.csv

# download the metadata
wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/metadata.csv
```

### Examine the data
```
adata
AnnData object with n_obs × n_vars = 73363 × 45768
```
- 73363 cells & 45768 genes

```
adata.X
array([[  0.,   0.,  32., ...,   0.,   0.,   0.],
       [  0.,   1.,  83., ...,   0.,   0.,   0.],
       [  0.,   0.,  37., ...,   0.,   0.,   0.],
       ...,
       [  0.,   0.,  50., ...,   0.,   0.,   0.],
       [  0.,   0.,  31., ...,   0.,   0.,   0.],
       [  0.,   0., 213., ...,   0.,   0.,   0.]], dtype=float32)
```
    + Since the values are integers, I am assuming that this is the raw count

```
adata.var
Empty DataFrame
Columns: []
Index: [0610005C13Rik, 0610006L08Rik, 0610007P14Rik, 0610009B22Rik, 0610009E02Rik, 0610009L18Rik, 0610009O20Rik, 0610010B08Rik, 0610010F05Rik, 0610010K14Rik, 0610011F06Rik, 0610012G03Rik, 0610025J13Rik, 0610030E20Rik, 0610031O16Rik, 0610037L13Rik, 0610038B21Rik, 0610039H22Rik, 0610039K10Rik, 0610040B10Rik, 0610040F04Rik, 0610040J01Rik, 0610043K17Rik, 1010001N08Rik, 1100001I12Rik, 1110001J03Rik, 1110002J07Rik, 1110002L01Rik, 1110002O04Rik, 1110004E09Rik, 1110004F10Rik, 1110006O24Rik, 1110007C09Rik, 1110008F13Rik, 1110008L16Rik, 1110008P14Rik, 1110012L19Rik, 1110015O18Rik, 1110017D15Rik, 1110018N20Rik, 1110019D14Rik, 1110020A21Rik, 1110025L11Rik, 1110028F11Rik, 1110028F18Rik, 1110032A03Rik, 1110032F04Rik, 1110034G24Rik, 1110036E04Rik, 1110037F02Rik, 1110038B12Rik, 1110038F14Rik, 1110046J04Rik, 1110051M20Rik, 1110054M08Rik, 1110057K04Rik, 1110057P08Rik, 1110058D11Rik, 1110058L19Rik, 1110059E24Rik, 1110059G10Rik, 1110065P20Rik, 1190002F15Rik, 1190002N15Rik, 1190003K10Rik, 1190005I06Rik, 1190007I07Rik, 1190028D05Rik, 1200014J11Rik, 1300002E11Rik, 1300017J02Rik, 1500002F19Rik, 1500004A13Rik, 1500009C09Rik, 1500009L16Rik, 1500011B03Rik, 1500011K16Rik, 1500012F01Rik, 1500012K07Rik, 1500015A07Rik, 1500015L24Rik, 1500015O10Rik, 1500017E21Rik, 1500026H17Rik, 1500035N22Rik, 1520401A03Rik, 1600002D24Rik, 1600002H07Rik, 1600002K03Rik, 1600010M07Rik, 1600012H06Rik, 1600014C10Rik, 1600014C23Rik, 1600014K23Rik, 1600015I10Rik, 1600019K03Rik, 1600020E01Rik, 1600023N17Rik, 1600025M17Rik, 1600027J07Rik, ...]

[45768 rows x 0 columns]
```
- Interestingly, the meta data has fewer cells (73347 cells as opposed to 73363 cells)
```
meta
               sample_name  donor_sex_id  ... cortical_layer_order cortical_layer_color
0        US-1250273_E1_S37             1  ...                   11              #7373FF
1        US-1250273_E2_S01             2  ...                   11              #7373FF
2        US-1250273_E2_S02             2  ...                   11              #7373FF
3        US-1250273_E2_S03             2  ...                   11              #7373FF
4        US-1250273_E2_S04             2  ...                   11              #7373FF
...                    ...           ...  ...                  ...                  ...
73342  SM-J39VT_S475_E1-50             1  ...                   17              #A7CC5C
73343  SM-J39VT_S476_E1-50             1  ...                   17              #A7CC5C
73344  SM-J39VT_S478_E1-50             1  ...                   17              #A7CC5C
73345  SM-J39VT_S479_E1-50             1  ...                   17              #A7CC5C
73346  SM-J39VT_S480_E1-50             1  ...                   17              #A7CC5C

[73347 rows x 57 columns]
```

- Subset
```
tissue_cell_ids = meta["sample_name"]
adata_subset = adata[tissue_cell_ids]
adata_subset
View of AnnData object with n_obs × n_vars = 73347 × 45768
```

- add meta
```
adata_subset.obs = meta
```

- Number of cell types
    + `class_label`: 3 cell types
    + `subclass_label`: 42 cell types
    + `cluster_label`: 382 cell types 

- Script for preprocessing: `qc_scrna_id4.py`
- Table recording the number of cells and genes
  
|Description       | Number of cells      | Number of genes       |
|  ---  |  ---  |  ---  |
| Original      | 73,347 cells      | 45,768 genes      |
| First filter      | 73,347 cells       | 42,532 genes       | 
| Second filter      | 73,347 cells       | 42,532 genes       | 
| Ngenes with ensembl      | 73,347 cells       | 16,141 genes       | 

## 5_Allen_Cortex.HippoCampus_Mouse_2020
### Download the data
```
pwd; date
/home/tphung/tphung_proj/Preprocessing_scRNA/data/5_Allen_Cortex.HippoCampus_Mouse_2020
Sun Jun 11 22:38:51 CEST 2023

wget https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/matrix.csv
--2023-06-11 22:38:55--  https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/matrix.csv
Resolving idk-etl-prod-download-bucket.s3.amazonaws.com (idk-etl-prod-download-bucket.s3.amazonaws.com)... 52.217.85.92, 52.216.178.67, 52.217.192.241, ...
Connecting to idk-etl-prod-download-bucket.s3.amazonaws.com (idk-etl-prod-download-bucket.s3.amazonaws.com)|52.217.85.92|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 72860251632 (68G) [application/vnd.ms-excel]
Saving to: ‘matrix.csv’

matrix.csv                              100%[============================================================================>]  67.86G  37.8MB/s    in 33m 37s
```
