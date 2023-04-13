# Documenting how to test the functions in `single_cell_helper_functions_v3.py`
- Notes on 2023-04-13: refactoring has been done only on 2 functions (`add_gene_names_human` and `add_gene_names_mouse`)
## add_gene_names_human function
1. Download human data
- Download
```
wget https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad
```
2. Testing to make sure that the original version of the function `add_gene_names_mouse` produces the same reults as the refactored version.
- human original version
```commandline
# need to make sure to replace the hard-coded path on lines 28
adata = anndata.read("Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
# gene_names_org = "symbol"
# gene_names_add = "Ensembl gene ID"
adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata,gene_names_org,gene_names_add)
adata_gene_converted.var.to_csv("original.out.csv", index=False, sep=",")
```
- human refactored version
```commandline
adata = anndata.read("Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
gene_names_fp = 'gene_names_human.txt'
gene_names_org = "symbol"
gene_names_add = "Ensembl gene ID"
adata.var.reset_index(inplace=True)
adata.var = adata.var.rename(columns={'index': gene_names_org})
adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add)
adata_gene_converted.var.to_csv("refactor.out.csv", index=False, sep=",")
```
- Results: the file `original.out.csv` and `refactor.out.csv` are the same. 

## add_gene_names_mouse function
1. Download mouse data
- Download the zip file for Primary Motor Area (MOp) from: http://portal.brain-map.org/atlases-and-data/rnaseq/mouse-aca-and-mop-smart-seq
- Downloaded on 2023-04-13
- Unzip
2. Format the mouse data into adata file format
- adata.var is the genes
- for testing, I'm just using the exon and not introns
```
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
adata = sc.read_csv("mouse_MOp_cells_2018-10-04_exon-matrix.csv")
adata.X.shape
#(45768, 4916)
#here, the rows are genes and the columns are cells. for consistency, we tranpose so that the rows are cells and the columns are genes
adata = adata.transpose()
adata.X.shape
#(4916, 45768)

# load in the gene name as adata.var & update the gene column to symbol to be compatible with the script singe_cell_helper_functions_v3.py
gene = pd.read_csv("mouse_MOp_cells_2018-10-04_genes-rows.csv")
gene.head()
            gene chromosome  entrez_id                   gene_name
0  0610005C13Rik          7      71661  RIKEN cDNA 0610005C13 gene
1  0610006L08Rik          7      76253  RIKEN cDNA 0610006L08 gene
2  0610007P14Rik         12      58520  RIKEN cDNA 0610007P14 gene
3  0610009B22Rik         11      66050  RIKEN cDNA 0610009B22 gene
4  0610009E02Rik          2  100125929  RIKEN cDNA 0610009E02 gene

gene.shape
(45768, 4)

gene=gene.rename(columns = {'gene':'symbol'})
gene.head()
          symbol chromosome  entrez_id                   gene_name
0  0610005C13Rik          7      71661  RIKEN cDNA 0610005C13 gene
1  0610006L08Rik          7      76253  RIKEN cDNA 0610006L08 gene
2  0610007P14Rik         12      58520  RIKEN cDNA 0610007P14 gene
3  0610009B22Rik         11      66050  RIKEN cDNA 0610009B22 gene
4  0610009E02Rik          2  100125929  RIKEN cDNA 0610009E02 gene

adata.var = gene

#format adata.obs
cells = pd.read_csv("mouse_MOp_cells_2018-10-04_samples-columns.csv")
adata.obs = cells

# save adata as h5ad
adata.write("mouse_MOp_cells_2018-10-04_exon.h5ad")
```

3. Testing to make sure that the original version of the function `add_gene_names_mouse` produces the same reults as the refactored version.
- mouse original version
```commandline
# need to make sure to replace the hard-coded path on lines 118 and 133
adata = anndata.read("mouse_MOp_cells_2018-10-04_exon.h5ad")
gene_names_org = "symbol"
gene_names_add = "Ensembl gene ID"
adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_mouse(adata,gene_names_org,gene_names_add)
adata_gene_converted.var.to_csv("test_original.csv", index=False, sep=",")
```

- mouse refactored version
```commandline
adata = anndata.read("mouse_MOp_cells_2018-10-04_exon.h5ad")
gene_names_mouse_fp = "gene_names_mouse.txt"
gene_names_human_fp = "gene_names_human.txt"
gene_names_org = "symbol"
gene_names_add = "Ensembl gene ID"
adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_mouse(adata, gene_names_mouse_fp, gene_names_human_fp, gene_names_org, gene_names_add)
adata_gene_converted.var.to_csv("test_refactored.csv", index=False, sep=",")
```
- Results: the file `test_original.csv` and `test_refactored.csv` are the same. 