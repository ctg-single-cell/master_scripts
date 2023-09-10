# This script runs cellex in order to prepare the input for cellect
# siletti et al. 
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2023-07-22

import scanpy as sc
import anndata
import pandas as pd
import argparse
import os
import sys
import numpy as np
import cellex
import scipy.sparse


def main(args):
    # read in the scRNAseq data
    adata_original = anndata.read(args.h5ad)
    ct_column = "supercluster_term" 

    # because the scRNAseq data contains genes that don't have ensemble, it might be cleaner to create a filtered adata where:
    # genes that were not convertible to ensembl
    adata_original.var_names_make_unique()
    gene_list = adata_original.var[adata_original.var["ensembl"].notnull()].index.values.tolist() #this obtains a list of gene in symbol that has an ensemble conversion
    adata = adata_original[:, gene_list].copy() #subset based on the genes (keep the genes if it was converted to ensemble successfully)

    # adata.X is the matrix in numpy format where row is cell ID and column is the gene
    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} rows (cells) and {adata_ncol} columns (genes)")


    print("Subset the metadata to have 1 column (cell_type) & rowname is the name of the cells")
    adata.obs["sample_name"] = adata.obs.index
    cell_colname = "sample_name" #example: sample_name
    metadata_temp = adata.obs[[cell_colname, ct_column]]
    print("Inspecting metadata dataframe after subsetting to 2 columns")
    print(metadata_temp.head())
    # convert cell_colname to index
    metadata = metadata_temp.set_index(cell_colname)
    print("Inspecting metadata dataframe after converting a column to index")
    print(metadata.head())
    metadata.index.name = None #removing the index name
    print("Inspecting final metadata dataframe")
    print(metadata.head())
    print(f"After modilyfing, metadata has {metadata.shape[0]} rows (cells) and {metadata.shape[1]} columns (attributes)")

    print("Transpose the matrix so that rows is genes and columns is cells")
    adata_t = adata.X.transpose()
    print(f"After transposing, adata has {adata_t.shape[0]} rows (genes) and {adata_t.shape[1]} columns (cells)")

    print("Making a dataframe for cellex input where rows is genes and columns is cells. Gene names are row index.")
    row_genes_name = np.asarray(adata.var["ensembl"])
    print(f"There are {len(row_genes_name)} genes.")
    col_cells_name = np.asarray(adata.obs['sample_name'])
    print(f"There are {len(col_cells_name)} cells.")
    cellex_data = pd.DataFrame.sparse.from_spmatrix(adata_t, index=row_genes_name, columns=col_cells_name)
    print("Inspect the cellex dataframe")
    print(cellex_data.iloc[0:5, 0:5])
    # cellex_data.to_csv(os.path.join(args.outdir, "cellex_data.csv"), index=True)
    # metadata.to_csv(os.path.join(args.outdir, "metadata.csv"), index=True)

    print("Now run cellex")
    eso = cellex.ESObject(data=cellex_data, annotation=metadata, verbose=True)
    eso.compute(verbose=True)
    prefixData = args.cellex_out_prefix
    dirOut = args.outdir
    eso.save_as_csv(file_prefix=prefixData, path=dirOut, verbose=True)

    # replace " " and "/" in the cell type name
    cellex_output = pd.read_csv(os.path.join(args.outdir, args.cellex_out_prefix + ".esmu.csv.gz"))
    cellex_output.columns = cellex_output.columns.str.replace(" ", "_")
    cellex_output.columns = cellex_output.columns.str.replace("/", "_")
    cellex_output.to_csv(os.path.join(args.outdir, args.cellex_out_prefix + ".esmu_fmt.csv"), index=False)


def parse_args():
    parser = argparse.ArgumentParser(description='Generate input for cellect using cellex and scRNAseq data')
    parser.add_argument('--h5ad', required=True, help='Path to the h5ad data. This is the file after the "preprocessing" step.')
    parser.add_argument('--cellex_out_prefix', required=True, help='Prefix of the output file')
    parser.add_argument('--outdir', required=True, help='Path to the output directory')
    return parser.parse_args()


main(parse_args())

