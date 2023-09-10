# This script is used to process the h5ad files for use in magma cell typing
# This is written for the Siletti dataset but #TODO is to generalize this to other datasets
# Date: 2023-07-24
# Author: Tanya Phung (t.n.phung@vu.nl)

import scanpy as sc
import anndata
import pandas as pd
import argparse
import os
import sys


def main(args):
    # read in the scRNAseq data
    adata_original = anndata.read(args.h5ad)

    if adata_original.var.index.name == "Approved symbol":
        adata_original.var.index.name = "symbol"

    adata_original.var.reset_index(inplace=True)

    gene_list = adata_original.var[adata_original.var["ensembl"].notnull()].index.values.tolist() #this obtains a list of gene in symbol that has an ensemble conversion
    adata = adata_original[:, gene_list].copy() #subset based on the genes (keep the genes if it was converted to ensemble successfully)

    # now, adata.X is the matrix in numpy format where row is cell ID and column is the gene
    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} rows (cells) and {adata_ncol} columns (genes)")

    # # set the cell id to be the index of the pandas adata.obs
    # adata.obs.set_index('CellID', inplace=True) #NOTE THAT the key CellID might differ with other datasets. 
    adata.obs.index.name = None

    adata.var.set_index('symbol', inplace=True)
    adata.write(args.out_h5ad)

def parse_args():
    parser = argparse.ArgumentParser(description='Process h5ad to use as inputs for MAGMA Celltyping')
    parser.add_argument('--h5ad', required=True, help='Path to the h5ad data. This is the file after the "preprocessing" step.')
    parser.add_argument('--out_h5ad', required=True, help='Name of the output file for h5ad.')
    return parser.parse_args()


main(parse_args())

