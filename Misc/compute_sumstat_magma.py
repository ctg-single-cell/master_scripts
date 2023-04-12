# This script returns mean, sd, cov, varom, and spec
# Row is the gene
# Column is the cell type
# Credit: Modified from Rachel's original script single_cell_MAGMA_files_v2.py
# Example usage:
# python compute_sumstat_magma.py --adata {input.adata} --ct_colname {params.ct_colname} --outdir {params.outdir} --genes_conversion {input.gene_conversion}

import scanpy as sc
import pandas as pd
import argparse
import os
import sys

def read_h5ad(data_path):
    """data is in .h5ad file format"""
    return sc.read_h5ad(data_path)

def read_csv(matrix_path, meta_path):
    "matrix is the count matrix in csv and meta is the meta file in csv format"
    adata = sc.read_csv(matrix_path)
    meta = pd.read_csv(meta_path)
    adata.obs = meta
    return adata



def convert_gene_names(genes_conversion, ctmatrix):
    """
    This function is used to convert (for now) from gene symbol to ensemble id
    :return:
    """
    # for now this is making a lot of assumptions but will need to make this more uniform
    with open(genes_conversion, "r") as f:
        genes_dict = {line.rstrip("\n").split()[1]: line.rstrip("\n").split()[0] for line in f}

    ctmatrix_rename = ctmatrix.rename(genes_dict, axis='columns')
    ctmatrix_t = ctmatrix.transpose()
    ctmatrix_t_rename = ctmatrix_t.rename(genes_dict, axis='rows')
    cond = ctmatrix.columns.isin(ctmatrix_rename.columns)
    ctmatrix_t_rename_drop = ctmatrix_t_rename.drop(ctmatrix_t[cond].index, inplace=False)
    return ctmatrix_t_rename_drop

def main(args):
    """
    :param adata:
    :param ct_colname: this is the column name that define the cell type from the metadata. Typically this could be cell_type, cluster_label, etc...
    :param outdir
    :param genes_conversion: this is the path to the gene converion file. For example: /home/tphung/data/reference/ENSG.genes.txt
    """

    data_path = args.data
    if data_path.endswith("h5ad"):
        adata = read_h5ad(data_path)
    elif data_path.endswith("csv"):
        if not args.meta:
            sys.exit("The input data is in csv file; therefore, the meta data file is required")
        else:
            adata =  read_csv(data_path, args.meta)
    ct_colname = args.ct_colname
    outdir = args.outdir
    genes_conversion = args.genes_conversion

    # define the available cell types in the dataset
    cts = adata.obs[ct_colname].dropna().unique()

    # first normalize
    # normalize counts per cell (in place)
    sc.pp.normalize_total(adata, target_sum=1e6)

    # define the genes in the dataset
    genes = adata.var_names

    # before log-transforming, we determine a filtering on a minimum of 1 transcript per million per celltype
    # start with computing the mean cpM per cell type
    means_cell_counts_pM = pd.DataFrame(data=None, index=cts, columns=genes, dtype=float)  # rows is the cell type; columns is the gene
    if data_path.endswith("h5ad"):
        for ct in cts:
            means_cell_counts_pM.loc[ct, :] = adata[adata.obs[ct_colname] == ct, :].X.mean(0)
    elif data_path.endswith("csv"):
        for ct in cts:
            means_cell_counts_pM.loc[ct, :] = adata.X[adata.obs[ct_colname] == ct, :].mean(0)

    # the low_filter matrix contains zeros for genes that have < 1 count pM in a certain cell type
    # this matrix will be used to multiply the same data when it is log-transformed
    low_filter = (1 - (means_cell_counts_pM < 1))  # 0 is <1, 1 if >1

    # specificity values
    # set NaNs (i.e. no expression in any of the cells/celltypes) in the mean matrix per cell type to zero
    spec_cell_counts_pM = ((means_cell_counts_pM.mul(low_filter)) / (means_cell_counts_pM.mul(low_filter)).sum(axis=0)).fillna(0)

    # log-transform the original data
    sc.pp.log1p(adata.X, base=2)

    # compute a measure per cell type (mean/sd/cov)
    means_cell_log_counts_pM = pd.DataFrame(data=None, index=cts, columns=genes, dtype=float)
    sds_cell_log_counts_pM = pd.DataFrame(data=None, index=cts, columns=genes, dtype=float)

    for ct in cts:
        Y = adata[adata.obs[ct_colname] == ct, :].to_df()
        Y.columns = genes
        means_cell_log_counts_pM.loc[ct, :] = Y.mean(0)
        sds_cell_log_counts_pM.loc[ct, :] = Y.std(0)

    cov_cell_log_counts_pM = sds_cell_log_counts_pM / means_cell_log_counts_pM
    varom_cell_log_counts_pM = sds_cell_log_counts_pM * cov_cell_log_counts_pM

    # filter out the values that have low counts per cell type
    means_cell_log_counts_pM = means_cell_log_counts_pM.mul(low_filter).fillna(0)
    sds_cell_log_counts_pM = sds_cell_log_counts_pM.mul(low_filter).fillna(0)
    cov_cell_log_counts_pM = (cov_cell_log_counts_pM.mul(low_filter)).fillna(0)
    varom_cell_log_counts_pM = (varom_cell_log_counts_pM.mul(low_filter)).fillna(0)

    # convert from gene symbol to ensemble gene
    means_cell_log_counts_pM.loc["Average"] = (means_cell_log_counts_pM.mean(axis=0))
    means_cell_log_counts_pM.index = [w.replace(' ', '_') for w in means_cell_log_counts_pM.index.values]
    means_cell_log_counts_pM.index = [w.replace('/', '_') for w in means_cell_log_counts_pM.index.values]
    means_cell_log_counts_pM_convert = convert_gene_names(genes_conversion, means_cell_log_counts_pM)

    sds_cell_log_counts_pM.index = [w.replace(' ', '_') for w in sds_cell_log_counts_pM.index.values]
    sds_cell_log_counts_pM.index = [w.replace('/', '_') for w in sds_cell_log_counts_pM.index.values]
    sds_cell_log_counts_pM_convert = convert_gene_names(genes_conversion, sds_cell_log_counts_pM)

    cov_cell_log_counts_pM.index = [w.replace(' ', '_') for w in cov_cell_log_counts_pM.index.values]
    cov_cell_log_counts_pM.index = [w.replace('/', '_') for w in cov_cell_log_counts_pM.index.values]
    cov_cell_log_counts_pM_convert = convert_gene_names(genes_conversion, cov_cell_log_counts_pM)

    varom_cell_log_counts_pM.index = [w.replace(' ', '_') for w in varom_cell_log_counts_pM.index.values]
    varom_cell_log_counts_pM.index = [w.replace('/', '_') for w in varom_cell_log_counts_pM.index.values]
    varom_cell_log_counts_pM_convert = convert_gene_names(genes_conversion, varom_cell_log_counts_pM)

    spec_cell_counts_pM.index = [w.replace(' ', '_') for w in spec_cell_counts_pM.index.values]
    spec_cell_counts_pM.index = [w.replace('/', '_') for w in spec_cell_counts_pM.index.values]
    spec_cell_counts_pM_convert = convert_gene_names(genes_conversion, spec_cell_counts_pM)

    means_cell_log_counts_pM_convert.to_csv(os.path.join(outdir, "means_cell_log_counts_pM_convert.tsv"), sep="\t", index=True, index_label="GENE") # save to a file)
    sds_cell_log_counts_pM_convert.to_csv(os.path.join(outdir, "sds_cell_log_counts_pM_convert.tsv"), sep="\t",
                                            index=True, index_label="GENE")  # save to a file)
    cov_cell_log_counts_pM_convert.to_csv(os.path.join(outdir, "cov_cell_log_counts_pM_convert.tsv"), sep="\t",
                                            index=True, index_label="GENE")  # save to a file)
    varom_cell_log_counts_pM_convert.to_csv(os.path.join(outdir, "varom_cell_log_counts_pM_convert.tsv"), sep="\t",
                                            index=True, index_label="GENE")  # save to a file)
    spec_cell_counts_pM_convert.to_csv(os.path.join(outdir, "spec_cell_log_counts_pM_convert.tsv"), sep="\t",
                                            index=True, index_label="GENE")  # save to a file)

def parse_args():
    parser = argparse.ArgumentParser(description='Generate input for magma gene property from scRNAseq data')
    parser.add_argument('--data', required=True, help='Path to the data. Could be either in h5ad format or csv format')
    parser.add_argument('--meta', required=False, help='Path to the meta file. Only needed if the data is in csv file format')
    parser.add_argument('--ct_colname', required=False, help='Column label of the cell type')
    parser.add_argument('--outdir', required=False, help='Path to the output directory')
    parser.add_argument('--genes_conversion', required=False,
                        help='Path to the gene conversion file')
    return parser.parse_args()


main(parse_args())

