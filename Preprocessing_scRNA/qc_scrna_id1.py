# This script outputs qc metrics for id1
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2023-04-17

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from matplotlib import pyplot as plt
import os
import single_cell_helper_functions_v3

def main():
    id = "id1" #change here per data
    ct_colname = "subclass_label" #change here per data
    base_dir = "/project/prjstphung/Preprocessing_scRNA/data" #change here if necessary

    matrix_path = os.path.join(base_dir, id, "matrix.csv") #change here per data
    meta_path = os.path.join(base_dir, id, "metadata.csv") #change here per data

    gene_names_org = "symbol"

    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    plot1 = os.path.join(base_dir, id, "plot1.jpeg")
    plot2 = os.path.join(base_dir, id, "plot2.jpeg")

    clean_adata_fp = os.path.join(base_dir, id, "id1.h5ad")

    # format table 2
    table2 = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2)

    # format table 3
    table3 = open(os.path.join(base_dir, id, "table3.csv"), "w")

    adata = sc.read_csv(matrix_path)
    meta = pd.read_csv(meta_path)
    adata.obs = meta

    print("Viewing the adata observations.", file=log)
    print(adata.obs, file=log)

    print("Viewing the adata variables.", file=log)
    print(adata.var, file=log)

    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} rows and {adata_ncol} columns.", file=log)
    original = ["Original", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(original), file=table2)

    # mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata.n_obs} rows and {adata.n_vars} columns.", file=log)
    # save to table 2
    first_filter = ["Filter#1", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(first_filter), file=table2)

    # second filtering based on mt percentage
    adata = adata[adata.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata.n_obs} rows and {adata.n_vars} columns.", file=log)
    # save to table 2
    second_filter = ["Filter#2", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(second_filter), file=table2)

    # plot highest expressed genes
    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.highest_expr_genes(adata, n_top=20, show=False)
        plt.savefig(plot1, bbox_inches="tight")

    # plot violin
    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, groupby=ct_colname, rotation=45, show=False)
        plt.savefig(plot2, bbox_inches="tight")

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    cts = adata.obs[ct_colname].dropna().unique()
    for ct in cts:
        ct_data = adata[adata.obs[ct_colname] == ct, :].X
        out = [ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)
    

        # convert to ensembl if gene_name_org is not ensembl
    if gene_names_org != "ensembl":
        gene_names_fp = 'code/conversion_files/gene_names_human.txt' #path to gene_names_human.txt
        gene_names_add = "Ensembl gene ID"
        adata.var.reset_index(inplace=True)
        adata.var = adata.var.rename(columns={'index': gene_names_org})
        adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add)
        ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
        print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log)
        # save to table 2
        gene_conversion = ["Ensembl_converted", str(adata_gene_converted.n_obs), str(ngenes_after_conversion)]
        print(",".join(gene_conversion), file=table2)

        # save 
        adata_gene_converted.write_h5ad(filename=clean_adata_fp)
        # print out the saved adata
        print(adata, file=log)
        print(adata.var, file=log) 
        print(adata.obs, file=log)       

    table2.close()
    table3.close()
    log.close()


main()

