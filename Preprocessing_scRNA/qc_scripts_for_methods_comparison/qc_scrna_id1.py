# This script outputs qc metrics for: 1_Allen_MCA_Human_2019
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2023-06-09
# Update on 2023-06-17 to incorporate Rachel's suggestions, namely: (1) rename cell type label to be more consistent, (2) save as sparse matrix
# Update on 2023-06-21: add 3 columns: keep_level_1, keep_level_2, keep_level_3 where TRUE mean keeping the cell and FALSE means removing the cell in downstream analyses because the number of cells for that cell type is less than 2

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from matplotlib import pyplot as plt
import os
import single_cell_helper_functions_v3
from scipy.sparse import csr_matrix

def main():
    # initialize
    id = "1_Allen_MCA_Human_2019" #change here per data
    base_dir = "/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data" #change here if necessary
    matrix_path = os.path.join(base_dir, id, "matrix.csv") #change here per data
    meta_path = os.path.join(base_dir, id, "metadata.csv") #change here per data
    gene_names_org = "symbol" #change here per data
    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    plot1 = os.path.join(base_dir, id, "plot1.jpeg")
    data_for_plot2_fn = os.path.join(base_dir, id, id + "_obs_for_plot2.csv")
    clean_adata_fp = os.path.join(base_dir, id, id + ".h5ad")

    # format table 2
    table2 = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2)

    # format table 3
    table3 = open(os.path.join(base_dir, id, "table3.csv"), "w")

    # read in adata
    adata = sc.read_csv(matrix_path)
    meta = pd.read_csv(meta_path)
    adata.obs = meta

    print("Viewing the adata observations.", file=log)
    print(adata.obs, file=log)

    print("Viewing the adata variables.", file=log)
    print(adata.var, file=log)

    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} cells and {adata_ncol} genes.", file=log)
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

    print(f"After first filtering (see documentation for definition), adata has {adata.n_obs} cells and {adata.n_vars} genes.", file=log)
    # save to table 2
    first_filter = ["Filter#1", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(first_filter), file=table2)

    # second filtering based on mt percentage
    adata = adata[adata.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata.n_obs} cells and {adata.n_vars} genes.", file=log)
    # save to table 2
    second_filter = ["Filter#2", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(second_filter), file=table2)

    # plot highest expressed genes
    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.highest_expr_genes(adata, n_top=20, show=False)
        plt.savefig(plot1, bbox_inches="tight")

    # rename cell type columns
    adata.obs.rename(columns={"class_label": "cell_type_level_1", "subclass_label": "cell_type_level_2", "cluster_label": "cell_type_level_3"}, inplace = True)

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset

    # create a list to store the cell types where the number of cells in that cell type is <2
    cts_level_1_rm = set({np.nan})
    cts_level_2_rm = set({np.nan})
    cts_level_3_rm = set({np.nan})

    cts_level1 = adata.obs["cell_type_level_1"].dropna().unique()
    l1_n_cells_with_annot = 0
    for ct in cts_level1:
        ct_data = adata[adata.obs["cell_type_level_1"] == ct, :].X
        out = ["Level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)
        l1_n_cells_with_annot += ct_data.shape[0]
        if ct_data.shape[0] < 2: #change here if necessary
            cts_level_1_rm.add(ct)
    l1_n_no_annot = ["Level_1", "unknown", str(adata.obs.shape[0]-l1_n_cells_with_annot)]
    print("|".join(l1_n_no_annot), file=table3)
    
    cts_level2 = adata.obs["cell_type_level_2"].dropna().unique()
    l2_n_cells_with_annot = 0
    for ct in cts_level2:
        ct_data = adata[adata.obs["cell_type_level_2"] == ct, :].X
        out = ["Level_2", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)
        l2_n_cells_with_annot += ct_data.shape[0]
        if ct_data.shape[0] < 2: #change here if necessary
            cts_level_2_rm.add(ct)
    l2_n_no_annot = ["Level_2", "unknown", str(adata.obs.shape[0]-l2_n_cells_with_annot)]
    print("|".join(l2_n_no_annot), file=table3)
    
    cts_level3 = adata.obs["cell_type_level_3"].dropna().unique()
    l3_n_cells_with_annot = 0
    for ct in cts_level3:
        ct_data = adata[adata.obs["cell_type_level_3"] == ct, :].X
        out = ["Level_3", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)
        l3_n_cells_with_annot += ct_data.shape[0]
        if ct_data.shape[0] < 2: #change here if necessary
            cts_level_3_rm.add(ct)
    l3_n_no_annot = ["Level_3", "unknown", str(adata.obs.shape[0]-l3_n_cells_with_annot)]
    print("|".join(l3_n_no_annot), file=table3)

    # add 3 columns to adata.obs (keep_level_1, keep_level_2, and keep_level_3)
    adata.obs['keep_level_1'] = [False if x in cts_level_1_rm else True for x in adata.obs['cell_type_level_1']]
    adata.obs['keep_level_2'] = [False if x in cts_level_2_rm else True for x in adata.obs['cell_type_level_2']]
    adata.obs['keep_level_3'] = [False if x in cts_level_3_rm else True for x in adata.obs['cell_type_level_3']]
    
    # convert to ensembl
    gene_names_fp = '/gpfs/work5/0/vusr0480/Preprocessing_scRNA/code/conversion_files/gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata.var.reset_index(inplace=True)
    adata.var = adata.var.rename(columns={'index': gene_names_org})
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log)
    adata_gene_converted.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2)

    # save 
    adata_gene_converted.X = csr_matrix(adata_gene_converted.X) #convert to sparse
    adata_gene_converted.write_h5ad(filename=clean_adata_fp)

    # save data to plot violin in R
    data_for_plot2 = adata_gene_converted.obs.reset_index()
    data_for_plot2.to_csv(data_for_plot2_fn, index=False)

    # print out the saved adata
    print(adata_gene_converted, file=log)
    print(adata_gene_converted.var, file=log) 
    print(adata_gene_converted.obs, file=log)
    print(adata_gene_converted.X, file=log)       

    table2.close()
    table3.close()
    log.close()


main()

