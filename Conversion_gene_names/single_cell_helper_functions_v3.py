#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:20:44 2022
Modified on Wed April 12 by Tanya Phung
Helper function for processing of single cell datasets
- add ensembl gene identifiers to a preprocessed dataset if other names are used
File with alternative human gene names gene_names_human.txt downloaded from genenames.org on 26/01/22
File with alternative human gene names gene_names_human_with_rat.txt downloaded from genenames.org on 11/4/22; now also including known RGD names
File with alternative mouse gene names gene_names_mouse.txt downloaded from http://www.informatics.jax.org/downloads/reports/index.html on 29/01/22
File with alternative rat gene names gene_names_rat.txt downloaded from https://download.rgd.mcw.edu/data_release/GENES.RAT.txt

Mouse - human homologs based on MGI ID from gene_names_human.txt
v2: rat conversion added
v3: Added previous symbols when matching symbols to other IDs ; added symbol as output name for human genes
@author: Rachel Brouwer
"""

import pandas as pd
import numpy as np
import re
import sys


def match_alias_previous(gene_names, merged, adata, alternate_symbol):
    """
    This function deals with alias and previous symbols.
    alternate_symbol is `Alias symbols` or `Previous symbols`
    """
    # create a list of all the alternate symbol:
    all_alternate_symbols = []
    for p in range(gene_names.shape[0]):
        all_alternate_symbols = all_alternate_symbols + list(re.split(', | ', gene_names.loc[p, alternate_symbol]))
    for gene in merged.index:
        # check first whether the gene name is in the all_alternate_symbols list; if so start further search
        if merged.loc[gene, "symbol"] != merged.loc[gene, 'Approved symbol']:
            if merged.loc[gene, 'symbol'] in all_alternate_symbols:
                print(merged.loc[gene, "symbol"])
                # alias is in the table, now we need to find out where
                found = []
                for p in range(gene_names.shape[0]):
                    if merged.loc[gene, 'symbol'] in list(re.split(', | ', gene_names.loc[p, alternate_symbol])):
                        found.append(p)
                if len(found) == 1:  # sometimes there are multiple other symbols that match, ignore these genes
                    print("found; approved symbol is:" + gene_names.loc[found[0], "Approved symbol"])
                    # now only add if this does not lead to duplicate entries:
                    if not gene_names.loc[found[0], 'NCBI Gene ID'] in merged['NCBI Gene ID'].tolist():
                        merged.iloc[gene, adata.var.shape[1]:] = gene_names.loc[found[0], :]
    return merged

def add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add):
    """ adata is an annotated data set with gene names (in some form) in the adata.var.
        adata must have a gene column (sometimes the gene names are stored as row index, one must convert this to a column before proceeding) and accepted gene column name are: {"Ensembl gene ID ; "HGNC ID"  ; "symbol" ; "NCBI Gene ID"}
        function add the ensembl gene names to adata.var
        current accepted names for gene_names_org = {"Ensembl gene ID ; "HGNC ID"  ; "symbol" ; "NCBI Gene ID"}
        current accepted names for gene_names_add = {"Ensembl gene ID ; "HGNC ID"  ; "NCBI Gene ID"}
        gene_names_fp is the path to the file used for conversion (i.e. conversion_files/gene_names_human.txt)
    """
    gene_names = pd.read_table(gene_names_fp)
    # remove Chromosome column (20230513 because the Siletti adata.var has a Chromosome column as well and this is creating 2 columns when merging)
    gene_names.drop(columns=['Chromosome'], inplace=True)
    print(f"INFO: The conversion file has {gene_names.shape[0]} rows and {gene_names.shape[1]} columns.")

    print(f"INFO: In the adata.var, there are {adata.var.shape[0]} rows and {adata.var.shape[1]} columns.")

    #check if gene_names_org exist in adata.var
    if gene_names_org not in adata.var.columns.tolist():
        sys.exit("value you specified for gene_names_org does not appear to exist in adata.var")

    conventions = ["Ensembl gene ID", "NCBI Gene ID", "HGNC ID"]
    if gene_names_org in conventions:
        print(f"INFO: User specifies that the gene name convention in the adata is {gene_names_org}.")
        gene_names = gene_names.dropna(subset=[gene_names_org])
        print(
            f"INFO: After removing NA in the gene column, the conversion file has {gene_names.shape[0]} rows and {gene_names.shape[1]} columns.")
        merged = pd.merge(adata.var, gene_names, how='left', on=gene_names_org)
        print(
            f"INFO: After merging adata.var and gene_names on {gene_names_org}, there are {merged.shape[0]} rows and {merged.shape[1]} columns.")
        merged.index = adata.var.index

    elif gene_names_org == "symbol": #we do symbol separately because of matching to Alias symbols and Previous symbols
        merged = pd.merge(adata.var, gene_names, how='left', left_on='symbol', right_on="Approved symbol")
        print(
            f"INFO: After merging adata.var and gene_names on {gene_names_org}, there are {merged.shape[0]} rows and {merged.shape[1]} columns.")

        # Check for Alias symbols - this could create duplicates if both the approved symbol and an alias is in the list, or if there are multiple aliases for different genes
        # replace nan by - in alias symbols
        gene_names['Alias symbols'] = gene_names['Alias symbols'].fillna("-")  # convert NaN to -
        merged_alias = match_alias_previous(gene_names, merged, adata, "Alias symbols")

        # we finally do the same for previous symbols
        gene_names['Previous symbols'] = gene_names['Previous symbols'].fillna("-")
        merged_alias_previous = match_alias_previous(gene_names, merged_alias, adata, "Previous symbols")

        merged_alias_previous.index = adata.var.index

    gene_key = {'Ensembl gene ID': 'ensembl', 'NCBI Gene ID': 'entrez_id', 'HGNC ID': 'hgnc'}
    if gene_names_add in conventions:
        merged_alias_previous.rename(columns={gene_names_add: gene_key[gene_names_add]}, inplace=True)
        outputvars = adata.var.columns.tolist() + [gene_key[gene_names_add]]
    elif gene_names_add == 'symbol':
        merged_alias_previous.rename(columns={"Approved symbol" : "symbol"}, inplace=True)
        outputvars = adata.var.columns.tolist() + [gene_key[gene_names_add]]

    adata.var = merged_alias_previous[outputvars]
    return adata


def add_gene_names_mouse(adata, gene_names_mouse_fp, gene_names_human_fp, gene_names_org, gene_names_add):
    """ adata is an annotated data set with gene names (in some form) in the adata.var
        first we need to add the MGI IDs (if not included) since these are used to match to human homologs
        current accepted names for gene_names_org = {"ENSMUSG" ; "symbol"}
        curren accepted names for gene_names_add = {"HGNC ID", "NCBI Gene ID", "Ensembl gene ID" }
    """
    print(f"INFO: In the adata.var, there are {adata.var.shape[0]} rows and {adata.var.shape[1]} columns.")
    gene_names_mouse = pd.read_table(gene_names_mouse_fp, header=None)
    gene_names_mouse.columns = ['MGI Marker Accession ID', 'Marker Symbol', 'Marker Name', 'cM Position',
                                'Chromosome', 'Ensembl Accession ID', 'Ensembl Transcript ID', 'Ensembl Protein ID', 'Feature Types', 'Genome Coordinate Start', 'Genome Coordinate End', 'Strand',
                                'BioTypes']
    print(f"INFO: The mouse conversion file has {gene_names_mouse.shape[0]} rows and {gene_names_mouse.shape[1]} columns.")
    # remove one gene that has two MGI and one ENSMUSG ID:
    gene_names_mouse = gene_names_mouse[gene_names_mouse['Ensembl Accession ID'] != 'ENSMUSG00000115016']
    print(f"INFO: After removing one gene that has teo MGI and one ENSMUSG ID, the conversion file has {gene_names_mouse.shape[0]} rows and {gene_names_mouse.shape[1]} columns.")

    if gene_names_org == "ENSMUSG":
        merged = pd.merge(adata.var, gene_names_mouse, how='left', left_on='ENSMUSG', right_on='Ensembl Accession ID')
        merged.index = adata.var.index

    elif gene_names_org == "symbol":
        # to avoid merging twice on genes that have an ambigous symbol remove all these rows from the merging data file
        gene_names_mouse = gene_names_mouse.loc[~gene_names_mouse['Marker Symbol'].duplicated(keep=False), :]
        print(
            f"INFO: After duplicate removal, the conversion file has {gene_names_mouse.shape[0]} rows and {gene_names_mouse.shape[1]} columns.")
        merged = pd.merge(adata.var, gene_names_mouse, how='left', left_on='symbol', right_on='Marker Symbol')
        merged.index = adata.var.index

    else:
        sys.exit("gene_names_org specified is not supported at this time. Please make sure that it's either ENSMUSG or symbol.")

    gene_names_human = pd.read_table(gene_names_human_fp)
    print(
        f"INFO: The human conversion file has {gene_names_human.shape[0]} rows and {gene_names_human.shape[1]} columns.")
    # remove duplicated MGI mouse IDs, we do not know which human gene they map to
    gene_names_human = gene_names_human.loc[~gene_names_human['Mouse genome database ID'].duplicated(keep=False)]
    print(
        f"INFO: After removing duplicated MGI mouse IDs, the human conversion file has {gene_names_human.shape[0]} rows and {gene_names_human.shape[1]} columns.")

    # now merge to the human homolog;
    merged = pd.merge(merged, gene_names_human[~gene_names_human['Mouse genome database ID'].isnull()], how='left',
                      left_on='MGI Marker Accession ID', right_on='Mouse genome database ID')

    gene_names_add_dict = {'NCBI Gene ID': 'entrez_id',
                           'HGNC ID': 'HGNC ID',
                           'Ensembl gene ID': 'ensembl'}

    merged.rename(columns={gene_names_add: gene_names_add_dict[gene_names_add]}, inplace=True)
    outputvars = adata.var.columns.tolist() + [gene_names_add_dict[gene_names_add]]
    # remove the genes that map to the same id
    duplicates = merged[merged[gene_names_add_dict[gene_names_add]].duplicated(keep=False)].index.values
    for gene in duplicates:
        merged.loc[gene, gene_names_add_dict[gene_names_add]] = np.nan
    adata.var = merged[outputvars]

    adata.var_rmna = adata.var.dropna(subset=gene_names_add_dict[gene_names_add])
    print(f"INFO: successfully converted {adata.var_rmna.shape[0]} gene.")
    return adata


def add_gene_names_rat(adata, gene_names_org, gene_names_add):
    """ adata is an annotated data set with gene names (in some form) in the adata.var
        first we need to add the RGD IDs (if not included) since these are used to match to human homologs
        function adds the entrez human gene names to adata.var
        current accepted names for gene_names_org = {"symbol"; "RGD ID" ; "Ensembl Rat ID"}
        curren accepted names for gene_names_add = {"HGNC ID", NCBI Gene ID"" }
    """

    gene_names_rat = pd.read_table('/project/prjsbrouwer2/single_cell/preprocessing/raw/conversion/gene_names_rat.txt',
                                   sep='\t', skiprows=84, low_memory=False)
    # this file contains a GENE_RGD_ID that we have to rewrite into RGD:# format
    gene_names_rat["RGD ID"] = "RGD:" + gene_names_rat['GENE_RGD_ID'].apply(str)
    #

    gene_names_human = pd.read_table(
        '/Users/rachel/Documents/single_cell/python_processing/raw/conversion/gene_names_human_incl_rat.txt',
        low_memory=False)
    # remove duplicated RGD rat IDs, we do not know which human gene they map to (this includes the NaNs)
    gene_names_human = gene_names_human.loc[
        ~gene_names_human['Rat genome database ID(supplied by RGD)'].duplicated(keep=False)]

    if gene_names_org == "Ensembl Rat ID":
        # first remove the rat data with duplicate Ensembl IDs (or missing, but that does not matter now as we won't be able to match)
        gene_names_rat = gene_names_rat[~gene_names_rat["ENSEMBL_ID"].duplicated(keep=False)]
        # now merge to the rat gene names to get the RGD values; then merge to the human data
        merged = pd.merge(adata.var, gene_names_rat, how='left', left_on='Ensembl Rat ID', right_on='ENSEMBL_ID')
        merged = pd.merge(merged,
                          gene_names_human[~gene_names_human['Rat genome database ID(supplied by RGD)'].isnull()],
                          how='left', left_on='RGD ID', right_on='Rat genome database ID(supplied by RGD)')
    if gene_names_org == "symbol":
        # first remove the rat data with duplicate symbols (or missing, but that does not matter now as we won't be able to match)
        gene_names_rat = gene_names_rat[~gene_names_rat["SYMBOL"].duplicated(keep=False)]
        merged = pd.merge(adata.var, gene_names_rat, how='left', left_on='symbol', right_on='SYMBOL')
        merged = pd.merge(merged,
                          gene_names_human[~gene_names_human['Rat genome database ID(supplied by RGD)'].isnull()],
                          how='left', left_on='RGD ID', right_on='Rat genome database ID(supplied by RGD)')

    # now merge to the human homolog;
    if gene_names_org == "RGD ID":
        merged = pd.merge(adata.var,
                          gene_names_human[~gene_names_human['Rat genome database ID(supplied by RGD)'].isnull()],
                          how='left', left_on='RGD ID', right_on='Rat genome database ID(supplied by RGD)')

    print(merged.columns)
    if gene_names_add == "NCBI Gene ID":
        merged.rename(columns={"NCBI Gene ID": 'entrez_id'}, inplace=True)
        print(merged.columns)
        outputvars = adata.var.columns.tolist() + ['entrez_id']
        print(outputvars)
        # remove the genes that map to the same entrez_id
        duplicates = merged[merged['entrez_id'].duplicated(keep=False)].index.values
        for gene in duplicates:
            merged.loc[gene, 'entrez_id'] = np.nan

    if gene_names_add == "HGNC ID":
        merged.rename(columns={"HGNC ID ", 'hgnc'}, inplace=True)
        outputvars = adata.var.columns.tolist() + ['hgnc']

    adata.var = merged[outputvars]

    return adata


def add_gene_names_human_Geschwind(adata):
    """ adata is an annotated data set with gene names (in some form) in the adata.var
        function add the ensembl gene names to adata.var
        This particular dataset provides gene names that are a combination of ENSEMBL IDs and Symbols
        Returns with added Entrez Gene IDs.
    """
    gene_names = pd.read_table('/project/prjsbrouwer2/single_cell/preprocessing/raw/conversion/gene_names_human.txt')

    adata.var['mixed_gene_name'] = adata.var.index
    merged = pd.merge(adata.var, gene_names, how='left', left_on='mixed_gene_name', right_on="Approved symbol")
    # we ignore previous symbols (probably previous for a reason)
    # but we look for Alias symbols
    # replace nan by - in alias symbols
    gene_names['Alias symbols'] = gene_names['Alias symbols'].fillna("-")
    # create a list of all the aliases:
    alias = []
    for p in range(gene_names.shape[0]):
        alias = alias + list(re.split(', | ', gene_names.loc[p, 'Alias symbols']))
    for gene in merged.index:
        # check first whether the gene name is in the aliases; if so start further search
        if merged.loc[gene, "mixed_gene_name"] != merged.loc[gene, 'Approved symbol']:
            if merged.loc[gene, 'mixed_gene_name'] in alias:
                print(merged.loc[gene, "mixed_gene_name"])
                # alias is in the table, now we need to find out where
                found = []
                for p in range(gene_names.shape[0]):
                    if merged.loc[gene, 'mixed_gene_name'] in list(
                            re.split(', | ', gene_names.loc[p, 'Alias symbols'])):
                        found.append(p)
                if len(found) == 1:  # sometimes there are multiple aliasses that match, ignore these genes
                    print("found; approved symbol is:" + gene_names.loc[found[0], "Approved symbol"])
                    # now only add if this does not lead to duplicate entries:
                    if not gene_names.loc[found[0], 'NCBI Gene ID'] in merged['NCBI Gene ID'].tolist():
                        merged.loc[gene, adata.var.shape[1]:] = gene_names.loc[found[0], :]

    # there are still empty lines where an ensembl ID was used
    for gene in merged.index:
        if 'ENSG' in merged.loc[gene, 'mixed_gene_name']:
            # only add when present in the gene_names file
            if merged.loc[gene, 'mixed_gene_name'] in gene_names['Ensembl gene ID'].tolist():
                merged.loc[gene, adata.var.shape[1]:] = gene_names.loc[np.where(
                    gene_names['Ensembl gene ID'] == merged.loc[gene, 'mixed_gene_name'])[0][0], :]
    merged.index = adata.var.index

    merged.rename(columns={'NCBI Gene ID': 'entrez_id'}, inplace=True)
    outputvars = adata.var.columns.tolist() + ['entrez_id']

    adata.var = merged[outputvars]
    # remove the genes that map to the same entrez_id
    duplicates = adata.var[adata.var['entrez_id'].duplicated(keep=False)].index.values
    for gene in duplicates:
        adata.var.loc[gene, 'entrez_id'] = np.nan

    return adata


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def compare_two_gene_selections(magmafile1, magmafile2):
    with open(magmafile1) as magmafile1:
        lines_file1 = [line.strip().split() for line in magmafile1]

    with open(magmafile2) as magmafile2:
        lines_file2 = [line.strip().split() for line in magmafile2]

    # we assume that the cell types are the same in both files, but in case it is not, we create a list of all the unique cell types in either file
    cts1 = [item[0] for item in lines_file1]
    cts2 = [item[0] for item in lines_file2]
    cts = np.unique(cts1 + cts2)

    for ct in cts:
        ct_list1_index = np.where([item[0] == ct for item in lines_file1])
        ct_list2_index = np.where([item[0] == ct for item in lines_file2])
        if len(ct_list1_index[0]) == 1 & len(ct_list2_index[0]) == 1:
            ct_list1 = lines_file1[ct_list1_index[0][0]][1:]
            ct_list2 = lines_file2[ct_list2_index[0][0]][1:]
            print(ct, len(ct_list1), len(ct_list2), len(intersection(ct_list1, ct_list2)))

    def compare_cell_selections(magmafile1):

        with open(magmafile1) as magmafile1:
            lines_file1 = [line.strip().split() for line in magmafile1]

        # we assume that the cell types are the same in both files, but in case it is not, we create a list of all the unique cell types in either file
        for i1 in range(len(lines_file1)):
            cts1 = lines_file1[i1][0]
            for i2 in range(i1):
                cts2 = lines_file1[i2][0]
                ct_list1 = lines_file1[i1][1:]
                ct_list2 = lines_file1[i2][1:]
                print(cts1, cts2, len(ct_list1), len(ct_list2), len(intersection(ct_list1, ct_list2)))

    def intersection_list(magmafile1):

        with open(magmafile1) as magmafile1:
            lines_file1 = [line.strip().split() for line in magmafile1]
        inter = intersection(lines_file1[0][1:], lines_file1[1][1:])
        for i in range(2, len(lines_file1)):
            inter = intersection(inter, lines_file1[i][i:])
        return (inter)

