#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:20:44 2022
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

def add_gene_names_human(adata,gene_names_org,gene_names_add):
    """ adata is an annotated data set with gene names (in some form) in the adata.var
        function add the ensembl gene names to adata.var
        current accepted names for gene_names_org = {"Ensembl gene ID ; "HGNC ID"  ; "symbol" ; "NCBI Gene ID"}
        current accepted names for gene_names_out = {"Ensembl gene ID ; "HGNC ID"  ; "NCBI Gene ID"}                                           
    """
    gene_names = pd.read_table('/project/prjsbrouwer2/single_cell/preprocessing/raw/conversion/gene_names_human.txt')    
    
    if gene_names_org == "Ensembl gene ID":
        gene_names = gene_names.dropna(subset=['Ensembl gene ID'])
        merged = pd.merge(adata.var,gene_names, how='left',on='Ensembl gene ID')
        merged.index = adata.var.index
    
    if gene_names_org == "NCBI Gene ID":
        gene_names = gene_names.dropna(subset=['NCBI Gene ID'])
        merged = pd.merge(adata.var,gene_names, how='left',on='NCBI Gene ID')
        merged.index = adata.var.index
  
    if gene_names_org == "HGNC ID":
        merged = pd.merge(adata.var,gene_names, how='left',on='HGNC ID')
    
    if gene_names_org == "symbol":
        merged = pd.merge(adata.var,gene_names, how='left',left_on='symbol',right_on="Approved symbol")
        # we ignore previous symbols (probably previous for a reason)
        # but we look for Alias symbols - this could create duplicates if both the approved symbol and an alias is in the list, or if there are multiple aliases for different genes
        # replace nan by - in alias symbols
        gene_names['Alias symbols'] = gene_names['Alias symbols'].fillna("-")
        # create a list of all the aliases:
        alias = []
        for p in range(gene_names.shape[0]):
            alias = alias + list(re.split(', | ',gene_names.loc[p,'Alias symbols'] ))
        for gene in merged.index:
            # check first whether the gene name is in the aliases; if so start further search
            if merged.loc[gene,"symbol"]!= merged.loc[gene,'Approved symbol']:
                if merged.loc[gene,'symbol'] in alias:
                    print(merged.loc[gene,"symbol"])
                    # alias is in the table, now we need to find out where
                    found = []
                    for p in range(gene_names.shape[0]):
                        if merged.loc[gene,'symbol'] in list(re.split(', | ',gene_names.loc[p,'Alias symbols'] )):
                            found.append(p)
                    if len(found) == 1: # sometimes there are multiple aliasses that match, ignore these genes
                        print("found; approved symbol is:" + gene_names.loc[found[0],"Approved symbol"])
                        # now only add if this does not lead to duplicate entries:
                        if not gene_names.loc[found[0],'NCBI Gene ID'] in merged['NCBI Gene ID'].tolist():
                            merged.loc[gene,adata.var.shape[1]:] = gene_names.loc[found[0],:]
                        
        # we finally do the same for previous symbols
                
        gene_names['Previous symbols'] = gene_names['Previous symbols'].fillna("-")           
        prev = []
        for p in range(gene_names.shape[0]):
            prev = prev + list(re.split(', | ',gene_names.loc[p,'Previous symbols'] ))
        for gene in merged.index:
            # check first whether the gene name is in the previous symbols; if so start further search
            if merged.loc[gene,"symbol"]!= merged.loc[gene,'Approved symbol']:
                if merged.loc[gene,'symbol'] in prev:
                    print(merged.loc[gene,"symbol"])
                    # alias is in the table, now we need to find out where
                    found = []
                    for p in range(gene_names.shape[0]):
                        if merged.loc[gene,'symbol'] in list(re.split(', | ',gene_names.loc[p,'Previous symbols'] )):
                            found.append(p)
                    if len(found) == 1: # sometimes there are multiple aliasses that match, ignore these genes
                        print("found; approved symbol is:" + gene_names.loc[found[0],"Approved symbol"])
                        # now only add if this does not lead to duplicate entries:
                        if not gene_names.loc[found[0],'NCBI Gene ID'] in merged['NCBI Gene ID'].tolist():
                            merged.loc[gene,adata.var.shape[1]:] = gene_names.loc[found[0],:]
                        
        merged.index = adata.var.index

    if gene_names_add == "Ensembl gene ID":
        merged.rename(columns={'Ensembl gene ID' : 'ensembl'},inplace=True)
        outputvars = adata.var.columns.tolist() + ['ensembl']
    if gene_names_add == "NCBI Gene ID":
        merged.rename(columns={'NCBI Gene ID' : 'entrez_id'},inplace=True)
        outputvars = adata.var.columns.tolist() + ['entrez_id']
    if gene_names_add == "HGNC ID":
        merged.rename(columns={'HGNC ID' : 'hgnc'},inplace=True)
        outputvars = adata.var.columns.tolist() + ['hgnc']
    if gene_names_add == "symbol":
        merged.rename(columns={"Approved symbol" : "symbol"},inplace=True)
        outputvars = adata.var.columns.tolist() + ['symbol']
    
    
    adata.var = merged[outputvars]
    return adata



def add_gene_names_mouse(adata,gene_names_org,gene_names_add):
    """ adata is an annotated data set with gene names (in some form) in the adata.var
        first we need to add the MGI IDs (if not included) since these are used to match to human homologs
        current accepted names for gene_names_org = {"ENSMUSG" ; "MGI" ; "symbol"}
        curren accepted names for gene_names_add = {"HGNC ID", "NCBI Gene ID", "Ensembl gene ID" }
    """
    gene_names_mouse = pd.read_table('/project/prjsbrouwer2/single_cell/preprocessing/raw/conversion/gene_names_mouse.txt',header=None)    
    gene_names_mouse.columns = ['MGI Marker Accession ID','Marker Symbol','Marker Name	','cM Position','Chromosome','Ensembl Accession ID','Ensembl Transcript ID','Ensembl Protein ID','Feature Types','Genome Coordinate Start','Genome Coordinate End','Strand','BioTypes']
    # remove one gene that has two MGI and one ENSMUSG ID:
    gene_names_mouse = gene_names_mouse[gene_names_mouse['Ensembl Accession ID'] != 'ENSMUSG00000115016']
     
    if gene_names_org == "ENSMUSG":
        merged = pd.merge(adata.var,gene_names_mouse, how='left',left_on='ENSMUSG',right_on='Ensembl Accession ID')
        merged.index = adata.var.index
  
    if gene_names_org == "symbol":
        # to avoid merging twice on genes that have an ambigous symbol remove all these rows from the merging data file
        gene_names_mouse = gene_names_mouse.loc[~gene_names_mouse['Marker Symbol'].duplicated(keep=False),:]
        merged = pd.merge(adata.var,gene_names_mouse, how='left',left_on='symbol',right_on='Marker Symbol')
        merged.index = adata.var.index
     
    gene_names_human = pd.read_table('/project/prjsbrouwer2/single_cell/preprocessing/raw/conversion/gene_names_human.txt')    
    # remove duplicated MGI mouse IDs, we do not know which human gene they map to
    gene_names_human = gene_names_human.loc[~gene_names_human['Mouse genome database ID'].duplicated(keep=False)]
    # now merge to the human homolog; 
    merged = pd.merge(merged,gene_names_human[~gene_names_human['Mouse genome database ID'].isnull()],how='left',left_on='MGI Marker Accession ID', right_on='Mouse genome database ID')

    
    print(merged.columns)
    if gene_names_add == "NCBI Gene ID":
        merged.rename(columns={"NCBI Gene ID" : 'entrez_id'},inplace=True)    
        print(merged.columns)
        outputvars = adata.var.columns.tolist() + ['entrez_id']
        print(outputvars)
        # remove the genes that map to the same entrez_id
        duplicates = merged[merged['entrez_id'].duplicated(keep=False)].index.values
        for gene in duplicates:
           merged.loc[gene,'entrez_id'] = np.nan
            
    if gene_names_add == "HGNC ID":
        #merged.rename(columns={'hgnc' ,"HGNC ID "},inplace =True)
        outputvars = adata.var.columns.tolist() + ["HGNC ID"]
        duplicates = merged[merged["HGNC ID"].duplicated(keep=False)].index.values
        for gene in duplicates:
           merged.loc[gene,"HGNC ID"] = np.nan
    
    if gene_names_add == "Ensembl gene ID":
        merged.rename(columns={'Ensembl gene ID' :'ensembl'},inplace =True)
        outputvars = adata.var.columns.tolist() + ["ensembl"]
        duplicates = merged[merged['ensembl'].duplicated(keep=False)].index.values
        for gene in duplicates:
           merged.loc[gene,'ensembl'] = np.nan
        
    adata.var = merged[outputvars]
    
    return adata


def add_gene_names_rat(adata, gene_names_org, gene_names_add):
    """ adata is an annotated data set with gene names (in some form) in the adata.var
        first we need to add the RGD IDs (if not included) since these are used to match to human homologs
        function adds the entrez human gene names to adata.var
        current accepted names for gene_names_org = {"symbol"; "RGD ID" ; "Ensembl Rat ID"} 
        curren accepted names for gene_names_add = {"HGNC ID", NCBI Gene ID"" }
    """
    
    gene_names_rat = pd.read_table('/project/prjsbrouwer2/single_cell/preprocessing/raw/conversion/gene_names_rat.txt',sep='\t',skiprows=84,low_memory=False)
    # this file contains a GENE_RGD_ID that we have to rewrite into RGD:# format
    gene_names_rat["RGD ID"] = "RGD:" + gene_names_rat['GENE_RGD_ID'].apply(str)
    # 
    
    
    gene_names_human = pd.read_table('/Users/rachel/Documents/single_cell/python_processing/raw/conversion/gene_names_human_incl_rat.txt',low_memory=False)    
    # remove duplicated RGD rat IDs, we do not know which human gene they map to (this includes the NaNs)
    gene_names_human = gene_names_human.loc[~gene_names_human['Rat genome database ID(supplied by RGD)'].duplicated(keep=False)]

    
    if gene_names_org == "Ensembl Rat ID":         
        # first remove the rat data with duplicate Ensembl IDs (or missing, but that does not matter now as we won't be able to match)
        gene_names_rat = gene_names_rat[~gene_names_rat["ENSEMBL_ID"].duplicated(keep=False)]
        # now merge to the rat gene names to get the RGD values; then merge to the human data
        merged = pd.merge(adata.var,gene_names_rat,how='left',left_on='Ensembl Rat ID', right_on='ENSEMBL_ID')
        merged = pd.merge(merged, gene_names_human[~gene_names_human['Rat genome database ID(supplied by RGD)'].isnull()],how='left',left_on='RGD ID', right_on='Rat genome database ID(supplied by RGD)')
    if gene_names_org == "symbol":
        # first remove the rat data with duplicate symbols (or missing, but that does not matter now as we won't be able to match)
        gene_names_rat = gene_names_rat[~gene_names_rat["SYMBOL"].duplicated(keep=False)]
        merged = pd.merge(adata.var,gene_names_rat,how='left',left_on='symbol', right_on='SYMBOL')
        merged = pd.merge(merged, gene_names_human[~gene_names_human['Rat genome database ID(supplied by RGD)'].isnull()],how='left',left_on='RGD ID', right_on='Rat genome database ID(supplied by RGD)')

    # now merge to the human homolog; 
    if gene_names_org == "RGD ID":
        merged = pd.merge(adata.var,gene_names_human[~gene_names_human['Rat genome database ID(supplied by RGD)'].isnull()],how='left',left_on='RGD ID', right_on='Rat genome database ID(supplied by RGD)')

    print(merged.columns)
    if gene_names_add == "NCBI Gene ID":
        merged.rename(columns={"NCBI Gene ID" : 'entrez_id'},inplace=True)    
        print(merged.columns)
        outputvars = adata.var.columns.tolist() + ['entrez_id']
        print(outputvars)
        # remove the genes that map to the same entrez_id
        duplicates = merged[merged['entrez_id'].duplicated(keep=False)].index.values
        for gene in duplicates:
           merged.loc[gene,'entrez_id'] = np.nan
            
    if gene_names_add == "HGNC ID":
        merged.rename(columns={"HGNC ID ",'hgnc' },inplace =True)
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
    merged = pd.merge(adata.var,gene_names, how='left',left_on='mixed_gene_name',right_on="Approved symbol")
    # we ignore previous symbols (probably previous for a reason)
    # but we look for Alias symbols
    # replace nan by - in alias symbols
    gene_names['Alias symbols'] = gene_names['Alias symbols'].fillna("-")
    # create a list of all the aliases:
    alias = []
    for p in range(gene_names.shape[0]):
        alias = alias + list(re.split(', | ',gene_names.loc[p,'Alias symbols'] ))
    for gene in merged.index:
        # check first whether the gene name is in the aliases; if so start further search
        if merged.loc[gene,"mixed_gene_name"]!= merged.loc[gene,'Approved symbol']:
            if merged.loc[gene,'mixed_gene_name'] in alias:
                print(merged.loc[gene,"mixed_gene_name"])
                # alias is in the table, now we need to find out where
                found = []
                for p in range(gene_names.shape[0]):
                    if merged.loc[gene,'mixed_gene_name'] in list(re.split(', | ',gene_names.loc[p,'Alias symbols'] )):
                        found.append(p)
                if len(found) == 1: # sometimes there are multiple aliasses that match, ignore these genes
                    print("found; approved symbol is:" + gene_names.loc[found[0],"Approved symbol"])
                    # now only add if this does not lead to duplicate entries:
                    if not gene_names.loc[found[0],'NCBI Gene ID'] in merged['NCBI Gene ID'].tolist():
                        merged.loc[gene,adata.var.shape[1]:] = gene_names.loc[found[0],:]
                    
    # there are still empty lines where an ensembl ID was used  
    for gene in merged.index:
        if 'ENSG' in merged.loc[gene,'mixed_gene_name']:
            # only add when present in the gene_names file
            if merged.loc[gene,'mixed_gene_name'] in gene_names['Ensembl gene ID'].tolist():
                merged.loc[gene,adata.var.shape[1]:] = gene_names.loc[np.where(gene_names['Ensembl gene ID'] == merged.loc[gene,'mixed_gene_name'])[0][0],:]         
    merged.index = adata.var.index

 
    merged.rename(columns={'NCBI Gene ID' : 'entrez_id'},inplace=True)
    outputvars = adata.var.columns.tolist() + ['entrez_id']

   
    adata.var = merged[outputvars]
    # remove the genes that map to the same entrez_id
    duplicates = adata.var[adata.var['entrez_id'].duplicated(keep=False)].index.values
    for gene in duplicates:
        adata.var.loc[gene,'entrez_id'] = np.nan

    return adata

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
    
def compare_two_gene_selections(magmafile1,magmafile2):
    
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
            print(ct,len(ct_list1),len(ct_list2),len(intersection(ct_list1,ct_list2)))
            
        
    def compare_cell_selections(magmafile1):
         
         with open(magmafile1) as magmafile1:
             lines_file1 = [line.strip().split() for line in magmafile1]

         # we assume that the cell types are the same in both files, but in case it is not, we create a list of all the unique cell types in either file
         for i1 in range(len(lines_file1)):
             cts1=lines_file1[i1][0]
             for i2 in range(i1):
                 cts2=lines_file1[i2][0]
                 ct_list1 = lines_file1[i1][1:]
                 ct_list2 = lines_file1[i2][1:]
                 print(cts1,cts2,len(ct_list1),len(ct_list2),len(intersection(ct_list1,ct_list2)))   
                 
    def intersection_list(magmafile1):
        
        with open(magmafile1) as magmafile1:
             lines_file1 = [line.strip().split() for line in magmafile1]
        inter = intersection(lines_file1[0][1:],lines_file1[1][1:])
        for i in range(2,len(lines_file1)):
            inter = intersection(inter,lines_file1[i][i:])
        return(inter)

