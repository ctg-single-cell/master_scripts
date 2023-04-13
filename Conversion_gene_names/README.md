# Gene names reference files

This directory contains the reference files that we use to go from one gene naming system to another. We chose HGNC (Human Genetics Nomenclature Committee) for human gene names, MGI (Mouse Genome Informatics) for mouse gene names and RGD (Rat Genome Database) for rat gene names as our primary resource. Since not all publicly available datasets use those gene naming conventions, it could be necessary to go from the original gene names to MGI/RGD, and then to human gene names. 

* gene_names_human.txt: For human gene names, we follow HGNC; file downloaded from (www.genenames.org) on 26/01/22. This file contains a multitude of alternative gene naming systems, including HGNC, Entrez (NCBI), Ensembl, symbol, symbol aliases, previous symbol. In addition this file contains the MGI identifier. This identifier is used as the link between human and mouse genes. We include only those that have a 1-1 correspondence. 
* gene_names_human_with_rat.txt: Same file as above, now also including RGD (Rat Genome Database) ID; file downloaded from (www.genenames.org) on 11/4/22; This ID is used as the link between human and rate genes. We include only those that have a 1-1 correspondence. 
* gene_names_mouse.txt: File containing alternative IDs for mouse genes; downloaded from http://www.informatics.jax.org/downloads/reports/index.html on 29/01/22.
* gene_names_rat.txt:  File containing alternative IDs for rat genes; downloaded from https://download.rgd.mcw.edu/data_release/GENES.RAT.txt on 3/5/22. 


