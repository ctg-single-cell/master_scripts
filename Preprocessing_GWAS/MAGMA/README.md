# Step 1: Preparing MAGMA on the input files

Set these directories on lisa to run:

```
magma=/home/brouwer2/sources/magma_1.10/magma
ref_input=/home/brouwer2/sources/magma_1.10/g1000_eur
phenos=/home/brouwer2/sources/sumstats/input_GWAS
gene_loc_dir=/home/brouwer2/sources/sumstats/magma_gene/gene_locations
```

See directory gene_locations for location files using different gene naming systems. 

```
gene_loc_all_entrez=${gene_loc_dir}/all_entrez_genes_37p13.loc 
gene_loc_protein_entrez=${gene_loc_dir}/protein_coding_entrez_genes_37p13.loc
gene_loc_protein_entrez_magma=${gene_loc_dir}/NCBI37.3.gene.loc # older version of gene list for comparison, distributed with MAGMA 
gene_loc_all_ensemble=${gene_loc_dir}/all_ensmbl_genes.loc
gene_loc_protein_ensemble=${gene_loc_dir}/ensmbl_protein_coding_genes.loc
```

The following needs to be done only once (creation of the annotation file)

```
$magma --annotate --snp-loc $ref_input/g1000_eur.bim --gene-loc $gene_loc_all_entrez --out ${gene_loc_dir}/all_entrez &
$magma --annotate --snp-loc $ref_input/g1000_eur.bim --gene-loc $gene_loc_protein_entrez --out ${gene_loc_dir}/protein_coding_entrez & 
$magma --annotate --snp-loc $ref_input/g1000_eur.bim --gene-loc $gene_loc_protein_entrez_magma --out ${gene_loc_dir}/protein_coding_entrez_magma &
$magma --annotate --snp-loc $ref_input/g1000_eur.bim --gene-loc $gene_loc_all_ensemble --out ${gene_loc_dir}/all_ensemble & 
$magma --annotate --snp-loc $ref_input/g1000_eur.bim --gene-loc $gene_loc_protein_ensemble --out ${gene_loc_dir}/protein_coding_ensemble &
```

Output of these scripts (.annot) files can be found in the gene_locations directory
--- 

Step 2: Creating gene-based files from GWAS summary statistics

The run_magma.do script contains the program call to MAGMA with GWAS specific settings. In an attempt not to rerun and make optimal use of resources, this script takes a "batch" parameter so that only part of the phenotypes can be run. 
Submit to the cluster using 
```
sbatch run_magma.do --export=batch="batch" 
```
where batch can be "first" or second"

*NOTE:* reformatting needs to be done for some of the files - see README in Preprocess_GWAS directory. i
*NOTE:* all MAGMA analyses were run using multi option - new release version 1.10 and are based on entrez IDs or ensemble IDs

