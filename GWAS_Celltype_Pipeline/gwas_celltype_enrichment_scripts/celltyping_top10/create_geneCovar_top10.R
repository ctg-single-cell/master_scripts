# This R script is to create geneCovar for magma cell typing
# for top 10 percent
# Notes on 2023-07-12: right now there is a script `create_geneCovar.R` to process for linear and `create_geneCovar_top10` for processing top10. Will need to combine into 1 script
# Change log:
# The original script is in /mnt/c/Users/tph205/Documents/ctg/codes/prioritize_ct_methods_comp/codes/magma_celltyping/create_geneCovar_top10.R but I copied it to this directory to process for siletti
library(orthogene)
library(dplyr)
library(EWCE)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser = ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--ctd", help ="Path to the ctd file. This is the output from script process_scrna_ewce.R")
parser$add_argument("--magma_gene_out", help ="Path to the magma gene out file (entrez)")
parser$add_argument("--out_basename", help ="Path to the basename of the output file.")
args = parser$parse_args()

# get arguments
ctd = args$ctd
magma = args$magma_gene_out
out_basename = args$out_basename

# source scripts from MAGMA.Celltyping
setwd("/gpfs/home6/tphung/software/MAGMA_Celltyping/R")

for (f in list.files(pattern="*.R")) {
  source(f)
}
for (f in list.files(pattern="*.r")) {
  source(f)
}

# get hgnc2entrez_orthogene
gene_map <- orthogene::all_genes(species = "human",
                                 method = "gprofiler",
                                 target = "ENTREZGENE_ACC",
                                 ensure_filter_nas = FALSE)
hgnc2entrez_orthogene <- gene_map |>
  dplyr::select(hgnc_symbol = Gene.Symbol,
                entrez = target) |>
  unique()

# load ctd file
ctd = EWCE::load_rdata(ctd)

ctd <- prepare_quantile_groups(ctd = ctd, input_species = "human", output_species = "human")

# initialize outputs
geneCovar_level1 = paste0(out_basename, "_level_1")

geneCovarFile_level1 <- create_top10percent_genesets_file(
  genesOutFile = magma,
  ctd = ctd,
  annotLevel = 1,
  ctd_species = "human",
  output_path = geneCovar_level1
)