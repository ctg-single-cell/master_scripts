library(EWCE)
library(anndata)
suppressPackageStartupMessages(library("argparse"))

set.seed(1234)

# create parser object
parser = ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--h5ad", help ="Path to the h5ad file.")
parser$add_argument("--outdir", help ="Path to the output file")
parser$add_argument("--file_prefix", help ="The prefix of the output file") #for example: 1_Allen_MCA_Human_2019
parser$add_argument("--groupName", help ="Name of group") #for example: all_entrez
args = parser$parse_args()

# get arguments
h5ad = args$h5ad
outdir = args$outdir
file_prefix = args$file_prefix
groupName = args$groupName

# get count from anndata
adata = read_h5ad(h5ad)

# Turn counts into matrix
# Because as.matrix(adata$X) returns a matrix where rows are the cells and columns are the genes. Therefore, we need to transpose it to have rows are the cells and columns are the genes.
count <- t(as.matrix(adata$X))
# print(ncol(count))
# print(nrow(count))

annotLevels = list(level1class=adata$obs$supercluster_term,
                   level2class=adata$obs$cell_type)

specificity_out <- EWCE::generate_celltype_data(
  exp = count,
  annotLevels = annotLevels,
  groupName = groupName,
  savePath = outdir,
  file_prefix = file_prefix,
  input_species = "human",
  dendrograms = FALSE
)