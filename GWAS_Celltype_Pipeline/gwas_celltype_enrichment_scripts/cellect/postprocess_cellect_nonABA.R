# This script processes outputs from cellect for each trait
# This is for the data that don't have levels
# Date: 2023-07-23

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--input", help ="Path to the input file.")
parser$add_argument("--trait", help ="Specify the trait.")
parser$add_argument("--scrnaseq_id", help ="Specify the id of the scrnaseq data.")
parser$add_argument("--outfile", help ="Specify the path to the output file.")

args = parser$parse_args()

# get arguments
input = args$input
trait = args$trait
scrnaseq_id = args$scrnaseq_id

cellect = read.csv(input)
cellect_trait = cellect %>% filter(gwas==trait)
cellect_trait_id = cellect_trait %>% 
  filter(specificity_id==scrnaseq_id)
cellect_trait_id$p_within = p.adjust(cellect_trait_id$pvalue, method = "bonferroni")
cellect_trait_id_sig = cellect_trait_id %>% filter(p_within<0.05)
cellect_trait_id_sig_clean = cellect_trait_id_sig %>% 
  select(annotation, p_within)
if (nrow(cellect_trait_id_sig_clean) > 0) {
cellect_trait_id_sig_clean$dataset = scrnaseq_id
cellect_trait_id_sig_clean$level = "level_1"
new_order = c("dataset", "level", "annotation", "p_within")
cellect_trait_id_sig_clean = cellect_trait_id_sig_clean[, new_order]
colnames(cellect_trait_id_sig_clean) = c("dataset", "level", "celltype", "padj")
} else {
  cellect_trait_id_sig_clean = data.frame(dataset=NA, level=NA, celltype=NA, padj=NA)
  }

write.table(cellect_trait_id_sig_clean, args$outfile, quote = F, row.names = F, sep = "\t")