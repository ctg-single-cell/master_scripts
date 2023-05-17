#!/usr/bin/env Rscript

# Date: 18-04-2023
# Author: Tanya Phung (t.n.phung@vu.nl)
# This script is used to compute adjusted pvalue after magma step 1
suppressPackageStartupMessages(library("argparse"))
library(data.table)
library(dplyr)
library(tidyr)

# example usage:
# Rscript calc_adj_pval_magma_step1.R --magma_outdir {/path/to/magma_outdir} # use bonferroni as a default
# Rscript calc_adj_pval_magma_step1.R --magma_outdir {/path/to/magma_outdir} --adjPmeth {correction_method} #Use other correction method other than bonferroni.

# create parser object
parser = ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--magma_outdir", help ="Path to the directory with results from running magma cell type step 1.")
parser$add_argument("--adjPmeth", default="bonferroni", help ="Adjusted method for multiple hypothesis testing. Default is bonferroni. \n
For accepted method, check ?p.adjust function in R.")

args = parser$parse_args()

filedir = args$magma_outdir
adjPmeth = args$adjPmeth

# Load all of the .gsa.out file for each trait that has been run against multiple scRNAseq dataset
files <- list.files(path=filedir, pattern="*.gsa.out.rmComment", recursive=FALSE)
print(files)

step1 <- data.frame() #combine all of the datasets
for(file in files){
  tmp = fread(paste0(filedir, file), data.table=F)
  ds = strsplit(file, "_magma.gsa.out")[[1]][1] #NOTE: this only work here when the file format is {nameofscRNA}_magma.gsa.out
  if("FULL_NAME" %in% colnames(tmp)){ #If variable FULL_NAME is in the dataframe, then convert VARIABLE to FULL_NAME
    tmp$VARIABLE <- tmp$FULL_NAME
    tmp <- tmp[,-ncol(tmp)]
  }
  tmp <- tmp[order(tmp$P),]
  tmp$ds <- ds
  tmp$P.adj.pds <- p.adjust(tmp$P, method=adjPmeth) #P.adj.pds is the adjusted p values for each dataset
  if(nrow(step1)==0){step1 <- tmp}
  else{step1 <- rbind(step1, tmp)}
}
step1$P.adj <- p.adjust(step1$P, method=adjPmeth)
tmp_out <- step1[,c("ds", "VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P", "P.adj.pds", "P.adj")]
colnames(tmp_out)[1:2] <- c("Dataset", "Cell_type")
write.table(tmp_out, paste0(filedir, "magma_celltype_step1.txt"), quote=F, row.names=F, sep="\t")

#return cell types that are significant within the same dataset of the same tissue
tmp_out_within = tmp_out %>%
  filter(P.adj.pds < 0.05)
write.table(tmp_out_within, paste0(filedir, "magma_celltype_step1_significant_within.txt"), quote=F, row.names=F, sep="\t")

#return cell types that are significant across datasets of the same tissue
tmp_out_across = tmp_out %>%
  filter(P.adj < 0.05)
write.table(tmp_out_across, paste0(filedir, "magma_celltype_step1_significant_across.txt"), quote=F, row.names=F, sep="\t")

rm(tmp_out)