# This script aims to get the output from ewce into a text file
library(EWCE)
ctd = "c://Users/tph205/Documents/ctg/temp/20230911_gwas_ct_pipeline/scrna/41_Siletti_Cerebellum.CBV_Human_2022/41_Siletti_Cerebellum.CBV_Human_2022_siletti.rda"
ctd = EWCE::load_rdata(ctd)

# specificity
specificity_df = as.data.frame(as.matrix(ctd[[1]]$specificity))

#replace space with _
names(specificity_df) <- gsub(" ", "_", names(specificity_df))

#convert rowname to a column
specificity_df$symbol <- row.names(specificity_df)

# convert gene symbol to ensemble
conversion_file = fread("c://Users/tph205/Documents/ctg/codes/master_scripts/Conversion_gene_names/conversion_files/gene_names_human.txt")

conversion_file_short = conversion_file %>%
  select("Approved symbol", "Ensembl gene ID")
colnames(conversion_file_short) = c("symbol", "GENE")

specificity_df_ensembl = merge(specificity_df, conversion_file_short, by = c("symbol"))

specificity_df_ensembl_fmt = specificity_df_ensembl %>%
  select(-symbol) %>%
  select(GENE, everything())
write.table(specificity_df_ensembl_fmt, "c://Users/tph205/Documents/ctg/temp/20230911_gwas_ct_pipeline/41_Siletti_Cerebellum.CBV_Human_2022/spec_ewce_linear.tsv", quote = F, sep = "\t", row.names = FALSE)

# mean
mean_df = as.data.frame(as.matrix(ctd[[1]]$mean_exp))

#replace space with _
names(mean_df) <- gsub(" ", "_", names(mean_df))

#convert rowname to a column
mean_df$symbol <- row.names(mean_df)

mean_df_ensembl = merge(mean_df, conversion_file_short, by = c("symbol"))

mean_df_ensembl_fmt = mean_df_ensembl %>%
  select(-symbol) %>%
  select(GENE, everything())

mean_df_ensembl_fmt$Average = rowMeans(mean_df_ensembl_fmt[,2:15])

write.table(mean_df_ensembl_fmt, "c://Users/tph205/Documents/ctg/temp/20230911_gwas_ct_pipeline/scrna/41_Siletti_Cerebellum.CBV_Human_2022/mean_ewce_linear.tsv", quote = F, sep = "\t", row.names = FALSE)


# spec_vanilla = read.table("c://Users/tph205/Documents/ctg/temp/20230911_gwas_ct_pipeline/41_Siletti_Cerebellum.CBV_Human_2022/spec_vanilla_linear.tsv", header = T)
# spec_cellex = read.table("c://Users/tph205/Documents/ctg/temp/20230911_gwas_ct_pipeline/41_Siletti_Cerebellum.CBV_Human_2022/spec_cellex_linear.tsv", header = T)
# 
# tmp = merge(specificity_df_ensembl, spec_vanilla, by = c("GENE"))
# merged = merge(tmp, spec_cellex, by = c("GENE"))