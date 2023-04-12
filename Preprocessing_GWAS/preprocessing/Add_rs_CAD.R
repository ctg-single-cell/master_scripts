library(utils)
setwd("/home/brouwer2/sources/sumstats/input_GWAS")

data <- read.table("GCST90132314_buildGRCh37.tsv",header =T)
rsIDS <- read.table("/home/brouwer2/sources/magma_1.10/g1000_eur/g1000_eur.bim",header=F)
names(rsIDS) <- c("chromosome","SNP","None","base_pair_location","a1","a2")
mdata <- merge(data,rsIDS,by=c("chromosome","base_pair_location"))
print(dim(mdata))
keep <-  which((toupper(as.character(mdata$effect_allele)) == as.character(mdata$a1) & toupper(as.character(mdata$other_allele)) == as.character(mdata$a2)) | (toupper(as.character(mdata$effect_allele)) == as.character(mdata$a2) & toupper(as.character(mdata$other_allele)) == as.character(mdata$a1)))
print(length(keep))
mdata <- mdata[keep,c("SNP","effect_allele","other_allele","odds_ratio","p_value","n")]
print(dim(mdata))
write.table(mdata,"CAD_reformatted",row.names=F,quote=F,col.names=T)
