library(utils)
setwd("/home/brouwer2/sources/sumstats/input_GWAS")

data <- read.table("PGCALZ2sumstatsExcluding23andMe.txt",header =T)
rsIDS <- read.table("/home/brouwer2/sources/magma_1.09a/g1000_eur/g1000_eur.bim",header=F)
names(rsIDS) <- c("chr","SNP","None","PosGRCh37","a1","a2")
mdata <- merge(data,rsIDS,by=c("chr","PosGRCh37"))
print(dim(mdata))
keep <-  which((as.character(mdata$testedAllele) == as.character(mdata$a1) & as.character(mdata$otherAllele) == as.character(mdata$a2)) | (as.character(mdata$testedAllele) == as.character(mdata$a2) & as.character(mdata$otherAllele) == as.character(mdata$a1)))
print(length(keep))
mdata <- mdata[keep,c("SNP","testedAllele","otherAllele","z","p","N")]
print(dim(mdata))
write.table(mdata,"AD_reformatted",row.names=F,quote=F,col.names=T)
