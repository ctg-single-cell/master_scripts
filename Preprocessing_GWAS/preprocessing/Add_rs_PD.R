library(utils)
library(stringr)

setwd("/home/brouwer2/sources/sumstats/input_GWAS")

data <- read.table("nallsEtAl2019_excluding23andMe_allVariants.tab",header =T)
rsIDS <- read.table("/home/brouwer2/sources/magma_1.09a/g1000_eur/g1000_eur.bim",header=F)
data[,c("chr","PosGRCh37")] <- str_split_fixed(str_split_fixed(data$SNP,"chr",2)[,2],":",2)[,c(1,2)]
data$chr <- as.integer(data$chr)
data$PosGRCh37 <- as.integer(data$PosGRCh37)

names(rsIDS) <- c("chr","rsID","None","PosGRCh37","a1","a2")
mdata <- merge(data,rsIDS,by=c("chr","PosGRCh37"))
print(dim(mdata))
keep <-  which((as.character(mdata$A1) == as.character(mdata$a1) & as.character(mdata$A2) == as.character(mdata$a2)) | (as.character(mdata$A1) == as.character(mdata$a2) & as.character(mdata$A2) == as.character(mdata$a1)))
print(length(keep))

mdata$Neff <- 2/(1/mdata$N_cases + 1/mdata$N_controls) 
print(mean(mdata$Neff))

mdata <- mdata[keep,c("rsID","b","se","A1","A2","p","Neff","N_cases","N_controls")]
print(dim(mdata))
write.table(mdata,"PD_reformatted",row.names=F,quote=F,col.names=T)


