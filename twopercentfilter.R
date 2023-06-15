#after filtering 2% 
#setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/simulations/")
#args = commandArgs(trailingOnly=TRUE)

vcf <- read.csv("checkcombos.csv", header=FALSE)

vcf$filterV2=ifelse((325 - 100) < vcf$V7 & vcf$V7 < (325 + 100), "keep", "filter")
vcf$filterV3=ifelse((238 - 100) < vcf$V14 & vcf$V14 < (238 + 100), "keep", "filter")
vcf$filterV4=ifelse((315 - 100) < vcf$V21 & vcf$V21 < (315 + 100), "keep", "filter")
vcf$filterV1=ifelse((249 - 100) < vcf$V28 & vcf$V28 < (249 + 100), "keep", "filter")

#all have to match
vcf=vcf[vcf$filterV2 == "keep" & vcf$filterV3 == "keep" & vcf$filterV4 == "keep" & vcf$filterV1 == "keep", ]
vcf=within(vcf, V1<-data.frame(do.call('rbind', strsplit(as.character(V1), 'congenital_plasma.', fixed=TRUE))))
vcf$final=sub("_[^_]+$", "", vcf$V1$X2)
vcf$final2=sub('.[^.]*$', '', vcf$final)



write.table(vcf$final2,"goodcombos2.csv",col.names=FALSE,row.names = FALSE,quote=FALSE)

