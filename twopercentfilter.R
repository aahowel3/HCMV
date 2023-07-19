#after filtering 2% 
#setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/simulations/")
#args = commandArgs(trailingOnly=TRUE)

#vcf <- read.csv("checkcombos.csv", header=FALSE, sep=","
vcf1 <- read.csv("congenital.plasma.100.check", header=FALSE, sep=",")
vcf2 <- read.csv("congenital.plasma.1000.check", header=FALSE, sep=",")
vcf3 <- read.csv("congenital.urine.100.check", header=FALSE, sep=",")
vcf4 <- read.csv("congenital.urine.1000.check", header=FALSE, sep=",")

vcf=cbind(vcf1,vcf2,vcf3,vcf4)
colnames(vcf)=paste("V",seq(1,28,by=1),sep="")

vcf$filterV2=ifelse((440 - 100) < vcf$V7 & vcf$V7 < (440 + 100), "keep", "filter")
vcf$filterV3=ifelse((305 - 100) < vcf$V14 & vcf$V14 < (305 + 100), "keep", "filter")
vcf$filterV4=ifelse((439 - 100) < vcf$V21 & vcf$V21 < (439 + 100), "keep", "filter")
vcf$filterV1=ifelse((341 - 100) < vcf$V28 & vcf$V28 < (341 + 100), "keep", "filter")

#all have to match
vcf=vcf[vcf$filterV2 == "keep" & vcf$filterV3 == "keep" & vcf$filterV4 == "keep" & vcf$filterV1 == "keep", ]
vcf=within(vcf, V1<-data.frame(do.call('rbind', strsplit(as.character(V1), 'congenital_plasma.', fixed=TRUE))))
vcf$final=sub(".sub", '', vcf$V1$X2)

vcf$final2=sub("_[^_]+$", "", vcf$final)
vcf$final3=sub('.[^.]*$', '', vcf$final2)




write.table(vcf$final3,"goodcombos2.csv",col.names=FALSE,row.names = FALSE,quote=FALSE)

