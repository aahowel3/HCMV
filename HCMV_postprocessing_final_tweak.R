###########remove <2% variants to feed into terbots script to caclc only with threhold vars
##setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/simulations/for_terbotscript/")
#vcf <- read.csv("congenital_plasma.DFE3.2.0e-06.0.0.01.0.38.without.1000.output.ms", header=FALSE)

library(dplyr)
args = commandArgs(trailingOnly=TRUE)

vcf <- read.csv(args[1], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=vcftotal %>%
  select(where(function(x) any(x %in% subsam)))

library(janitor)
vcftotal2=vcftotal2 %>%
  row_to_names(row_number = 1)

bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]

bionly=vcftotal2[ , -which(names(vcftotal2) %in% bad)]
trionly=vcftotal2[ , which(names(vcftotal2) %in% bad)]
#trionly2=trionly[duplicated(colnames(trionly)), fromLast = TRUE, ]
#trionly2=trionly[, grep("\\.1$", colnames(trionly))]

bionly=as.data.frame(lapply(bionly, as.numeric))
bionly2percent = bionly[,colSums(bionly[-1,]) >= 2]
bionly2percent = bionly2percent[,colSums(bionly2percent[-1,]) <= 98]

con_plasma_100_bi=length(bionly2percent)

vcf <- read.csv(args[2], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=vcftotal %>%
  select(where(function(x) any(x %in% subsam)))

library(janitor)
vcftotal2=vcftotal2 %>%
  row_to_names(row_number = 1)

bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]

bionly=vcftotal2[ , -which(names(vcftotal2) %in% bad)]
trionly=vcftotal2[ , which(names(vcftotal2) %in% bad)]
#trionly2=trionly[duplicated(colnames(trionly)), fromLast = TRUE, ]
#trionly2=trionly[, grep("\\.1$", colnames(trionly))]

bionly=as.data.frame(lapply(bionly, as.numeric))
bionly2percent = bionly[,colSums(bionly[-1,]) >= 20]
bionly2percent = bionly2percent[,colSums(bionly2percent[-1,]) <= 980]

con_plasma_1000_bi=length(bionly2percent)


vcf <- read.csv(args[3], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=vcftotal %>%
  select(where(function(x) any(x %in% subsam)))

library(janitor)
vcftotal2=vcftotal2 %>%
  row_to_names(row_number = 1)

bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]

bionly=vcftotal2[ , -which(names(vcftotal2) %in% bad)]
trionly=vcftotal2[ , which(names(vcftotal2) %in% bad)]
#trionly2=trionly[duplicated(colnames(trionly)), fromLast = TRUE, ]
#trionly2=trionly[, grep("\\.1$", colnames(trionly))]

bionly=as.data.frame(lapply(bionly, as.numeric))
bionly2percent = bionly[,colSums(bionly[-1,]) >= 2]
bionly2percent = bionly2percent[,colSums(bionly2percent[-1,]) <= 98]

con_urine_100_bi=length(bionly2percent)

vcf <- read.csv(args[4], header=FALSE)
vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))
vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)
###wihtout looking at unique only
#trialleics are twice as likely to stay or get removed since theyre represented as TWO columns with the same number
tranfor=as.data.frame(unique(t(vcftotal[1,])))
colnames(tranfor)="V1"
invar=runif(23500-length(unique(t(vcftotal[1,]))))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)
#but then when you subsample by the numbers remainnig will salvage both columns 
vcftotal2=vcftotal %>%
  select(where(function(x) any(x %in% subsam)))


library(janitor)
vcftotal2=vcftotal2 %>%
  row_to_names(row_number = 1)

bad=colnames(vcftotal2)[duplicated(colnames(vcftotal2))]

bionly=vcftotal2[ , -which(names(vcftotal2) %in% bad)]
trionly=vcftotal2[ , which(names(vcftotal2) %in% bad)]
#trionly2=trionly[duplicated(colnames(trionly)), fromLast = TRUE, ]
#trionly2=trionly[, grep("\\.1$", colnames(trionly))]

bionly=as.data.frame(lapply(bionly, as.numeric))
bionly2percent = bionly[,colSums(bionly[-1,]) >= 20]
bionly2percent = bionly2percent[,colSums(bionly2percent[-1,]) <= 980]

con_urine_1000_bi=length(bionly2percent)



ifelse(((431 - 100) <= con_plasma_100_bi &  con_plasma_100_bi <= (431 + 100)) & 
         ((304 - 100) <= con_plasma_1000_bi & con_plasma_1000_bi <= (304 + 100)) & 
         ((431 - 100) <= con_urine_100_bi & con_urine_100_bi <= (431 + 100)) & 
         ((335 - 100) <= con_urine_1000_bi & con_urine_1000_bi <= (335 + 100)), paste0("good", " ", args[1], " ", con_plasma_100_bi, " ",con_plasma_1000_bi, " ",con_urine_100_bi, " ", con_urine_1000_bi), 
	paste0("bad", " ", args[1], " ", con_plasma_100_bi, " ",con_plasma_1000_bi, " ",con_urine_100_bi, " ", con_urine_1000_bi))







