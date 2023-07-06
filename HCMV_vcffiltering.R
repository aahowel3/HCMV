
####################
vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_B103plasma6mo.vcf")

df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df <- merge(df1,df2, by="row.names")
df$total = as.numeric(df$AO) + as.numeric(df$RO)
#drop sites less than 100X - this is not downsampling this is just filtering for coverage
df=df[df$total > 100,]

#cut down to 100 frequency
df$RO_down <- round(1000*(df$RO/df$total))
df$AO_down <- round(1000*(as.numeric(df$AO)/df$total))

df=df[df$RO_down <= 980,]
df=df[df$AO_down <= 980,]

#special case for hanlpasma28d - var with ref/. geno called but no location bc not real var 
#drop to match vcf length of vcfR read in
#df=df[!is.na(df$Row.names),]

#in new vcf - replace RO values here with RO from col above 
df$RO_down <- paste(";RO=", df$RO_down, sep="")
df$RO_down <- paste(df$RO_down, ";", sep="")

#in new vcf - replace AO values here with AO from col above 
df$AO_down <- paste(";AO=", df$AO_down, sep="")
df$AO_down <- paste(df$AO_down, ";", sep="")

#easier packacge to write to
library(vcfR)
library(stringr)
vcf <- read.vcfR("475_sam_dupl_rm_no5_HAN1plasma28d.vcf")

vcf <- read.vcfR("qual0_475_sam_dupl_rm_no5_nomulti_B103plasma6mo.vcf")
vcf <- read.vcfR("qual0_475_sam_dupl_rm_no5_nomulti_B103urine6mo.vcf")


#1 - cutdown to length of filter file above 
vcf@fix=vcf@fix[vcf@fix[,'POS'] %in% df$start,]
#this works - but it does not change the rowID - no way to link @gt lines and @fix lines 
#apparently you dont need to drop the corresponding @gt lines just writes it out correct - cool 
vcf@fix[,'INFO']=str_replace_all(vcf@fix[,'INFO'], ";RO=[[:digit:]]+\\;", as.character(df$RO_down))
vcf@fix[,'INFO']=str_replace_all(vcf@fix[,'INFO'], ";AO=[[:digit:]]+\\;", as.character(df$AO_down))

#will output it compressed regardless of filename
write.vcf(vcf, file="good_variants_B103urine6mo.1000.vcf.gz")


#reutrn param combinations within an order of mag of patient sample value

vcf <- read.csv("B103_plasmasim100.csv", header=FALSE)

#after filtering 2% 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/simulations/")
vcf <- read.csv("check_paramcombos.csv", header=FALSE)
vcf=within(vcf, V1<-data.frame(do.call('rbind', strsplit(as.character(V1), ',', fixed=TRUE))))
vcf=within(vcf, V2<-data.frame(do.call('rbind', strsplit(as.character(V2), ',', fixed=TRUE))))
vcf=within(vcf, V3<-data.frame(do.call('rbind', strsplit(as.character(V3), ',', fixed=TRUE))))
vcf=within(vcf, V4<-data.frame(do.call('rbind', strsplit(as.character(V4), ',', fixed=TRUE))))

vcf$V2$X7=as.numeric(vcf$V2$X7)
vcf$V3$X7=as.numeric(vcf$V3$X7)
vcf$V4$X7=as.numeric(vcf$V4$X7)
vcf$V1$X7=as.numeric(vcf$V1$X7)


vcf$filterV2=ifelse((325 - 100) < vcf$V2$X7 & vcf$V2$X7 < (325 + 100), "keep", "filter")
vcf$filterV3=ifelse((238 - 100) < vcf$V3$X7 & vcf$V3$X7 < (238 + 100), "keep", "filter")
vcf$filterV4=ifelse((315 - 100) < vcf$V4$X7 & vcf$V4$X7 < (315 + 100), "keep", "filter")
vcf$filterV1=ifelse((249 - 100) < vcf$V1$X7 & vcf$V1$X7 < (249 + 100), "keep", "filter")

#all have to match
vcf=vcf[vcf$filterV2 == "keep" & vcf$filterV3 == "keep" & vcf$filterV4 == "keep" & vcf$filterV1 == "keep", ]
vcf=within(vcf, V1$X1<-data.frame(do.call('rbind', strsplit(as.character(V1$X1), 'congenital_plasma.', fixed=TRUE))))
write.table(vcf$V1$X1$X2,"goodcombos2.csv",col.names=FALSE,row.names = FALSE,quote=FALSE)



#what if we only care about the 100x sampled ones 
vcf=vcf[vcf$filterV2 == "keep" & vcf$filterV4 == "keep", ]
vcf=vcf[vcf$filterV3 == "keep" & vcf$filterV5 == "keep", ]

#fuck  just look at em all lined up
vcf$l = vcf$V1$X1
vcf$l1 = vcf$V1$X3
vcf$l2 = vcf$V1$X6

vcf$new1=vcf$V1$X7
vcf$new2=vcf$V2$X7
vcf$new3=vcf$V3$X7
vcf$new4=vcf$V4$X7
vcf2=vcf[c("l","l1","l2","new1","new2","new3","new4")]

write.csv(vcf2,"goodcombos2.csv")


###########remove >2% variants to feed into terbots script to caclc only with threhold vars
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/simulations/for_terbotscript/")
vcf <- read.csv("congenital_plasma.DFE1.0.0.0.1.0.38.with.100.output.ms", header=FALSE)

vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))

vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)

df2=as.data.frame(lapply(vcftotal, as.numeric))
dd1 = df2[,colSums(df2[-1,]) > 2]

