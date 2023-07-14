#100x downsampling and 2% frequency removal
#input and output a vcf so it will go back into terbots sumstatscript cleanly

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/patient_data/")
library('VariantAnnotation')
vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_HAN1urine6d.vcf")
vcf <- readVcf("475_sam_dupl_rm_no5_HAN1plasma28d.vcf")

vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_B103plasma1week.vcf")
vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_B103plasma6mo.vcf")
vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_B103urine1week.vcf")
vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_B103urine6mo.vcf")

########################
#compare number of sites and dsit in filtered v. unfiltered
####does not plot very well ignore
#as far as I can tell the full 23500 sites vs the 3300 viable sites are evenly dist
vcf <- readVcf("475_sam_dupl_rm_no5_plasma6mo.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df3 <- merge(df1,df2, by="row.names")

df2=df[df$QUAL > 0,]
#there is no debating the qual filter that stays
df3=df[df$QUAL > 0,]
#actually lets take a pause on qual filter - if its "any prob at all the sites poly" we dont really want it poly rgith?
#not yt - want to know viable sites we can sample from period 

#SPLIT AO
for(i in 1:length(df3$AO)){
  df3$AO[i] <- paste(df3$AO[i],collapse=" ")

}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$AO=gsub( "\\c", "", as.character(df3$AO))
df3$AO=gsub( "\\(", "", as.character(df3$AO))
df3$AO=gsub( "\\)", "", as.character(df3$AO))
df3$AO=gsub( " ", "", as.character(df3$AO))
df3$AO=gsub( ":", ",", as.character(df3$AO))

library(tidyverse)
#need missing entries to be NA stringploit was just putting first value in third col??? 
df3 = df3 %>% 
  separate(col = AO, sep=",",into= str_c('AO',1:3 )) # color1 through color7 column names


#SPLIT SAR
for(i in 1:length(df3$AO)){
  df3$AO[i] <- paste(df3$AO[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$AO=gsub( "\\c", "", as.character(df3$AO))
df3$AO=gsub( "\\(", "", as.character(df3$AO))
df3$AO=gsub( "\\)", "", as.character(df3$AO))
df3$AO=gsub( " ", "", as.character(df3$AO))
df3$AO=gsub( ":", ",", as.character(df3$AO))

#SPLIT SAF
df3 = df3 %>% 
  separate(col = AO, sep=",",into= str_c('AO',1:3 )) # color1 through color7 column names
#SPLIT SAR
for(i in 1:length(df3$AO)){
  df3$AO[i] <- paste(df3$AO[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$AO=gsub( "\\c", "", as.character(df3$AO))
df3$AO=gsub( "\\(", "", as.character(df3$AO))
df3$AO=gsub( "\\)", "", as.character(df3$AO))
df3$AO=gsub( " ", "", as.character(df3$AO))
df3$AO=gsub( ":", ",", as.character(df3$AO))

library(tidyverse)
#need missing entries to be NA stringploit was just putting first value in third col??? 
df3 = df3 %>% 
  separate(col = AO, sep=",",into= str_c('AO',1:3 )) # color1 through color7 column names






#ok now figure out how to check if any of the trialleic sites are >2%
#act first question is how many sites are trialliec period
dfr=df3[df3$AO1 != "NA" & is.na(df3$AO2) == FALSE & is.na(df3$AO3) == FALSE, ]

dfr$AO1freq= as.numeric(dfr$AO1) / as.numeric(dfr$DP)
dfr$AO2freq= as.numeric(dfr$AO2) / as.numeric(dfr$DP)
dfr$AO3freq= as.numeric(dfr$AO3) / as.numeric(dfr$DP)

dfx=dfr[dfr$AO1freq > 0.02 | dfr$AO2freq > 0.02 | dfr$AO3freq > 0.02, ]

dfx=dfr[dfr$AO1freq > 0.02 & dfr$AO2freq > 0.02 & dfr$AO3freq > 0.02,]
        
        | dfr$AO3freq > 0.02, ]


#so using df3 as the true number of visible sites - how many of them are polymorphic? 
#this is where we can then use qual etc 
df=df3[df3$QUAL > 0,]
df=df3[df3$DP > 100,]

#cut down to 100 frequency
df$RO_down <- round(100*(df$RO/df$DP))
df$AO1_down <- round(100*(as.numeric(df$AO1)/df$DP))
df$AO2_down <- round(100*(as.numeric(df$AO2)/df$DP))
df$AO3_down <- round(100*(as.numeric(df$AO3)/df$DP))


df=df[df$RO_down <= 98,]
df=df[df$AO1_down <=98 & df$AO1_down > 1 | df$AO2_down <=98 & df$AO2_down > 1 | df$AO3_down <=98 & df$AO3_down > 1,] 
df=df[complete.cases(df$Row.names),]













#check and see how much each decreases the number of sites
#first check how many sites are triallelic
#mapping quality ≥30, base quality ≥20, remove sites with strand bias. 
#check yet ANOTHER intermediate
#filtered on SRR SSSSSS and qual but NOT on multoalele yet 
vcf <- readVcf("qual0_475_sam_dupl_rm_no5.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df <- merge(df1,df2, by="row.names")


#df$total = as.numeric(df$AO) + as.numeric(df$RO)
#drop sites less than 100X - this is not downsampling this is just filtering for coverage
#df=df[df$total > 100,]

df=within(df, Row.names<-data.frame(do.call('rbind', strsplit(as.character(Row.names), ':', fixed=TRUE))))
df5=within(df, Row.names$X2<-data.frame(do.call('rbind', strsplit(as.character(Row.names$X2), '_', fixed=TRUE))))

vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_B103plasma6mo.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df <- merge(df1,df2, by="row.names")
df=within(df, Row.names<-data.frame(do.call('rbind', strsplit(as.character(Row.names), ':', fixed=TRUE))))
df6=within(df, Row.names$X2<-data.frame(do.call('rbind', strsplit(as.character(Row.names$X2), '_', fixed=TRUE))))

d1new=as.data.frame(df5$Row.names$X2$X1)
colnames(d1new)=c("V1")
d1new$V2="plasma6mo_unfiltered"

d2new=as.data.frame(df6$Row.names$X2$X1)
colnames(d2new)=c("V1")
d2new$V2="plasma6mo_filtered"

library(lattice)
V6=rbind(d1new,d2new)
V6$V1=as.numeric(V6$V1)
dotplot(as.numeric(V1)~V2, data=V6)

ggplot(data=V6, aes(x=V2, y=V1)) +
  geom_point(size=0.000001)

################this section is nothign ignore this

####################
vcf <- readVcf("qual0_475_sam_dupl_rm_no5_nomulti_B103plasma6mo.vcf")

df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df <- merge(df1,df2, by="row.names")
df$total = as.numeric(df$AO) + as.numeric(df$RO)
#drop sites less than 100X - this is not downsampling this is just filtering for coverage
df=df[df$total > 100,]

#cut down to 100 frequency
df$RO_down <- round(100*(df$RO/df$total))
df$AO_down <- round(100*(as.numeric(df$AO)/df$total))

df=df[df$RO_down <= 98,]
df=df[df$AO_down <= 98,]

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




#plot good_variants.vcf BP vs DP, MQ vs DP, QUAL vs. BP manually  
#vcfR plot packages does not work on freebayes TAGs 
library('VariantAnnotation')
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/patient_data/")
vcf <- readVcf("good_variants_B103plasma1week.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data <- merge(df1,df2, by="row.names")

vcf <- readVcf("good_variants_B103plasma6mo.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data2 <- merge(df1,df2, by="row.names")

vcf <- readVcf("good_variants_B103urine1week.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data3 <- merge(df1,df2, by="row.names")

vcf <- readVcf("good_variants_B103urine6mo.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data4 <- merge(df1,df2, by="row.names")

vcf <- readVcf("good_variants_HAN1plasma6d.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data5 <- merge(df1,df2, by="row.names")

vcf <- readVcf("good_variants_HAN1plasma28d.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data6 <- merge(df1,df2, by="row.names")

vcf <- readVcf("good_variants_HAN1urine6d.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data7 <- merge(df1,df2, by="row.names")

vcf <- readVcf("good_variants_HAN1urine25d.output.vcf")
df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
data8 <- merge(df1,df2, by="row.names")


#histogram of cocverages - graphs populate left to right top to bottom
par(mfrow=c(3,8))
#depth
plot(data$start, data$DP, pch = 19, xlab="position", ylab="read depth",cex.axis = 2)
title("B103 plasma 1 week")
plot(data2$start, data2$DP, pch = 19, xlab="", ylab="",cex.axis = 2)
title("B103 plasma 6 months")
plot(data3$start, data3$DP, pch = 19, xlab="", ylab="",cex.axis = 2)
title("B103 urine 1 week")
plot(data4$start, data4$DP, pch = 19, xlab="", ylab="",cex.axis = 2)
title("B103 urine 6 months")
plot(data5$start, data5$DP, pch = 19, xlab="", ylab="",cex.axis = 2)
title("HANChild1 plasma 6 days")
plot(data6$start, data6$DP, pch = 19, xlab="", ylab="",cex.axis = 2)
title("HANChild1 plasma 28 days")
plot(data7$start, data7$DP, pch = 19, xlab="", ylab="",cex.axis = 2)
title("HANChild1 urine 6 days")
plot(data8$start, data8$DP, pch = 19, xlab="", ylab="",cex.axis = 2)
title("HANChild1 urine 28 days")

#mapping quality 
plot(data$start, data$MQM, pch = 19, xlab="position", ylab="mapping quality", col="blue",cex.axis = 2)
points(data$start, data$MQMR, pch = 19, col="red") 
plot(data2$start, data2$MQM, pch = 19, xlab="", ylab="", col="blue",cex.axis = 2)
points(data2$start, data2$MQMR, pch = 19, col="red") 
plot(data3$start, data3$MQM, pch = 19, xlab="", ylab="", col="blue",cex.axis = 2)
points(data3$start, data3$MQMR, pch = 19, col="red") 
plot(data4$start, data4$MQM, pch = 19, xlab="", ylab="", col="blue",cex.axis = 2)
points(data4$start, data4$MQMR, pch = 19, col="red") 
plot(data5$start, data5$MQM, pch = 19, xlab="", ylab="", col="blue",cex.axis = 2)
points(data5$start, data5$MQMR, pch = 19, col="red") 
plot(data6$start, data6$MQM, pch = 19, xlab="", ylab="", col="blue",cex.axis = 2)
points(data6$start, data6$MQMR, pch = 19, col="red") 
plot(data7$start, data7$MQM, pch = 19, xlab="", ylab="", col="blue",cex.axis = 2)
points(data7$start, data7$MQMR, pch = 19, col="red") 
plot(data8$start, data8$MQM, pch = 19, xlab="", ylab="", col="blue",cex.axis = 2)
points(data8$start, data8$MQMR, pch = 19, col="red") 

#base quality
plot(data$start, data$QA, pch = 19, xlab="position", ylab="base quality (Phred)", col="blue",cex.axis = 2)
points(data$start, data$QR, pch = 19, col="red") 
plot(data2$start, data2$QA, pch = 19, xlab="n", ylab="", col="blue",cex.axis = 2)
points(data2$start, data2$QR, pch = 19, col="red") 
plot(data3$start, data3$QA, pch = 19, xlab="n", ylab="", col="blue",cex.axis = 2)
points(data3$start, data3$QR, pch = 19, col="red") 
plot(data4$start, data4$QA, pch = 19, xlab="n", ylab="", col="blue",cex.axis = 2)
points(data4$start, data4$QR, pch = 19, col="red") 
plot(data5$start, data5$QA, pch = 19, xlab="n", ylab="", col="blue",cex.axis = 2)
points(data5$start, data5$QR, pch = 19, col="red")
plot(data6$start, data6$QA, pch = 19, xlab="n", ylab="", col="blue",cex.axis = 2)
points(data6$start, data6$QR, pch = 19, col="red")
plot(data7$start, data7$QA, pch = 19, xlab="n", ylab="", col="blue",cex.axis = 2)
points(data7$start, data7$QR, pch = 19, col="red")
plot(data8$start, data8$QA, pch = 19, xlab="n", ylab="", col="blue",cex.axis = 2)
points(data8$start, data8$QR, pch = 19, col="red")


#ofiicla official final final figure
library(dplyr)
library(tidyr)
library(ggplot2)
#have to filter by congenital and plasma and 100 vs 1000X sampling 
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/")
data=read.csv("allfinalsnps_andcombos.csv",header=FALSE)
data2 = data %>%
  extract(V1, into = c("First", "Second"), "^([^.]+)\\.(.*)")

data2$Second <- sub(".[^.]+$", "", data2$Second)

data3 = data2 %>%
  separate(Second, into = c("combo", "suffix"), sep = "\\.(?=[^.]+$)")

data3$compartment_sampling <- paste(data3$First,data3$suffix)


plot1=data3[data3$First == "congenital_plasma" & data3$suffix == "100", ]
plot2=data3[data3$First == "congenital_plasma" & data3$suffix == "1000", ]
plot3=data3[data3$First == "congenital_urine" & data3$suffix == "100", ]
plot4=data3[data3$First == "congenital_urine" & data3$suffix == "1000", ]



ggplot() +
  geom_bar(data=plot1, aes(combo), position="dodge") +
  geom_jitter(data=data3, aes(x=combo, y=V7, col=compartment_sampling)) + 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) 


#jeff wants figs seperates and dotplots to be a box/whisker plot
ggplot() +
  geom_bar(data=plot1, aes(combo), position="dodge") + 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) 



# grouped boxplot
ggplot(data3, aes(x=combo, y=V7, fill=compartment_sampling)) + 
  geom_boxplot() + 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  + geom_hline(yintercept=20)

#single boxplots

#P100 is 325 - p1000 is 238 - u100 is 315 and u1000 is 249
#change this for each yintercept 
ggplot(plot4, aes(x=combo, y=V7)) + 
  geom_boxplot() + 
  theme(axis.text = element_text(face="bold"), # make the axis labels bold to make it more readable
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1), # make the x-axis labels vertical for readability
        axis.title.x = element_blank(), # next two lines are to remove axis titles
        axis.title.y = element_blank()) +
  geom_hline(yintercept=249) +
  ggtitle("Urine Compartment 1000 Samples")


