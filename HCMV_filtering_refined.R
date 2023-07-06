########################
library(tidyverse)
library('VariantAnnotation')

setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/patient_data/")
#compare number of sites and dsit in filtered v. unfiltered
####does not plot very well ignore
#as far as I can tell the full 23500 sites vs the 3300 viable sites are evenly dist
vcf <- readVcf("475_sam_dupl_rm_no5_plasma6mo.vcf")
vcf <- readVcf("475_sam_dupl_rm_no5_urine6mo.vcf")

df1 = as.data.frame(info(vcf))
df2=as.data.frame(rowRanges(vcf))
df3 <- merge(df1,df2, by="row.names")


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

#need missing entries to be NA stringploit was just putting first value in third col??? 
df3 = df3 %>% 
  separate(col = AO, sep=",",into= str_c('AO',1:3 )) # color1 through color7 column names


#SPLIT SAR
for(i in 1:length(df3$SAR)){
  df3$SAR[i] <- paste(df3$SAR[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$SAR=gsub( "\\c", "", as.character(df3$SAR))
df3$SAR=gsub( "\\(", "", as.character(df3$SAR))
df3$SAR=gsub( "\\)", "", as.character(df3$SAR))
df3$SAR=gsub( " ", "", as.character(df3$SAR))
df3$SAR=gsub( ":", ",", as.character(df3$SAR))

df3 = df3 %>% 
  separate(col = SAR, sep=",",into= str_c('SAR',1:3 )) # color1 through color7 column names

#SPLIT SAF
for(i in 1:length(df3$SAF)){
  df3$SAF[i] <- paste(df3$SAF[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$SAF=gsub( "\\c", "", as.character(df3$SAF))
df3$SAF=gsub( "\\(", "", as.character(df3$SAF))
df3$SAF=gsub( "\\)", "", as.character(df3$SAF))
df3$SAF=gsub( " ", "", as.character(df3$SAF))
df3$SAF=gsub( ":", ",", as.character(df3$SAF))

df3 = df3 %>% 
  separate(col = SAF, sep=",",into= str_c('SAF',1:3 )) # color1 through color7 column names


#SPLIT RPR reads placed right
for(i in 1:length(df3$RPR)){
  df3$RPR[i] <- paste(df3$RPR[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$RPR=gsub( "\\c", "", as.character(df3$RPR))
df3$RPR=gsub( "\\(", "", as.character(df3$RPR))
df3$RPR=gsub( "\\)", "", as.character(df3$RPR))
df3$RPR=gsub( " ", "", as.character(df3$RPR))
df3$RPR=gsub( ":", ",", as.character(df3$RPR))

df3 = df3 %>% 
  separate(col = RPR, sep=",",into= str_c('RPR',1:3 )) # color1 through color7 column names

#SPLIT RPRL reads placed left 
for(i in 1:length(df3$RPL)){
  df3$RPL[i] <- paste(df3$RPL[i],collapse=" ")
  
}
#what an absolutele pain in the ass to get it into a stringsplit format
df3$RPL=gsub( "\\c", "", as.character(df3$RPL))
df3$RPL=gsub( "\\(", "", as.character(df3$RPL))
df3$RPL=gsub( "\\)", "", as.character(df3$RPL))
df3$RPL=gsub( " ", "", as.character(df3$RPL))
df3$RPL=gsub( ":", ",", as.character(df3$RPL))

df3 = df3 %>% 
  separate(col = RPL, sep=",",into= str_c('RPL',1:3 )) # color1 through color7 column names



#so using df3 as the true number of visible sites - how many of them are polymorphic? 
#this is where we can then use qual etc 
df=df3[df3$DP > 100,]
df=df[df$QUAL > 0,]



#cut down to 100 frequency
df$RO_down <- round(100*(df$RO/df$DP))
df$AO1_down <- round(100*(as.numeric(df$AO1)/df$DP))
df$AO2_down <- round(100*(as.numeric(df$AO2)/df$DP))
df$AO3_down <- round(100*(as.numeric(df$AO3)/df$DP))

#get rid of sites where ref is greater than 98 - if 99 other sites cant be greater than 2%
df=df[df$RO_down <= 98,]
#get rid of sites where no alt is at least greater than 2% but also less than 99
df=df[df$AO1_down <=98 & df$AO1_down > 1 | df$AO2_down <=98 & df$AO2_down > 1 | df$AO3_down <=98 & df$AO3_down > 1,] 
df=df[complete.cases(df$Row.names),]

###
#if were saying the site is truly a MAF (greater than 2%) then it HAS to be supported by more than zero reads 
#df <- df[(df$RO_down >= 2 & df$SRF > 0 & df$SRR > 0) | (df$AO1_down >= 2 & df$SAR1 > 0 & df$SAF1 > 0) | (df$AO2_down >= 2 & df$SAR2 > 0 & df$SAF2 > 0) | 
#           (df$AO3_down >= 2 & df$SAR3 > 0 & df$SAF3 > 0), ]
#df=df[complete.cases(df$Row.names),]

##
##THIS - THIS IS THE RIGHT IDEA 
#NEEDS TO BE AT LEAST TWO SITES MEETTHING THESE CONDITIONS - NOT OROROR 
#OTHERWISE ITS NOT TRULY POLYMORPHIC 
#SOMETHINGS BEING WEIRD WITH THE NA THO
df=as.data.frame(df, stringsAsFactors = FALSE)
#perfect! just throws a warning not an error!!
df[is.na(df)] <- 0

###whoops whoops cant od this 
#all those sites still need to be supported by at least one read one ither sites and at that site 
##tris=subset(df5,  rowSums(df5[65:68] >= 2) == 3)
#PREV CHUNK WAS FOR TOTAL THIS NOW BREAKS DOWN BY BI/TRI/QUAD ALLELEIC
# == INSTEAD OF >=
df5 <- df[rowSums(cbind((df$RO_down >= 2 & df$SRR > 0 & df$SRF > 0), (df$AO1_down >= 2 & df$SAR1 > 0 & df$SAF1 > 0), 
                        (df$AO2_down >= 2 & df$SAR2 > 0 & df$SAF2 > 0),  
                        (df$AO3_down >= 2 & df$SAR3 > 0 & df$SAF3 > 0))) == 3,]

df5=df5[complete.cases(df5$Row.names),]

#need to do the same thing but with RPR and RPR and RPL
#there doesnt seem to be a reference equivalent for that though just the alternaties 
###yesssssss
#IF its not the ref allele thats being considered - then the alts MUST be supported 

df5=subset(df5, ((RO_down >= 2) & rowSums(cbind((df5$AO1_down >= 2 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                         (df5$AO2_down >= 2 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                         (df5$AO3_down >= 2 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 1) | 
             ((RO_down < 2) & rowSums(cbind((df5$AO1_down >= 2 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                                             (df5$AO2_down >= 2 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                                              (df5$AO3_down >= 2 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 2)) 
df5=df5[complete.cases(df5$Row.names),]



###for 1000x downsampling 
#so using df3 as the true number of visible sites - how many of them are polymorphic? 
#this is where we can then use qual etc 
df=df3[df3$DP > 100,]
df=df[df$QUAL > 0,]



#cut down to 100 frequency
df$RO_down <- round(1000*(df$RO/df$DP))
df$AO1_down <- round(1000*(as.numeric(df$AO1)/df$DP))
df$AO2_down <- round(1000*(as.numeric(df$AO2)/df$DP))
df$AO3_down <- round(1000*(as.numeric(df$AO3)/df$DP))


df=df[df$RO_down <= 980,]
df=df[df$AO1_down <=980 & df$AO1_down >= 20 | df$AO2_down <=980 & df$AO2_down >= 20 | df$AO3_down <=980 & df$AO3_down >= 20,] 
df=df[complete.cases(df$Row.names),]

df=as.data.frame(df, stringsAsFactors = FALSE)
#perfect! just throws a warning not an error!!
df[is.na(df)] <- 0

###yesssssss
df5 <- df[rowSums(cbind((df$RO_down >= 20 & df$SRR > 0 & df$SRF > 0), (df$AO1_down >= 20 & df$SAR1 > 0 & df$SAF1 > 0), 
                        (df$AO2_down >= 20 & df$SAR2 > 0 & df$SAF2 > 0),  
                        (df$AO3_down >= 20 & df$SAR3 > 0 & df$SAF3 > 0))) == 4,]
df5=df5[complete.cases(df5$Row.names),]


df5=subset(df5, ((RO_down >= 20) & rowSums(cbind((df5$AO1_down >= 20 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                                                (df5$AO2_down >= 20 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                                                (df5$AO3_down >= 20 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 3) | 
             ((RO_down < 20) & rowSums(cbind((df5$AO1_down >= 20 & df5$RPR1 > 0 & df5$RPL1 > 0), 
                                            (df5$AO2_down >= 20 & df5$RPR2 > 0 & df5$RPL2 > 0),  
                                            (df5$AO3_down >= 20 & df5$RPR3 > 0 & df5$RPL3 > 0))) >= 4)) 
df5=df5[complete.cases(df5$Row.names),]





