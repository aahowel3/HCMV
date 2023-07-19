###########subsample to 99.7 to match empirical data
#setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/simulations/for_terbotscript/")

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
tranfor=as.data.frame(t(vcftotal))
invar=runif(23500-length(vcftotal))
whole=c(tranfor$V1,invar)
subsam=sample(whole,23430)

vcftotal2=vcftotal %>%
  select(where(function(x) any(x %in% subsam)))


#check dist of points on the line
#library(lattice)
#library(ggplot2)
#tranfor=as.data.frame(t(vcftotal))
#tranfor$lab="position"
#tranfor$V1=as.numeric(tranfor$V1)
#ggplot(data=tranfor, aes(x=lab, y=V1)) +
#  geom_point(size=0.00000001)




df2=as.data.frame(lapply(vcftotal2, as.numeric))
x1=as.data.frame(df2[1,])
x2=round(df2[-1,], 1)

write( "//", file = args[2], append = T)
cat( "segsites:", length(x2), file = args[2], append = T)
x1[1]=paste("\npositions:", x1[1])
write.table( x1,
             file = args[2],
             append = T,
             sep = " ",
             row.names = F,
             col.names = F,
             na="",
             quote = F)
write.table( x2,
             file = args[2],
             append = T,
             sep = "",
             row.names = F,
             col.names = F,
             na="",
             quote = F)

