###########remove >2% variants to feed into terbots script to caclc only with threhold vars
#setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/simulations/for_terbotscript/")

args = commandArgs(trailingOnly=TRUE)
vcf <- read.csv(args[1], header=FALSE)

vcftail=tail(vcf, -3)
vcftail2=data.frame(do.call("rbind", strsplit(as.character(vcftail$V1), "", fixed = TRUE)))

vcfhead=as.data.frame(vcf[3,])
vcfhead2=data.frame(do.call("rbind", strsplit(as.character(vcfhead$'vcf[3, ]'), " ", fixed = TRUE)))
vcfhead2$X1 <- NULL
names(vcfhead2) <- names(vcftail2) 
vcftotal=rbind(vcfhead2, vcftail2)

df2=as.data.frame(lapply(vcftotal, as.numeric))
dd1 = df2[,colSums(df2[-1,]) > 20]
dd2 = dd1[,colSums(dd1[-1,]) < 980]
x1=as.data.frame(dd2[1,])
x2=round(dd2[-1,], 1)


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


#####cannot get positions and table to print on same line to save my life manually editing it


