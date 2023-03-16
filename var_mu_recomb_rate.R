#generate variable mutation and recombination rate map for slim simulations
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/HCMV_slim/")

pos=seq(1, 23500, 500)
recomb_rate=runif(length(pos), min = 9.8e-8,max = 9.8e-6)
df1<-data.frame(pos,recomb_rate)

mu_rate=runif(length(pos), min = 2e-8, max = 2e-6)
df2<-data.frame(pos,mu_rate)

write.table(df1, "recomb.var.txt",row.names = FALSE, quote=FALSE,sep="\t",col.names=FALSE)
write.table(df2, "mu.var.txt",row.names = FALSE, quote=FALSE,sep="\t",col.names=FALSE)


#create params file for slims
DFE = c("DFE1,0.5,0.24,0.12,0.13","DFE2,0.1,0.7,0.1,0.1",
        "DFE3,0.1,0.1,0.1,0.7")
recomb = c("fix,9.8e-6","fix,9.8e-7","fix,9.8e-8","var,0")
mu = c("fix,2e-6","fix,2e-7","fix,2e-8","var,0")
progeny = c(0.01,0.07,0.1)
burnin= c("with","without")
gr= c(0.38,1)


df=expand.grid(DFE, recomb, mu, progeny, burnin, gr)
df_burnin=expand.grid(DFE, recomb, mu, progeny, gr)

df=within(df, Var1<-data.frame(do.call('rbind', strsplit(as.character(Var1), ',', fixed=TRUE))))
df_burnin=within(df_burnin, Var1<-data.frame(do.call('rbind', strsplit(as.character(Var1), ',', fixed=TRUE))))


df=within(df, Var2<-data.frame(do.call('rbind', strsplit(as.character(Var2), ',', fixed=TRUE))))
df_burnin=within(df_burnin, Var2<-data.frame(do.call('rbind', strsplit(as.character(Var2), ',', fixed=TRUE))))


df=within(df, Var3<-data.frame(do.call('rbind', strsplit(as.character(Var3), ',', fixed=TRUE))))
df_burnin=within(df_burnin, Var3<-data.frame(do.call('rbind', strsplit(as.character(Var3), ',', fixed=TRUE))))

write.csv(df, "params_congenital.txt",row.names = FALSE,quote=c(1,6,8,11))
write.csv(df_burnin, "params_burnin.txt",row.names = FALSE, quote=c(1,6,8))


