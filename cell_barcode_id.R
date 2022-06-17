.libPaths("/ihome/tbruno/arc85/Rlibs_Oct_2019")
library(Seurat)

cli <- commandArgs(trailingOnly=TRUE)
args <- strsplit(cli,"=",fixed=T)

dat <- Read10X(args[[1]])
write.table(colnames(dat),file=paste(args[[2]],"_whitelist.csv",sep=""),quote=F,row.names=F,col.names=F,sep=",")
