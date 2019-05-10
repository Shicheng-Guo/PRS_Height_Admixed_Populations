#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) #arg1 is vcf name, arg2 is panel.txt name

ind.filename <- paste(args[1], ".ind", sep = "")
ind <- read.table(ind.filename, as.is=TRUE)
panel <- read.table(args[2], as.is=TRUE)

pmap<-panel[,2]
names(pmap)<-panel[,1]

ind[,3]<-pmap[ind[,1]]

ind.renamed.filename <- paste(args[1], "_renamed_pops.ind", sep = "")
write.table(ind, ind.renamed.filename, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
