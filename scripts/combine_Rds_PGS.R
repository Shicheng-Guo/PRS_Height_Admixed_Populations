#!/usr/bin/env Rscript
############################
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
suppressMessages(library("optparse"))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
options(scipen=10)
options(digits=10)


print(args[1])
print(args[2])
home="~/height_prediction/"
dir<-paste0(home, args[1], "/", args[2], "/")
args2<- args[-1][-1]
vector('list', length(args2))-> res_list
names(res_list)<-args2
#print(names(res_list))
for(I in 1:length(args2)){
	#print(args[I])
	readRDS(paste0(dir, "output/PGS_", args[2], "_", args2[I], '.Rds'))->res_list[[I]]
	print(paste0(args2[I], ' done'))
}

saveRDS(res_list, paste0(dir, "output/PGS_", args[2],".Rds"))
#End

