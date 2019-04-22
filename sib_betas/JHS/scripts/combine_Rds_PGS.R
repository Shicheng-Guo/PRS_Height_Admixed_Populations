#!/usr/bin/env Rscript
############################
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
library("optparse")
library(data.table)
library(dplyr)
options(scipen=10)
options(digits=10)


currentdir<-paste0(getwd(), "/")
parentdir<-gsub("scripts/", "")

#print(args[1])
vector('list', length(args))-> res_list
names(res_list)<-args
#print(names(res_list))
for(I in 1:length(args)){
	#print(args[I])
	readRDS(paste0(parentdir, '/output/PGS_JHS_', args[I], '.Rds'))->res_list[[I]]
	print(paste0(args[I], ' done'))
}

saveRDS(res_list,paste0(parentdir, 'output/PGS_JHS.Rds'))
#End

