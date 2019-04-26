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


print(args[1])
print(args[2])
print(args[3])
home="~/height_prediction/"
dir<-paste0(home, args[2], "/", args[3], "/")
	test<-vector('list', 22)
	#print(args[I])
	for(X in 1:22){
	fread(paste0('zcat ', dir, 'output/hei_SNPs_chr', X, '.vcf.gz'))-> vcf
	setDT(vcf)
	colnames(vcf)[3]<-'MarkerName'
	setkey(vcf, MarkerName)
	readRDS(paste0(dir, 'prunned_1kg/LD_prunned_hei_chr', X, '_',args[1], '.Rds'))[['keep']]-> bla
	data.table(MarkerName=bla)-> test[[X]]
	colnames(test[[X]])<-'MarkerName'
	setkey(test[[X]], MarkerName)
	vcf[test[[X]], nomatch=0]-> test[[X]]
	}
	#do.call(rbind, test)-> testB
	names(test)<-seq(1:22)
	print(paste0(args[1], ' done'))


saveRDS(test,paste0(dir, 'output/hei_', args[1], '.Rds'))
#End
