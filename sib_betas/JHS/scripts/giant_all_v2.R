#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
############################
#Load packages
library("optparse")
library(data.table)
library(dplyr)
library(biomaRt)
library(vcfR)
#************************************************************************
# read in gwas data and list of 697 SNPs highlighted in the wood paper **
#************************************************************************
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat('Chromosome is ', args[1], '\n')

dir<-paste0(getwd(), "/")
dir1<-"/project/mathilab/bbita/gwas_admix/height_prediction/sib_betas/"
fread(paste0(dir1, 'input/sib_beta_gwas_p.txt'))-> ukb_height #read in GWAS summary statistics for height from the Uk Biobank
ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]
ukb_height[order(CHR,POS)]-> ukb_height
ukb_height[CHR==args[1]]-> ukb_height
setkey(ukb_height, MarkerName)
ukb_height[order(CHR,POS)]-> ukb_height

#*************************
# select UKB  SNPs *
#*************************
fread(paste0('zcat ', dir, 'output/hei_SNPs_chr', args[1],'.vcf.gz'))-> dt
colnames(dt)[1:4]<-c('CHR','POS', 'MarkerName', 'REF')
setkey(dt, CHR, POS)
unique(dt)-> dt
dt[REF %in% c("A","C","G","T")]-> dt
setkey(ukb_height, CHR, POS)
setkey(dt, CHR, POS)
dt[ukb_height, nomatch=0]-> hei # 

#*************************
# Store in Rds format ****
#*************************
saveRDS(hei, file=paste0(dir, 'output/hei_chr', args[1], '.Rds'))

#*****
#END *
#*****

