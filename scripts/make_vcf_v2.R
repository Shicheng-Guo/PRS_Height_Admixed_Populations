#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
############################
#Load packages
library("optparse")
library(data.table)
library(dplyr)
library(biomaRt)
#************************************************************************
# read in gwas data and list of 697 SNPs highlighted in the wood paper **
#************************************************************************
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat('Chromosome is ', args[1], '\n')

home="/home/bbita/height_prediction/"
dir=paste0(home,  args[2], "/", args[3],"/scripts/")
parentdir<-gsub("/scripts", "", dir)
if(args[2]=='sib_betas'){
        fread(paste0(home, args[2], "/input/sib_beta_gwas_p.txt"))-> ukb_height #read in betas from sibling analyses
} else if (args[2]=='gwas'){
        fread(paste0('zcat ', home, args[2], "/input/50.assoc.tsv.gz"))-> ukb_height #read in GWAS summary statistics for height from the Uk Biobank
        ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
        ukb_height[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
        ukb_height[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
        ukb_height[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height
}
ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]
ukb_height[order(CHR,POS)]-> ukb_height

#*************************
# select UKB  SNPs *
#*************************
fread(paste0('zcat ', parentdir, 'output/hei_SNPs_chr', args[1],'.vcf.gz'))-> dt
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
saveRDS(hei, file=paste0(parentdir, 'output/hei_chr', args[1], '.Rds'))

#*****
#END *
#*****

