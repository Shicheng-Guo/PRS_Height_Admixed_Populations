#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
#Load packages
library("optparse")
library(data.table)
library(dplyr)
library(biomaRt)
#**************************************************************************
#read in gwas data and list of 697 SNPs highlighted in the wood paper *****
#**************************************************************************
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
#************************************************************
#using bcftools, get only relevant positions in vcf files
#************************************************************

print('checkpoint 1')

paste0('mkdir -m 777 ', parentdir, args[1], '/')-> smtg
#try(system(smtg))
#for(X in 1:22){ #generates temp files
#	system(paste0('touch ', parentdir, args[1],"/temp_chr", X, ".txt"))
#        ukb_height[CHR==X]-> a
#        setkey(a, POS)
#        unique(a)-> a
#        fwrite(as.data.frame(a[, .(CHR,POS)][order(CHR,POS)][,V3:=paste(CHR,POS,sep="\t")][,V3]), file=paste0(parentdir, args[1],"/temp_chr", X, ".txt"), quote=F, col.names=F, row.names=F)
#        print(X)
#        paste0(parentdir, args[1],"/temp_chr", X, ".txt")-> b
#        system(paste0('sort ', b, " |uniq > ", parentdir, args[1], "/tmp"))
#        system(paste0('mv ', parentdir, args[1], '/tmp ',b ))
#}

#system(paste0(home,'scripts/master.sh ', parentdir, args[1], " ", args[3])) #generates temp2 files
#print('so far so good')
#Sys.sleep(5)
#system(paste0(home,'scripts/combine_temp.bash ', parentdir, " ", args[3])) #generates vcf files
#Sys.sleep(60)
this<-paste0('for K in {12..22}; do bsub -M 40000 -o ', parentdir, 'logs/logthis -e ' , parentdir, 'logs/logthis Rscript --vanilla ', home, 'scripts/make_vcf_v2.R $K ', args[2], ' ' , args[3], '; done') #generates .Rds files
system(this)
andthis<-paste0('for K in {3..11}; do bsub -M 90240 -o ', parentdir, 'logs/logandthis -e ', parentdir, 'logs/logandthis Rscript --vanilla ', home, 'scripts/make_vcf_v2.R $K ', args[2], ' ' , args[3], '; done') #generates .Rds files
system(andthis)
andthis2<-paste0('for K in {1..2}; do bsub -M 90240 -o ', parentdir, 'logs/logandthis2 -e ', parentdir, 'logs/logandthis2 Rscript --vanilla ', home,'scripts/make_vcf_v2.R $K ', args[2], ' ' , args[3], '; done') #generates .Rds files
system(andthis2)
#******
#END **
#******
