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
#args<-c("ALL","F")

fread('/project/mathilab/bbita/gwas_admix/height_prediction/sib_betas/input/sib_beta_gwas_p.txt')-> ukb_height #read in GWAS summary statistics for height from the Uk Biobank
ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]

ukb_height[order(CHR,POS)]-> ukb_height

#************************************************************
#using bcftools, get only relevant positions in vcf files
#************************************************************

print('checkpoint 1')

#if skip is true, this will be skipped:
dir<-paste0(getwd(), "/")
paste0('mkdir -m 777 ', dir, args[1], '/')-> smtg
try(system(smtg))
for(X in 1:22){
	system(paste0('touch ', dir, args[1],"/temp_chr", X, ".txt"))
        ukb_height[CHR==X]-> a
        setkey(a, POS)
        unique(a)-> a
        fwrite(as.data.frame(a[, .(CHR,POS)][order(CHR,POS)][,V3:=paste(CHR,POS,sep="\t")][,V3]), file=paste0(dir, args[1],"/temp_chr", X, ".txt"), quote=F, col.names=F, row.names=F)
        print(X)
        paste0(dir, args[1],"/temp_chr", X, ".txt")-> b
        system(paste0('sort ', b, " |uniq > ", dir, args[1], "/tmp"))
        system(paste0('mv ', dir, args[1], '/tmp ',b ))
}

system(paste0('bash ', dir, 'scripts/master_sib.sh ', dir, args[1]))
print('so far so good')
Sys.sleep(10)
system(paste0('bash ', dir, 'scripts/combine_temp_sib.bash ', dir))
Sys.sleep(60)
this<-paste0('for K in {12..22}; do bsub -o ', dir, 'logs/logthis -e' , dir, 'logs/logthis -M 12000 Rscript --vanilla ', dir, 'scripts/giant_all_v2.R $K; done')
system(this)
andthis<-paste0('for K in {3..11}; do bsub -M 10240 -o ', dir, 'logs/logandthis -e ', dir, 'logs/logandthis Rscript --vanilla ', dir, 'scripts/giant_all_v2.R $K; done')
system(andthis)
andthis2<-paste0('for K in {1..2}; do bsub -M 15240 -o ', dir, 'logs/logandthis2 -e ', dir, 'logs/logandthis2 Rscript --vanilla ',dir,'scripts/giant_all_v2.R $K; done')
system(andthis2)
#******
#END **
#******
