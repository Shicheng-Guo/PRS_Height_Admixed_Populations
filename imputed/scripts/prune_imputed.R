#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
library(data.table)
library(rbgen)
library(dplyr)
options(scipen=10)
options(digits=10)

cat(args[1], '\n')
cat(args[2], '\n')
#read in summ stats file
ukb<-fread('zcat ~/height_prediction/gwas/input/50_raw_filtered.txt.gz', fill=T)
ukb[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
ukb[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
ukb[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
ukb[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb
gc()
#prune
p_thresh<-as.numeric(args[2])
W<-as.numeric(args[1])
vec_all<-vector('list',22)
for(chr in 22:1){
	tmp<-ukb[CHR==chr]
	tmp$POS<-as.numeric(tmp$POS)
	tmp$CHR<-as.numeric(tmp$CHR)
	if(args[3]=='HRS'){
		bim<-fread(paste0('/project/mathilab/data/HRS/data/imputed/HRS_AFR_imputed_chr', chr, '.bim'))
	} else if (args[3]=='UKB'){
		bim<-fread(paste0('/project/mathilab/data/UKB/imputed/ukb_imp_chr', chr,'_afr.bim'))
	}
	colnames(bim)<-c('CHR', 'MarkerName', 'V3','POS', 'REF', 'ALT')
	bim$POS<-as.numeric(bim$POS)
	bim$CHR<-as.numeric(bim$CHR)
        setkey(bim, CHR,POS)
        setkey(tmp, CHR,POS)
	bim[tmp, nomatch=0]-> tmp
	tmp[order(p)]->tmp
	tmp[p<=p_thresh]-> tmp
	tmp[order(p)]-> tmp
	tmp[1,POS]->p1
	tmp[, Dist:=abs(POS-p1)]
	vec2<-data.frame(CHR=NA,POS=NA,Allele1=NA, Allele2=NA, MarkerName=NA, b=NA)
	while(nrow(tmp[p<=p_thresh])>1){
		cat('Starting another round\n')
       		tmp[order(p)]-> tmp
        	ind<-tmp[1,.(CHR,POS,Allele1,Allele2,MarkerName,b)]
        	vec2<-rbind(vec2,ind)
		cat(paste0('Index SNP is ', ind$MarkerName, "\n"))
        	tmp[order(p)]-> tmp
        	tmp[, Dist:=abs(POS-ind[,POS])]
        #physcial window size
        	cat(paste0('Physical window of ', W, 'kb around ', ind$MarkerName, ' removed\n'))
        	tmp[Dist>abs(W)]-> tmp
        	cat(paste0(nrow(tmp), ' SNPs left to check\n'))
        	}
	vec2<-vec2[-1]
	cat('Chr ', chr, ' done\n')
	ukb[CHR!=chr]-> ukb
	if(args[3]=='HRS'){
		fwrite(as.data.frame(vec2[,5]),paste0('~/height_prediction/imputed/output/chr_', chr,'_', W,'_',p_thresh, '.txt'), sep="\t",col.names=F)
	} else if(args[3]== 'UKB'){
		fwrite(as.data.frame(vec2[,5]),paste0('~/height_prediction/imputed/output/UKB_chr_', chr,'_', W,'_',p_thresh, '.txt'), sep="\t",col.names=F)
	}	
	vec_all[[chr]]<-vec2
	remove(vec2)
}
if(args[3]=='HRS'){
	saveRDS(vec_all,file=paste0('~/height_prediction/imputed/output/vec_all_', W, '_', p_thresh, '.Rds'))
} else if(args[3]=='UKB'){
	saveRDS(vec_all,file=paste0('~/height_prediction/imputed/output/UKB_vec_all_', W, '_', p_thresh, '.Rds'))
}
#

