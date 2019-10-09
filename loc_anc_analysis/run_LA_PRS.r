#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
old <- Sys.time()
#args<-c('phys_100000_0.0005', 21)
## Load libraries
library(optparse)
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(parallel)
options(scipen=999)
source('~/height_prediction/scripts/mclapply2.R')
#
#1. WHI data. Genotype and phenotype
#2. Local ancestry info
#3. For each PRS variant, record whether it is Afr or Eur
#3. For EUR, use ukb betas. 
#4. For AFR, use UKB_afr imputed betas
#5. Calculate PRS

###########################################
#args<-'phys_100000_0.0005'
###plink betas for AFR
cat('checkpoint number 1\n')
chr<-args[2]
plink<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/association_v3.Res.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'
plink<-plink[CHR==chr]
#
plink2<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/test3.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","A1","TEST","OBS_CT","PLINK", "SE","T_STAT", "UNADJ") #plink is BETA
plink2<-plink2[CHR==chr]
gc()
setkey(plink, MarkerName, CHR, UNADJ)
setkey(plink2, MarkerName, CHR, UNADJ)
plink[plink2, nomatch=0]-> final_plink
remove(plink, plink2)
gc()
final_plink$POS<-as.numeric(final_plink$POS)
gc()
final_plink$CHR<-as.numeric(final_plink$CHR)
gc()
final_plink[order(CHR,POS)] -> final_plink
select(final_plink, CHR, MarkerName, POS, REF, ALT, A1, PLINK)-> final_plink

###read in PRS SNPs for WHI
hei<-readRDS(paste0('~/height_prediction/gwas/WHI/output/hei_', args[1], '_v2.Rds'))[[chr]]
##for each chromosome chr
what <- paste0("~/height_prediction/input/WHI/WHI_b37_strand_include_kgCY_chr", chr)
ind<-read.table(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr, ".phind"), as.is=TRUE)
geno <- read_fwf(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr, ".phgeno"), fwf_widths(rep(1,NROW(ind))))
ancestry <- read_fwf(paste0(what, "_rfmix_out.0.Viterbi.txt"), fwf_empty(paste0(what, "_rfmix_out.0.Viterbi.txt")))
snp<-fread(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr,  ".phsnp"))
colnames(snp)<-c('i.MarkerName', 'CHR', 'V3', 'POS',  'REF', 'ALT')
#Remove reference samples
samples.to.include <- ind[,3]=="Case"
geno <- geno[,samples.to.include]
ind <- ind[samples.to.include,]

if(!all(dim(geno)==dim(ancestry))){
    stop("Genotype and ancestry matrices different sizes")
}

geno<-geno[which(snp$POS %in% hei$POS),]
ancestry<-ancestry[which(snp$POS %in% hei$POS),]
snp[POS %in% hei$POS]-> snp1
final_plink$PLINK<-scale(final_plink$PLINK)
plink1<-final_plink[CHR==chr]
setkey(plink1, CHR, POS)
hei1<-hei [MarkerName %in% snp1$i.MarkerName]
setkey(hei1, CHR, POS)
gc()
cat('checkpoint number 2\n')
plink1[hei1,]-> plink_prun
plink_prun %>% select(contains("0_")) %>% as.data.table-> plink_prun2

#prepare data
colnames(geno)<-paste0(rep(colnames(plink_prun2),each=2), c("_A", "_B"))
colnames(ancestry)<-paste0(rep(colnames(plink_prun2),each=2), c("_A", "_B"))

#write a function to generate a properly formatted input matrix for the PRS function

Separate<-function(X, X2, col){ #col is the colname
	temp<-X2 %>% separate(col,  paste0(col, c("_A", "_B"))) %>% select(contains(col)) %>% as.data.table
	return(temp)
}

cat('checkpoint number 3\n')
#system.time(input_MATRIX<-mclapply2(names(plink_prun2), function(Y) Separate(X=plink_prun, X2=plink_prun2, col=Y))) #test for 100 samples takes 42 seconds
#saveRDS(input_MATRIX, file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr, '_temp.Rds'))
input_MATRIX<-readRDS(file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr, '_temp.Rds'))
#write a function to calculate the ancestry-specific PRS
input<-do.call(cbind, input_MATRIX)
input2 <- sapply(input, as.numeric)
input2<-as.data.table(as.data.frame(input2))
LA_PRS<- function(X=plink_prun, X2="0_710880_A"){
	dataA<-cbind(select(X, CHR, MarkerName, i.MarkerName, POS, REF, ALT,A1, Allele1, Allele2, b, p, N, PLINK),  Set=unlist(input2[,get(X2)]))
	anc<-select(ancestry, contains(X2))
	afrA<-which(anc[,1]==1)
	eurA<-which(anc[,1]==2)
	dataAF<-dataA[afrA,]
	dataAF[,Set:=as.numeric(dataAF[,Set])]
	dataAF[, PRS_part:=ifelse(Allele1==ALT, Set*PLINK, ifelse(Allele1==REF, abs(1-Set)*PLINK, NA))]
	dataEU<-dataA[eurA,]
	dataEU[,Set:=as.numeric(dataEU[,Set])]
	dataEU[, PRS_part:=ifelse(Allele1==ALT, Set*b, ifelse(Allele1==REF, abs(1-Set)*b, NA))]
	res<-sum(rbind(dataAF,dataEU)$PRS_part, na.rm=T)
	return(res)
}
LA_PRS() #test
colnames(input2)-> I
system.time(mclapply2(1:length(I), function(Z) LA_PRS(X2=I[Z]))-> test) #
saveRDS(test, file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr, '_prs.Rds'))
cat('checkpoint number 4\n')
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
