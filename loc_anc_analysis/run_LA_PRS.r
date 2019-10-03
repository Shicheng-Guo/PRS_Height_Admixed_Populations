#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}

## Load libraries
library(optparse)
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
options(scipen=999)
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
plink<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/association_v3.Res.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'
N<-8816*2
#
plink2<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR_imputed/test3.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","A1","TEST","OBS_CT","PLINK", "SE","T_STAT", "UNADJ") #plink is BETA
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
#chr<-22
chr<-args[2]
hei<-readRDS(paste0('~/height_prediction/gwas/WHI/output/hei_', args[1], '_v2.Rds'))[[chr]]
##for each chromosome chr
#for(chr in 1:22){
what <- paste0("~/height_prediction/input/WHI/WHI_b37_strand_include_kgCY_chr", chr)
ind<-read.table(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr, ".phind"), as.is=TRUE)
geno <- read_fwf(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr, ".phgeno"), fwf_widths(rep(1,NROW(ind))))
ancestry <- read_fwf(paste0(what, "_rfmix_out.0.Viterbi.txt"), fwf_empty(paste0(what, "_rfmix_out.0.Viterbi.txt")))
snp<-fread(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr,  ".phsnp"))
colnames(snp)<-c('i.MarkerName', 'CHR', 'V3', 'POS',  'REF', 'ALT')
#Remove reference samples
samples.to.include <- ind[,3]=="Case"
geno <- geno[,samples.to.include]
#geno<-as.data.frame(geno)[,samples.to.include]
#setDT(geno)
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

#for each individual
res_all<-vector('list', length(colnames(plink_prun2)))
names(res_all)<-colnames(plink_prun2)
#two haplotypes per individuals
counter<-0
for(i in colnames(plink_prun2)){
#i=colnames(plink_prun2)[1]
res_all[[i]]<-vector('list', 2)
names(res_all[[i]])<-c("A","B")
i2=paste0(i, '_A')
colnames(geno)<-paste0(rep(colnames(plink_prun2),each=2), c("_A", "_B"))
colnames(ancestry)<-paste0(rep(colnames(plink_prun2),each=2), c("_A", "_B"))
temp<-plink_prun2 %>% separate(i,  paste0(i, c("_A", "_B"))) %>% select(contains(i)) %>% as.data.table
data<-cbind(select(plink_prun, CHR, MarkerName, i.MarkerName, POS, REF, ALT,A1, Allele1, Allele2, b, p, N, PLINK),  Set=temp[,get(i2)]) #the real ALT (or ALLELE 1, always) is the effect allele. In plink, A1 is the one for which beta is estimated, which tends to be the reference allele. REF ALT here mean nothing, they came from WHI.
#assing positions to AFR or EUR
haps<-select(geno, contains(i)) #why does this not match temp? nO idea.
anc<-select(ancestry, contains(i))
afrA<-which(anc[,1]==1)
eurA<-which(anc[,1]==2)
afrB<-which(anc[,2]==1)
eurB<-which(anc[,2]==2)
#
data1<-data[afrA,]
data1[,Set:=as.numeric(data1[,Set])]
data1[, PRS_part:=ifelse(Allele1==ALT, Set*PLINK, ifelse(Allele1==REF, abs(1-Set)*PLINK, NA))]
data2<-data[eurA,]
data2[,Set:=as.numeric(data2[,Set])]
data2[, PRS_part:=ifelse(Allele1==ALT, Set*b, ifelse(Allele1==REF, abs(1-Set)*b, NA))]
res_all[[i]][['A']]<-sum(rbind(data1,data2)$PRS_part)
#hap B
i2=paste0(i, '_B')
data<-cbind(select(plink_prun, CHR, MarkerName, i.MarkerName, POS, REF, ALT,A1, Allele1, Allele2, b, p, N, PLINK),  Set=temp[,get(i2)])
data1<-data[afrB,]
data1[,Set:=as.numeric(data1[,Set])]
data1[, PRS_part:=ifelse(Allele1==ALT, Set*PLINK, ifelse(Allele1==REF, abs(1-Set)*PLINK, NA))]
data2<-data[eurB,]
data2[,Set:=as.numeric(data2[,Set])]
data2[, PRS_part:=ifelse(Allele1==ALT, Set*b, ifelse(Allele1==REF, abs(1-Set)*b, NA))]
res_all[[i]][['B']]<-sum(rbind(data1,data2)$PRS_part)
#cat(i, "done\n")
counter<-counter+1
cat(counter, '\n')
}

save(res_all, file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr, '.Rds')
