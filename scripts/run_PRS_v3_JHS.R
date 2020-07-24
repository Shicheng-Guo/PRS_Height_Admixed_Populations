#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
old <- Sys.time()
## Load libraries
library(optparse)
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(parallel)
options(scipen=999)
#Overview
#1. WHI data: Genotype and phenotype, local ancestry (per haplotype)
#2. Merge with UKB_afr (imputed) GWAS output. Effect alleles for UKB_AFR gwas are the ALT allele. The same is true for the UKB GWAS.
#2. For each PRS variant, record whether it is Afr or Eur
#3. For EUR, use ukb betas. 
#4. For AFR, use UKB_afr imputed betas weighted by alpha + UKB beta weighted by 1-alpha
#5. Calculate PRS for each haplotype
#6. Add PRS from two haplotypes that compose one individual. 
#7. This script runs each chromosome separely. The script (run_partial_r2.r) combines these outputs and evaluates prediction power
###########################################
#args<-c('phys_100000_0.0005', 20) # arguments should be: pruned set, chromosome#comment this line
chr<-args[2]
###read in PRS SNPs for WHI
hei<-readRDS(paste0('~/height_prediction/gwas/JHS/output/hei_', args[1], '_v2.Rds'))[[chr]]
select(hei, -c(QUAL, FILTER,INFO, FORMAT, N, SE))-> hei #remove some useless columns
gc()
colnames(hei)[3:4]<-c('REF_hei', 'ALT_hei') #these are not true ref/alt columns. They are called this way because of the plink to vcf conversion which does not always preserve ref/alt states.
##for each chromosome, read in .phgeno, .phsnp, .phind, and local ancestry files:
what <- paste0("/project/mathilab/data/JHS/local_ancestry/RFMix/JHS_b37_strand_kgCY_chr", chr)
ind<-read.table(paste0('/project/mathilab/data/JHS/data/phased/hapi-ur/JHS_b37_strand_kgCY_chr', chr, ".phind"), as.is=TRUE)
geno <- read_fwf(paste0('/project/mathilab/data/JHS/data/phased/hapi-ur/JHS_b37_strand_kgCY_chr', chr, ".phgeno"), fwf_widths(rep(1,NROW(ind))))
gc()
ancestry <- read_fwf(paste0(what, "_rfmix_out.0.Viterbi.txt.gz"), fwf_empty(paste0(what, "_rfmix_out.0.Viterbi.txt.gz")))
gc()
snp<-fread(paste0('/project/mathilab/data/JHS/data/phased/hapi-ur/JHS_b37_strand_kgCY_chr', chr,  ".phsnp"))
colnames(snp)<-c('MarkerName', 'CHR', 'V3', 'POS',  'REF_snp', 'ALT_snp') #these are not true ref/alt columns.
gc()
#Remove reference samples from snp and geno files
samples.to.include <- ind[,3]=="Unknown"
geno <- geno[,samples.to.include]
ind <- ind[samples.to.include,]

if(!all(dim(geno)==dim(ancestry))){
    stop("Genotype and ancestry matrices different sizes")
}

#name columns of geno and ancestry datasets (individual subject IDs)
names(geno)<-gsub(':', '_', ind$V1)
names(ancestry)<-gsub(':', '_', ind$V1)

#need to make sure all datasets have same set of SNPs (and in the same order)
hei<-unique(hei, by='POS')
snp1<-snp[which(snp$POS %in% hei$POS),] #not all SNPs from the pruned sets are present in the .phsnp files.
geno1<-geno[which(snp$POS %in% hei$POS),]
anc1<-ancestry[which(snp$POS %in% hei$POS),]
#
colnames(hei)[2]<-'MarkerName'
hei %>% select(-contains("0_")) %>% as.data.table-> hei #Allele1 is always the ALT allele. If we check hei[ALT_hei==Allele1] we will see that that's not the case for all positions, showing that ref_hei/alt_hei are not reliable for ref/alt info.
merge(snp1, hei, by=c('CHR', 'POS', 'MarkerName'), sort=F)->hei1 #another test here: hei1[REF_snp==ALT_hei] all lines are true
gc()
cat('checkpoint number 2\n')
if(nrow(snp1)!=nrow(hei1)){
	snp2<-snp1[-which(!(snp1$POS %in% hei1$POS)),]
	geno2<-geno1[-which(!(snp1$POS %in% hei1$POS)),]
	anc2<-anc1[-which(!(snp1$POS %in% hei1$POS)),]
	hei2<-hei1[-which(!(snp1$POS %in% hei1$POS)),]
} else{
	snp2<-snp1
        geno2<-geno1
        hei2<-hei1
        anc2<-anc1	
}

cat('checkpoint number 3\n')

#
#in 'hei' Allele1 is the effect allele, which is always ALT. The columns ALT and REF mean nothing. 
#in geno, 1=REF, 0=ALT
hei2[, AlM:=ifelse(Allele1==ALT_snp, "YES", "NO")] #if YES, 0=REF and 1=ALT. If NO, 1=REF, 0=ALT, thus 1-Geno
#in the geno file 0 is ALT and 1 is REF
LA_PRS<- function(X=geno2, X2=1, Y=anc2){
dataA<-cbind(data.table(Anc=unlist(Y[,X2]), hei2))[, Geno:=ifelse(AlM=='YES', (1-unlist(X[,X2])), unlist(X[,X2]))][,PRS_part:=b*Geno] ##if YES, 0=REF and 1=ALT. If NO, 1=REF, 0=ALT, thus 1-Geno
eurA<-dataA[Anc==2]
res<-sum(eurA$PRS_part, na.rm=T)
return(res)
}

system.time(RES<-lapply(1:ncol(geno2), function(I) LA_PRS(X2=I))) ### 
names(RES)<-colnames(geno)
lapply(seq(from=1, to=length(RES), by=2), function(I) (RES[[I]]+RES[[I+1]]))-> RES2  #combine two chr from each individual
names(RES2)<-unique(gsub("_A", "", gsub("_B", "", colnames(geno2))))
saveRDS(RES2, file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_prs_JHS_EuR.Rds'))
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


