#!/usr/bin/env Rscript

############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
old <- Sys.time()
#args<-c('phys_100000_0.0005', 22, 0.10)
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
###plink betas for AFR
#args<-c('phys_100000_0.0005', 22, 0.1)
chr<-args[2]
final_plink<-as.data.table(readRDS('~/height_prediction/loc_anc_analysis/output/final_plink.Rds'))[CHR==chr]
###read in PRS SNPs for WHI
hei<-readRDS(paste0('~/height_prediction/gwas/WHI/output/hei_', args[1], '_v2.Rds'))[[chr]]
select(hei, -c(i.MarkerName, QUAL, FILTER,INFO, FORMAT, N, SE))-> hei #remove some useless columns
##for each chromosome chr
what <- paste0("~/height_prediction/input/WHI/WHI_b37_strand_include_kgCY_chr", chr)
ind<-read.table(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr, ".phind"), as.is=TRUE)
geno <- read_fwf(paste0('/project/mathilab/data/WHI/data/phased/hapi-ur/WHI_b37_strand_include_kgCY_chr', chr, ".phgeno"), fwf_widths(rep(1,NROW(ind))))
gc()
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

names(geno)<-gsub(':', '_', ind$V1)
names(ancestry)<-gsub(':', '_', ind$V1)


snp[POS %in% hei$POS][POS %in% final_plink$POS]-> snp1
geno<-geno[which(snp1$POS %in% hei$POS),]
ancestry<-ancestry[which(snp1$POS %in% hei$POS),]
#final_plink$PLINK<-scale(final_plink$PLINK)
setkey(final_plink, CHR, POS)
hei1<-hei[MarkerName %in% snp1$i.MarkerName]
setkey(hei1, CHR, POS)
gc()
cat('checkpoint number 2\n')
final_plink[hei1,nomatch=0]-> plink_prun
plink_prun$PLINK<-scale(plink_prun$PLINK)
plink_prun$b<-scale(plink_prun$b) #mnot sure if this is needed
plink_prun %>% select(-contains("0_")) %>% as.data.table-> plink_prun
setDT(geno)
setDT(ancestry)

cat('checkpoint number 3\n')
	
#plink_prun[,AlleleMatch:=ifelse(Effect_Allele_plink==ALT, "TRUE", ifelse(Effect_Allele_plink==REF, "FALSE", NA))]
plink_prun[,AlM_eur:=ifelse(Allele1==ALT, "TRUE", ifelse(Allele1==REF, "FALSE", NA))]
plink_prun[,AlM_afr:=ifelse(Effect_Allele_plink==ALT, "TRUE", ifelse(Effect_Allele_plink==REF, "FALSE", NA))]
LA_PRS<- function(X=geno, X2=10, Y=ancestry, alpha=as.numeric(args[3])){
dataA<-cbind(data.table(Geno=unlist(X[,.SD,.SDcols=X2]), Anc=unlist(Y[,.SD,.SDcols=X2]), plink_prun))
afrA<-dataA[Anc==1]
eurA<-dataA[Anc==2]
if(nrow(afrA)>=1){
afrA[,PRS_part:=ifelse(AlM_afr=='TRUE' & AlM_eur=='TRUE',((alpha*PLINK*Geno)+((1-alpha)*b*Geno)), ifelse(AlM_afr=='TRUE' & AlM_eur=='FALSE',((alpha*PLINK*Geno)+((1-alpha)*b*abs(Geno-1))),ifelse(AlM_afr=='FALSE' & AlM_eur=='TRUE',((alpha*PLINK*abs(1-Geno))+((1-alpha)*b*Geno)),ifelse(AlM_afr=='FALSE' & AlM_eur=='FALSE', (((alpha*PLINK)*abs(1-Geno))+((1-alpha)*b*abs(1-Geno))), NA))))]
}
if(nrow(eurA)>=1){
eurA[, PRS_part:=ifelse(AlM_eur=='TRUE', b*Geno, ifelse(AlM_eur=='FALSE', b*abs(Geno-1),NA))]
}
res<-sum(bind_rows(afrA,eurA)$PRS_part, na.rm=T)
 return(res)
}

system.time(RES<-lapply(1:ncol(geno), function(I) LA_PRS(X2=I))) ### for chr22
names(RES)<-colnames(geno)
lapply(seq(from=1, to=length(RES), by=2), function(I) RES[[I]]+RES[[I+1]])-> RES2  #combine two chr from each individual
names(RES2)<-unique(gsub("_A", "", gsub("_B", "", colnames(geno))))
saveRDS(RES2, file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_', args[3], '_prs.Rds'))
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
