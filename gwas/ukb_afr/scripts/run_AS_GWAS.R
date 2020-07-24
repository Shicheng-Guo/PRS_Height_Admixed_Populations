#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
#Load packages
library("optparse")
## Rough outline of ancestry-specific GWAS analysis.
library(readr)
library(data.table)
library(dplyr)
options(scipen=999)
#args<-22
#what <- paste0("~/height_prediction/gwas/ukb_afr/input/UKB_kgCY_chr", args[1])
what<-paste0("/project/mathilab/data/UKB/phased/hapi-ur/UKB_kgCY_chr", args[1])
#Copy the test files:
#>scp pmacs:/project/mathilab/data/WHI/data/WHI_phenotypes.txt .

#system(paste0('cp /project/mathilab/data/UKB/phased/hapi-ur/UKB_kgCY_chr', args[1],'.ph* ~/height_prediction/input/ukb_afr/'))
#system(paste0('cp /project/mathilab/data/UKB/local_ancestry/RFMix/UKB_kgCY_chr', args[1], '_rfmix_out.0.Viterbi.txt.gz ~/height_prediction/input/ukb_afr/'))

#Load genotype and ancestry data. 
ind<-read.table(paste0(what, ".phind"), as.is=TRUE)
geno <- read_fwf(paste0(what, ".phgeno"), fwf_widths(rep(1,NROW(ind))))
ancestry <- read_fwf(paste0("/project/mathilab/data/UKB/local_ancestry/RFMix/UKB_kgCY_chr", args[1], "_rfmix_out.0.Viterbi.txt.gz"), fwf_empty(paste0("/project/mathilab/data/UKB/local_ancestry/RFMix/UKB_kgCY_chr", args[1], "_rfmix_out.0.Viterbi.txt.gz")))
snp<-fread(paste0(what, ".phsnp"))
#Remove reference samples #Nov.16: don't remember what this was for but it doesn't apply to this data.
#samples.to.include <- ind[,3]=="Case"
exclude.samples<-grep("\\bNA",ind$V1)
ind <- ind[-exclude.samples,]
geno<-geno[,-exclude.samples]
#setDT(ancestry)
#ancestry[, paste0("X", exclude.samples):=NULL]

if(!all(dim(geno)==dim(ancestry))){
    stop("Genotype and ancestry matrices different sizes")
}


## Load and merge phenotypes
pheno<-read.table("~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt", as.is=TRUE, header=TRUE)
rownames(pheno) <- pheno$ID
ind <- ind[,c(1,2)]
colnames(ind) <- c("HapID", "Sex")
#ind$SUBJID<-gsub("_[AB]", "", gsub("0:", "", ind[,1]))
ind$SUBJID<-gsub("_[AB]", "",gsub("[0-9]+:", "", ind[,1]))
ind$HEIGHTX <- pheno[as.character(ind$SUBJID),"Height"]
ind$AGE <- pheno[as.character(ind$SUBJID),"Age"]
ind$AGE2 <- ind$AGE^2
setDT(ind)


#read in PCs for ukb_afr
pcs_A<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/tmp')
pcs_A<-pcs_A[,.(SUBJID, PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
pcs_A<-arrange(pcs_A, SUBJID) %>% as.data.table
pcs_A[,SUBJID:=gsub("[0-9]+_", "", pcs_A$SUBJID)]
pcs_A<-pcs_A[, HapID:=paste0(SUBJID, ":", SUBJID, "_A")]
pcs_B<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/tmp')
pcs_B<-pcs_B[,.(SUBJID, PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
pcs_B<-arrange(pcs_B, SUBJID) %>% as.data.table
pcs_B[,SUBJID:=gsub("[0-9]+_", "", pcs_B$SUBJID)]
pcs_B<-pcs_B[, HapID:=paste0(SUBJID, ":", SUBJID, "_B")]
pcs<-setDT(arrange(rbind(pcs_A,pcs_B), SUBJID))
setkey(pcs, HapID, SUBJID)
pcs<-pcs[which(pcs$SUBJID %in% ind$SUBJID),]# remove 4 liens that don't match
setkey(ind, HapID, SUBJID)
setkey(pcs, HapID, SUBJID)
ind<-ind[pcs]

lapply(1:ncol(ancestry), function(X) table(ancestry[,X]))-> test
table(unlist(lapply(test, function(X) X[1]>X[2])))  #POP1 is AFR
#
#FALSE  TRUE
# 1106  3604
geno[,-which(is.na(ind$HEIGHTX))]-> geno #tremove people without height info
ancestry[,-which(is.na(ind$HEIGHTX))]-> ancestry
ind[-which(is.na(ind$HEIGHTX)),]-> ind

N.snps <- dim(geno)[1]
N.haps <- dim(geno)[2]

#Results will be the regression betas 
betas <- data.frame(ALL=rep(0, N.snps), POP1=rep(0,N.snps), POP2=rep(0, N.snps), SNP_ID=snp$V1, POS=snp$V4, CHR=snp$V2, Tstat_all=rep(NA,N.snps), Tstat_pop1=rep(NA,N.snps), Tstat_pop2=rep(NA,N.snps), SE_all=rep(NA, N.snps), SE_pop1=rep(NA, N.snps), SE_pop2=rep(NA,N.snps)) #I think POP1 is AFR

test<-lapply(1:ncol(ancestry), function(X) table(ancestry[,X]))
do.call(rbind,test[which(unlist(lapply(test, function(X) length(X)==2)))])-> test2

for(i in 1:N.snps){
    cat(paste0("\r", i))
    #ALL GWAS
    gt <- unlist(geno[i,])
    if(sum(gt)==0 ||sum(gt)==length(gt)){ #if all genotypes are the same, the mdoel won't produce a beta for gt
	betas[i,1]<-NA
	betas[i,7]<-NA
	} else{	 
    		model.all <- lm(ind$HEIGHTX~gt+ind$Sex+ind$AGE+ind$AGE2+ind$PC1+ind$PC2+ind$PC3+ind$PC4+ind$PC5+ind$PC6+ind$PC6+ind$PC7+ind$PC8+ind$PC9+ind$PC10)
    		betas[i,1] <- summary(model.all)$coefficients['gt','Estimate']
    		#Ancestry specific GWAS
	
		betas[i,7]<-coef(summary(model.all))[, "t value"][['gt']]
	}
	include.1 <- unlist(ancestry[i,])==1
    	include.2 <- unlist(ancestry[i,])==2
	if(sum(gt[include.1])==0 || sum(gt[include.1])==length(gt[include.1])){
	       betas[i,2]<-NA
	       betas[i,8]<-NA
	} else {
		model.1 <- lm(ind$HEIGHTX~gt+ind$Sex+ind$AGE+ind$AGE2+ind$PC1+ind$PC2+ind$PC3+ind$PC4+ind$PC5+ind$PC6+ind$PC6+ind$PC7+ind$PC8+ind$PC9+ind$PC10, subset=include.1)
		betas[i,2] <- summary(model.1)$coefficients['gt','Estimate']
		betas[i,8]<-coef(summary(model.1))[, "t value"][['gt']]
	}
	if(sum(gt[include.2])==0 || sum(gt[include.2])==length(gt[include.2])){
        	betas[i,3]<-NA
		betas[i,9]<-NA
	} else {
    		model.2 <- lm(ind$HEIGHTX~gt+ind$Sex+ind$AGE+ind$AGE2+ind$PC1+ind$PC2+ind$PC3+ind$PC4+ind$PC5+ind$PC6+ind$PC6+ind$PC7+ind$PC8+ind$PC9+ind$PC10, subset=include.2)
    		betas[i,3] <- summary(model.2)$coefficients['gt','Estimate']
		betas[i,9]<-coef(summary(model.2))[, "t value"][['gt']]
	}
}
write.table(betas, paste0("~/height_prediction/gwas/ukb_afr/output/AS_Beta_chr",  args[1],"example.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)


