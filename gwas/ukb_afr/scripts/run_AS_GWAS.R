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
#args<-22
what <- paste0("~/height_prediction/gwas/ukb_afr/input/UKB_kgCY_chr", args[1])

#Copy the test files:
#>scp pmacs:/project/mathilab/data/WHI/data/WHI_phenotypes.txt .

system(paste0('cp /project/mathilab/data/UKB/phased/hapi-ur/UKB_kgCY_chr', args[1],'.ph* ~/height_prediction/gwas/ukb_afr/input/'))
system(paste0('cp /project/mathilab/data/UKB/local_ancestry/RFMix/UKB_kgCY_chr', args[1], '_rfmix_out.0.Viterbi.txt.gz ~/height_prediction/gwas/ukb_afr/input/'))

#Load genotype and ancestry data. 
ind<-read.table(paste0(what, ".phind"), as.is=TRUE)
geno <- read_fwf(paste0(what, ".phgeno"), fwf_widths(rep(1,NROW(ind))))
ancestry <- read_fwf(paste0(what, "_rfmix_out.0.Viterbi.txt.gz"), fwf_empty(paste0(what, "_rfmix_out.0.Viterbi.txt.gz")))
snp<-fread(paste0(what, ".phsnp"))
#Remove reference samples #Nov.16: don't remember what this was for but it doesn't apply to this data.
#samples.to.include <- ind[,3]=="Case"
exclude.samples<-grep("\\bNA",ind$V1)
ind <- ind[-exclude.samples,]
setDT(geno)
geno[,colnames(geno)[exclude.samples]:=NULL]
#setDT(ancestry)
#ancestry[, paste0("X", exclude.samples):=NULL]

if(!all(dim(geno)==dim(ancestry))){
    stop("Genotype and ancestry matrices different sizes")
}

## Load and merge phenotypes
pheno <- read.table("UKB_AFR_pheno.txt", as.is=TRUE, header=TRUE)
rownames(pheno) <- pheno$ID
ind <- ind[,c(1,2)]
colnames(ind) <- c("HapID", "Sex")
#ind$SUBJID<-gsub("_[AB]", "", gsub("0:", "", ind[,1]))
ind$SUBJID<-gsub("_[AB]", "",gsub("[0-9]+:", "", ind[,1]))
ind$HEIGHTX <- pheno[as.character(ind$SUBJID),"Height"]
ind$AGE <- pheno[as.character(ind$SUBJID),"Age"]
ind$AGE2 <- ind$AGE^2

N.snps <- dim(geno)[1]
N.haps <- dim(geno)[2]

#Results will be the regression betas 
betas <- data.frame(ALL=rep(0, N.snps), POP1=rep(0,N.snps), POP2=rep(0, N.snps), SNP_ID=snp$V1, POS=snp$V4, CHR=snp$V2) #I think POP1 is AFR

#read in PCs for ukb_afr

#
#
#
#
#
#
#
#
#

#Now iterate over SNPS
for(i in 1:N.snps){
    cat(paste0("\r", i))
    
    #ALL GWAS
    gt <- unlist(geno[i,])
    model.all <- lm(ind$HEIGHTX~gt+ind$Sex+ind$AGE+ind$AGE2)
    betas[i,1] <- summary(model.all)$coefficients[12]

    #Ancestry specific GWAS

    include.1 <- unlist(ancestry[i,])==1
    include.2 <- unlist(ancestry[i,])==2
    model.1 <- lm(ind$HEIGHTX~gt+ind$Sex+ind$AGE+ind$AGE2, subset=include.1)
    model.2 <- lm(ind$HEIGHTX~gt+ind$Sex+ind$AGE+ind$AGE2, subset=include.2)
    betas[i,2] <- summary(model.1)$coefficients[12]
    betas[i,3] <- summary(model.2)$coefficients[12]
}

write.table(betas, paste0("/project/mathilab/bbita/gwas_admix/new_height/ukb_afr/AS_Beta_chr",  args[1],"example.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

