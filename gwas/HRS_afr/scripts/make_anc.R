#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)

}
#Load packages
library("optparse")
library(data.table)
library(readr)

#add ancestry
what<-paste0("/project/mathilab/data/HRS/data/phased/hapi-ur/HRS_AFR_b37_strand_include_kgCY_chr", args[1]) #Pop1 is AFR
#what <- paste0("~/height_prediction/input/HRS_afr/HRS_AFR_b37_strand_include_kgCY_chr", args[1])  #pop1 is AFR
#ancestry <- read_fwf(paste0(what, "_rfmix_out.0.Viterbi.txt.gz"), fwf_empty(paste0(what, "_rfmix_out.0.Viterbi.txt.gz")))
ancestry<- read_fwf(paste0("/project/mathilab/data/HRS/local_ancestry/RFMix/HRS_AFR_b37_strand_include_kgCY_chr", args[1], "_rfmix_out.0.Viterbi.txt.gz"), fwf_empty(paste0("/project/mathilab/data/HRS/local_ancestry/RFMix/HRS_AFR_b37_strand_include_kgCY_chr", args[1], "_rfmix_out.0.Viterbi.txt.gz")))
gc()

ind<-read.table(paste0(what, ".phind"), as.is=TRUE)
#Remove reference samples
exclude.samples<-grep("\\bNA",ind$V1)
ind <- ind[-exclude.samples,]

indiv<-vector('list', ncol(ancestry))
for(I in 1:ncol(ancestry)){
  #      indiv[[I]]<-(sum(ancestry[[1]][,I]==1)+ sum(ancestry[[2]][,I]==1)+ sum(ancestry[[3]][,I]==1)+sum(ancestry[[4]][,I]==1)+sum(ancestry[[5]][,I]==1)+sum(ancestry[[6]][,I]==1)+sum(ancestry[[7]][,I]==1)
#+sum(ancestry[[8]][,I]==1)+sum(ancestry[[9]][,I]==1)+sum(ancestry[[10]][,I]==1)+sum(ancestry[[11]][,I]==1)+sum(ancestry[[12]][,I]==1)+sum(ancestry[[13]][,I]==1)+sum(ancestry[[14]][,I]==1)+sum(ancestry[[15
#]][,I]==1)+sum(ancestry[[16]][,I]==1)+sum(ancestry[[17]][,I]==1)+sum(ancestry[[18]][,I]==1)+sum(ancestry[[19]][,I]==1)+sum(ancestry[[20]][,I]==1)+sum(ancestry[[21]][,I]==1)+ sum(ancestry[[22]][,I]==1))/su
#m(unlist(lapply(1:22, function(X) nrow(ancestry[[X]]))))
indiv[[I]]<-(sum(ancestry[,I]==1))/nrow(ancestry)
}

indiv2<-vector('list', ncol(ancestry)/2)

i=0
for(I in seq(from=1,to=length(indiv), by=2)){
        i=i+1
        indiv2[[i]]<-sum(indiv[[I]],indiv[[I+1]])/2
        cat(i,"\r")
}

#admixture (obsolete)
#anc_WHI<-cbind(fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.2.Q'), fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.fam')[,V2])
#colnames(anc_WHI)<-c("AFR_ANC","EUR_ANC","SUBJID")
anc_HRS<-data.table(AFR_ANC=unlist(indiv2), EUR_ANC=1-unlist(indiv2), SUBJID=as.integer(unique(gsub("_[A|B]", "",gsub("[0-9]+:","",ind[,1])))), CHR=args[1])
fwrite(anc_HRS, file=paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr',args[1], '.txt'), sep="\t")
