#!/usr/bin/env Rscript
library(data.table)
library(readr)
ancestry<-vector('list', 22)
#add ancestry
for(chr in 1:22){
        what <- paste0("~/height_prediction/input/WHI/WHI_b37_strand_include_kgCY_chr", chr)  #pop1 is AFR
        ancestry[[chr]] <- read_fwf(paste0(what, "_rfmix_out.0.Viterbi.txt"), fwf_empty(paste0(what, "_rfmix_out.0.Viterbi.txt")))
        gc()
}

ind<-read.table(paste0(what, ".phind"), as.is=TRUE)
#Remove reference samples
samples.to.include <- ind[,3]=="Case"
ind <- ind[samples.to.include,]

indiv<-vector('list', ncol(ancestry[[22]]))
for(I in 1:ncol(ancestry[[22]])){
        cat(I,"\r")
        indiv[[I]]<-(sum(ancestry[[1]][,I]==1)+ sum(ancestry[[2]][,I]==1)+ sum(ancestry[[3]][,I]==1)+sum(ancestry[[4]][,I]==1)+sum(ancestry[[5]][,I]==1)+sum(ancestry[[6]][,I]==1)+sum(ancestry[[7]][,I]==1)
+sum(ancestry[[8]][,I]==1)+sum(ancestry[[9]][,I]==1)+sum(ancestry[[10]][,I]==1)+sum(ancestry[[11]][,I]==1)+sum(ancestry[[12]][,I]==1)+sum(ancestry[[13]][,I]==1)+sum(ancestry[[14]][,I]==1)+sum(ancestry[[15
]][,I]==1)+sum(ancestry[[16]][,I]==1)+sum(ancestry[[17]][,I]==1)+sum(ancestry[[18]][,I]==1)+sum(ancestry[[19]][,I]==1)+sum(ancestry[[20]][,I]==1)+sum(ancestry[[21]][,I]==1)+ sum(ancestry[[22]][,I]==1))/su
m(unlist(lapply(1:22, function(X) nrow(ancestry[[X]]))))
}

indiv2<-vector('list', ncol(ancestry[[22]])/2)

i=0
for(I in seq(from=1,to=length(indiv), by=2)){
        i=i+1
        indiv2[[i]]<-sum(indiv[[I]],indiv[[I+1]])/2
        cat(i,"\r")
        }

#admixture (obsolete)
#anc_WHI<-cbind(fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.2.Q'), fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.fam')[,V2])
#colnames(anc_WHI)<-c("AFR_ANC","EUR_ANC","SUBJID")
anc_WHI<-data.table(AFR_ANC=unlist(indiv2), EUR_ANC=1-unlist(indiv2), SUBJID=as.integer(unique(gsub("_[A|B]", "",gsub("0:","",ind[,1])))))
anc_WHI[,SUBJID:=paste0("0_", SUBJID)]
fwrite('~/height_prediction/input/WHI/rfmix_anc.txt')
