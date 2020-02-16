#!/usr/bin/env Rscript
library(data.table)
library(stats)
fread('My_Pheno.txt')-> pheno
fread('tmp')-> tmp


setkey(tmp, SUBJID)
pheno[,SUBJID:=paste0(ID,"_", ID)]
setkey(pheno, SUBJID)
pheno[,ID:=NULL]
cbind(FID=pheno$SUBJID, pheno)[,.(FID,SUBJID,Sex, Height,Age)]-> pheno
pheno[tmp, nomatch=0][,Age2:=Age^2]-> final


colnames(final)[1]<-"FID"
colnames(final)[2]<-"IID"

final[complete.cases(final),]-> final2
final2[, Res.Height:=summary(lm(final2$Height~final2$Age))$residuals]
fwrite(final2, file="Pheno.txt", sep="\t")
