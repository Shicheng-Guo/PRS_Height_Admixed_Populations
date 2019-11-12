#!/usr/bin/env Rscript
library(data.table)
library(stats)
fread('My_Pheno.txt')-> pheno
#fread('tmp')-> tmp

#setkey(tmp, SUBJID)
pheno[,SUBJID:=ID]
setkey(pheno, SUBJID)
pheno[,ID:=NULL]
cbind(FID=pheno$SUBJID, pheno)[,.(FID,SUBJID,Sex, Height,Age)]-> pheno
#pheno[tmp, nomatch=0][,Age2:=Age^2]-> final


colnames(pheno)[1]<-"FID"
colnames(pheno)[2]<-"IID"

pheno[complete.cases(pheno),]-> final
#final2[, Res.Height:=summary(lm(final2$Height~final2$Age))$residuals]
fwrite(final, file="Pheno.txt", sep="\t")
