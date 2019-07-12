#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#**************************************
#*      CALCULATE POLYGENIC SCORES   **
#**************************************
source('~/height_prediction/scripts/PolygenicScore_v2.R')
library("optparse")
library(data.table)
library(dplyr)
library(biomaRt)
library(parallel)
options(scipen=10)
options(digits=10)


PolScore2<- function(data=vcf2){
        #readRDS(paste0(home, panel, "/", panel2,  '/output/hei_', tag, '_v2.Rds'))-> hei
        data-> hei2
        samps<-colnames(hei2)[14:ncol(hei2)]
        
        hei2[ALT==Allele1]-> temp1
        hei2[REF==Allele1]-> temp2 #im ignoring the other two rows for now
        vector('list', length(samps))-> temp_list
        names(temp_list)<-samps
        cat('Number of samples is', length(samps), '\n')
        if(nrow(temp1)>0 & nrow(temp2)>0){
                matrix(nrow=nrow(temp1)+nrow(temp2), ncol=length(samps))-> my_matrix
                colnames(my_matrix)<-samps
                rownames(my_matrix)<-c(temp1[,MarkerName],temp2[,MarkerName])
                b1<-temp1[,b]
                b2<-temp2[,b]
                counter<-0
	                for(i  in samps){
                        my_matrix[which(temp1[,i, with=F]=="0/0"),i]<-0
                        my_matrix[which(temp1[,i, with=F]=="1/1"),i]<-2
                        my_matrix[which(temp1[,i, with=F]=="1/0"),i]<-1
                        my_matrix[which(temp1[,i, with=F]=="0/1"),i]<-1

                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="0/0"),i]<-2
                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="1/1"),i]<-0
                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="1/0"),i]<-1
                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="0/1"),i]<-1
                        counter<-counter+1
                        cat(counter,'\r')
                }
                 apply(my_matrix*c(b1,b2), 2, sum, na.rm=T)-> res
                cat('Finished first for loop\n')
                } else if (nrow(temp1)>0){
                matrix(nrow=nrow(temp1), ncol=length(samps))-> my_matrix
                colnames(my_matrix)<-samps
                rownames(my_matrix)<-temp1[,MarkerName]
                b1<-temp1[,b]
		              counter<-0
                for(i in samps){
                        my_matrix[which(temp1[,i, with=F]=="0/0"),i]<-0
                        my_matrix[which(temp1[,i, with=F]=="1/1"),i]<-2
                        my_matrix[which(temp1[,i, with=F]=="1/0"),i]<-1
                        my_matrix[which(temp1[,i, with=F]=="0/1"),i]<-1
                }
                 apply(my_matrix*b1, 2, sum, na.rm=T)-> res
        cat('Finished second  for loop\n')
        } else if (nrow(temp2)>0){
                matrix(nrow=nrow(temp2), ncol=length(samps))-> my_matrix
                colnames(my_matrix)<-samps
                rownames(my_matrix)<-temp2[,MarkerName]
                b2<-temp2[,b]
                for(i in samps){

                        my_matrix[which(temp2[,i, with=F]=="0/0"),i]<-2
                        my_matrix[which(temp2[,i, with=F]=="1/1"),i]<-0
                        my_matrix[which(temp2[,i, with=F]=="1/0"),i]<-1
                        my_matrix[which(temp2[,i, with=F]=="0/1"),i]<-1
        }
         apply(my_matrix*b2, 2, sum, na.rm=T)-> res
        }
        cat('Finished third  for loop\n')
#acollect sample names for this population from 1000G data.
        return(res)
}

############
W<-as.numeric(args[1])
P<-as.numeric(args[2])
vec_all<-readRDS(paste0('~/height_prediction/imputed/output/vec_all_',W, '_', P,'.Rds'))


all_prs<-vector('list', 2)
names(all_prs)<-c('HRS_afr', 'HRS_eur')
all_prs[['HRS_afr']]<-vector('list',22)
all_prs[['HRS_eur']]<-vector('list',22)

for(chr in 1:22){
	betas<-vec_all[[chr]]
	vcf<-fread(paste0('zcat output/HRS_afr_', W, '_',P, '_chr_', chr, '.vcf.gz'))
	colnames(vcf)[1]<-"CHR"
	setkey(betas, CHR,POS)
	setkey(vcf, CHR, POS)
	vcf2<-betas[vcf,nomatch=0]
	all_prs[['HRS_afr']][[chr]]<-PolScore2(data=vcf2)
	remove(vcf,vcf2)
	vcf<-fread(paste0('zcat output/HRS_eur_', W, '_', P, '_chr_', chr, '.vcf.gz'))
        colnames(vcf)[1]<-"CHR"
        setkey(betas, CHR,POS)
        setkey(vcf, CHR, POS)
        vcf2<-betas[vcf,nomatch=0]
        all_prs[['HRS_eur']][[chr]]<-PolScore2(data=vcf2)
	cat(chr, ' done\n')
	}

#combine

samps<-names(all_prs[['HRS_afr']][[1]]) #sum PGS across chromosomes.
samps2<-names(all_prs[['HRS_eur']][[1]])
PGS2<-vector('list',2)
names(PGS2)<-names(all_prs)
PGS2[[1]]<-vector('list', length(samps))
PGS2[[2]]<-vector('list', length(samps2))
names(PGS2[[1]])<-samps
names(PGS2[[2]])<-samps2
for (S in samps){
	I<-'HRS_afr'
        sum(all_prs[[I]][[1]][[S]],all_prs[[I]][[2]][[S]],all_prs[[I]][[3]][[S]],all_prs[[I]][[4]][[S]],all_prs[[I]][[5]][[S]],all_prs[[I]][[6]][[S]],all_prs[[I]][[7]][[S]],all_prs[[I]][[8]][[S]],all_prs[[I]][[9]][[S]],all_prs[[I]][[10]][[S]],all_prs[[I]][[11]][[S]],all_prs[[I]][[12]][[S]],all_prs[[I]][[13]][[S]],all_prs[[I]][[14]][[S]],all_prs[[I]][[15]][[S]],all_prs[[I]][[16]][[S]],all_prs[[I]][[17]][[S]],all_prs[[I]][[18]][[S]],all_prs[[I]][[19]][[S]],all_prs[[I]][[20]][[S]],all_prs[[I]][[21]][[S]],all_prs[[I]][[22]][[S]])->PGS2[[I]][[S]]
        cat(paste0(S, ' done\n'))
        }

for(S in samps2){
	I<-'HRS_eur'
        sum(all_prs[[I]][[1]][[S]],all_prs[[I]][[2]][[S]],all_prs[[I]][[3]][[S]],all_prs[[I]][[4]][[S]],all_prs[[I]][[5]][[S]],all_prs[[I]][[6]][[S]],all_prs[[I]][[7]][[S]],all_prs[[I]][[8]][[S]],all_prs[[I]][[9]][[S]],all_prs[[I]][[10]][[S]],all_prs[[I]][[11]][[S]],all_prs[[I]][[12]][[S]],all_prs[[I]][[13]][[S]],all_prs[[I]][[14]][[S]],all_prs[[I]][[15]][[S]],all_prs[[I]][[16]][[S]],all_prs[[I]][[17]][[S]],all_prs[[I]][[18]][[S]],all_prs[[I]][[19]][[S]],all_prs[[I]][[20]][[S]],all_prs[[I]][[21]][[S]],all_prs[[I]][[22]][[S]])->PGS2[[I]][[S]]
        cat(paste0(S, ' done\n'))
        }

saveRDS(all_prs, file=paste0('~/height_prediction/imputed/output/PGS1', W, '_', P, '.Rds'))
saveRDS(PGS2, file=paste0('~/height_prediction/imputed/output/PGS2', W, '_', P, '.Rds'))
