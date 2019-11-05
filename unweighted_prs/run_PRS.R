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
#args<-c("LD_prun","LD_250000_0.01_0.5")
cat(args)
options(scipen=999)
#options(digits=10)
home='~/height_prediction/'

hei<-readRDS(file=paste0(home, args[1], '/', args[2],'/output/hei_', args[3], '_v2.Rds'))
lapply(1:22, function(X) hei[[X]][,b:=ifelse(b>0, 1, ifelse(b<0, -1,0))]) #this line is crucial for the unweighted PRS. Everything else is the same as with the weighted one.

cat('checkpoint 3\n')
PGS<-vector('list',22)
names(PGS)<-c(1:22)
for (CR in 1:22){
        print(paste0("Chromosome is ", CR))
        try(PolScore2(CHR=CR, panel=args[1], panel2=args[2], tag=args[3]))-> PGS[[CR]]
        cat(paste0(CR, '  done\n'))
}

samps<-names(PGS[[1]]) #sum PGS across chromosomes.
PGS2<-vector('list', length(samps))
names(PGS2)<-samps
for (S in samps){
        sum(PGS[[1]][[S]],PGS[[2]][[S]],PGS[[3]][[S]],PGS[[4]][[S]],PGS[[5]][[S]],PGS[[6]][[S]],PGS[[7]][[S]],PGS[[8]][[S]],PGS[[9]][[S]],PGS[[10]][[S]],PGS[[11]][[S]],PGS[[12]][[S]],PGS[[13]][[S]],PGS[[14]][[S]],PGS[[15]][[S]],PGS[[16]][[S]],PGS[[17]][[S]],PGS[[18]][[S]],PGS[[19]][[S]],PGS[[20]][[S]],PGS[[21]][[S]],PGS[[22]][[S]])->PGS2[[S]]
        cat(paste0(S, ' done\n'))
        }
saveRDS(PGS2, file=paste0(home, 'unweighted_prs/output/PGS_', args[2], '_', args[3], '.Rds'))

#TheEnd
