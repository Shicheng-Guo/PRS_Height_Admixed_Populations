#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
#Load packages #########################
library("optparse")
library(data.table)
library(dplyr)
library(asbio)
########################################
home<-"/home/bbita"
dir<-"height_prediction/strat_prs/scripts"
source(paste0(home, "/", dir,'/PRS_calc.R'))
#partial R2 function
source(paste0(home, "/", dir,'/Rsq_R2.R'))
########################################
#####################################################################################################
#### First part: calculate recombination rates for PRS SNPs and divide them into quantiles ##########
#####################################################################################################
## Read in betas and recombination maps
rate.dist<-as.numeric(args[1])
w_map<-args[2]
PRS<-vector('list', 22)
#afr<-fread('/project/mathilab/bbita/gwas_admix/new_height/ukb_afr_betas_100000_0.0005.txt') #need to fix this path
#CHR,POS,MarkerName,i.MarkerName,REF,ALT,Allele1,Allele2,b,SE,p,N
afr<-do.call(rbind, readRDS(paste0('~/height_prediction/', args[3],'/HRS_afr/output/hei_phys_100000_0.0005_v2.Rds')))[,.(CHR,POS,MarkerName,i.MarkerName,REF,ALT,Allele1,Allele2,b,p)]
lapply(1:22, function(chr) fread(paste0('zcat /project/mathilab/data/maps/hm2/hm2/genetic_map_GRCh37_chr', chr,'.txt.gz'))[,CHR:=gsub("chr","",Chromosome)][, Chromosome:=NULL])-> rec #need to fix this path
for(chr in 1:22){colnames(rec[[chr]])<-c('POS', 'RATE_cM_Mb', 'MAP_cM', 'CHR')}
lapply(1:22, function(chr) fread(paste0('zcat /project/mathilab/data/maps_b37/maps_chr.', chr, '.gz')))-> maps #need to fix this path
for(chr in 1:22){colnames(maps[[chr]])[1]<-"POS"}
lapply(1:22, function(chr) setkey(rec[[chr]],POS))
lapply(1:22, function(chr) setkey(maps[[chr]],POS))
lapply(1:22, function(chr) maps[[chr]][rec[[chr]], nomatch=0])-> map
lapply(1:22, function(chr) map[[chr]][,POS2:=POS])
remove(maps)
cat('loading done\n')
chrs<-vector('list',22)
for (chr in 22:1){
        cat('starting chr', chr, '\n')
        a<-seq(from=map[[chr]]$POS[1]-(rate.dist/2), to=map[[chr]]$POS[nrow(map[[chr]])], by=rate.dist)
        b<-c(a[-1],a[length(a)]+rate.dist)
        if(w_map=="AA"){
                AA.rate <- approxfun(map[[chr]]$POS, map[[chr]]$AA_Map, rule=2)
                e<-data.table(POS=a, POS2=b, Win=paste0(a,"|",b))
                setkey(e, POS, POS2)
                setkey(map[[chr]], POS, POS2)
                f<-foverlaps(map[[chr]], e)
                f[, diff:=AA.rate(f$POS2)-AA.rate(f$POS)]
                chrs[[chr]]<-f
                } else if (w_map=="CEU"){
                CEU.rate<-approxfun(map[[chr]]$POS, map[[chr]]$CEU_LD, rule=2)
                e<-data.table(POS=a, POS2=b, Win=paste0(a,"|",b))
                setkey(e, POS, POS2)
                setkey(map[[chr]], POS, POS2)
                f<-foverlaps(map[[chr]], e)
                f[, diff:=CEU.rate(f$POS2)-CEU.rate(f$POS)]
                chrs[[chr]]<-f
                }
}
do.call(rbind, chrs)-> f  #combine into one data.table with all chromosomes
f[,Quantile:=cut(diff, breaks=quantile(diff), na.rm=T, include.lowest=T)]
afr[, POS1:=POS-(rate.dist/2)][,POS2:=POS+(rate.dist/2)]

#Stratify genome into 4 quantiles of recombination rate
if(w_map=="CEU"){
        afr[, diff:=CEU.rate(POS2)-CEU.rate(POS1), by=CHR]
        } else if (w_map=="AA"){
         afr[, diff:=AA.rate(POS2)-AA.rate(POS1), by=CHR]
}

#Split PRS SNPs according to these quantiles
#lev<-quantile(f$diff)
#q1<-afr[diff>=lev[1] & diff<lev[2]]
#q2<-afr[diff>=lev[2] & diff<lev[3]]
#q3<-afr[diff>=lev[3] & diff<lev[4]]
#q4<-afr[diff>=lev[4]]
#cat('checkpoint\n')
afr[,Quantile:=cut(diff, breaks=quantile(diff), na.rm=T, include.lowest=T)]
q1<-split(afr, by='Quantile')[[1]]
q2<-split(afr, by='Quantile')[[2]]
q3<-split(afr, by='Quantile')[[3]]
q4<-split(afr, by='Quantile')[[4]]
#calculate PRS for each of these quantiles
hei<-lapply(1:22, function(chr) readRDS(paste0('~/height_prediction/', args[3], '/HRS_afr/output/hei_phys_100000_0.0005_v2.Rds'))[[chr]]) 
hei<-do.call(rbind, hei)
hei[MarkerName %in% q1$MarkerName]-> hei_q1
hei[MarkerName %in% q2$MarkerName]-> hei_q2
hei[MarkerName %in% q3$MarkerName]-> hei_q3
hei[MarkerName %in% q4$MarkerName]-> hei_q4
prs<-vector('list', 4)

names(prs)<-c("q1","q2","q3", "q4")
prs[['q1']]<-PolScore2(hei2=hei_q1)
cat('q1 done\n')
prs[['q2']]<-PolScore2(hei2=hei_q2)
cat('q2 done\n')
prs[['q3']]<-PolScore2(hei2=hei_q3)
cat('q3 done\n')
prs[['q4']]<-PolScore2(hei2=hei_q4)
cat('q4 done\n')
PRS<-prs
remove(prs)
saveRDS(PRS,file=paste0("~/height_prediction/strat_prs/output/prs_HRS_afr_", args[3], '_', rate.dist, "_", w_map, "_v2.Rds")) #store results
obj<-c(nrow(hei_q1), nrow(hei_q2), nrow(hei_q3), nrow(hei_q4))
saveRDS(obj, file=paste0("~/height_prediction/strat_prs/output/Nr_SNPs_HRS_afr_",args[3], '_', rate.dist, "_", w_map, "_v2.Rds")) #store results

#Make a list
PRS2<-vector('list', length(PRS[['q1']]))
names(PRS2)<-names(PRS[['q1']])
for(J in names(PRS2)){
        PRS2[[J]]<-data.table(q1=PRS[['q1']][[J]],q2=PRS[['q2']][[J]],q3=PRS[['q3']][[J]], q4=PRS[['q4']][[J]])
        cat(J, '\r')
}
do.call(rbind,PRS2)-> PRS2 #combine into one data.table
rownames(PRS2)<-names(PRS[[1]])
saveRDS(PRS2, file=paste0('~/height_prediction/strat_prs/output/PRS2_HRS_afr_',args[3], '_',rate.dist,"_", w_map, '_v2.Rds')) #store results


#########################################################
#### Second part: calculate partial R2 for PRS. #########
#########################################################

#read in phenotype data
fread('~/height_prediction/input/HRS_afr/HRS_AFR_phenotypes.txt', fill=T)[,V5:=NULL]-> Pheno_HRS_afr
#a partial R2 function
source('~/height_prediction/strat_prs/scripts/Rsq_R2.R')

##############

#add PGS to Pheno table in order to be able to make multiple analyses
#Pheno_HRS_afr[,ID:=paste0(ID, "_", ID)]
as.character(Pheno_HRS_afr$ID)-> Pheno_HRS_afr$ID
Pheno_HRS_afr$ID<-paste0(Pheno_HRS_afr[,ID],"_",Pheno_HRS_afr[,ID])
setkey(Pheno_HRS_afr, ID)

#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/HRS_afr/rfmix_anc_chr', X, '.txt'))))

anc_HRS_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

anc_HRS_afr[,ID:=paste0(SUBJID,"_",SUBJID)][,SUBJID:=NULL]
as.character(anc_HRS_afr$ID)-> anc_HRS_afr$ID
setkey(anc_HRS_afr, ID)


PGS2_HRS_afr<-vector('list', 4)
names(PGS2_HRS_afr)<-c("q1","q2","q3","q4")

for (I in names(PGS2_HRS_afr)){
        data.table(ID=rownames(PRS2), PGS=PRS2[,get(I)])-> PGS2_HRS_afr[[I]]
        setkey(PGS2_HRS_afr[[I]], ID)
        PGS2_HRS_afr[[I]][Pheno_HRS_afr, nomatch=0]-> PGS2_HRS_afr[[I]]
        PGS2_HRS_afr[[I]][anc_HRS_afr, nomatch=0]-> PGS2_HRS_afr[[I]]
        PGS2_HRS_afr[[I]][,AGE2:=AGE^2]
	PGS2_HRS_afr[[I]][,HEIGHT:=HEIGHT*100]
	PGS2_HRS_afr[[I]][which(!is.na(PGS2_HRS_afr[[I]][,HEIGHT])),]-> PGS2_HRS_afr[[I]]
}

#run linear models
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX, X))-> lm0_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS, X))-> lm1_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~AGE, X))-> lm2_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~AGE2, X))-> lm3_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~EUR_ANC, X))-> lm4_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS+AGE, X))-> lm5_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~PGS+AGE2, X))-> lm6_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2, X))-> lm7_HRS_afr
lapply(PGS2_HRS_afr, function(X) lm(HEIGHT~SEX+AGE+AGE2+PGS, X))-> lm8_HRS_afr

#Get partial R2, i.e, the proportion of variation in height explained by the PRS
partial_r2_HRS_afr<-lapply(1:length(PGS2_HRS_afr), function(X) partial.R2(lm7_HRS_afr[[X]], lm8_HRS_afr[[X]])) 
names(partial_r2_HRS_afr)<- names(PGS2_HRS_afr)


saveRDS(partial_r2_HRS_afr,file=paste0('~/height_prediction/strat_prs/output/part_R2_HRS_afr_',args[3], "_", rate.dist,  "_", w_map, '_v2.Rds')) #store results

