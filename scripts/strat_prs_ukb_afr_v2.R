#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
#Load packages #########################
suppressPackageStartupMessages({
library("optparse")
library(data.table)
library(dplyr)
library(asbio)
library(boot)
})
options(scipen=999)
########################################
#calculate PRS function
home<-"/home/bbita"
dir<-"height_prediction/strat_prs/scripts"
source(paste0('scripts/PRS_calc.R'))
#partial R2 function
source(paste0('scripts/Rsq_R2.R'))
########################################
#####################################################################################################
#### First part: calculate recombination rates for PRS SNPs and divide them into quantiles ##########
#####################################################################################################
## Read in betas and recombination maps

#TEST
#args<-c(6000, 'AA', 'gwas')

rate.dist<-as.numeric(args[1])
w_map<-args[2]
PRS<-vector('list', 22)
afr<-do.call(rbind, readRDS(paste0('~/height_prediction/', args[3], '/ukb_afr/output/hei_phys_100000_0.0005_v2.Rds')))[,.(CHR,POS,MarkerName,MarkerName,REF,ALT,Allele1,Allele2,b,p)]
if(args[2]=='AA'){
lapply(1:22, function(chr) fread(paste0('zcat /project/mathilab/data/maps/hm2/hm2/genetic_map_GRCh37_chr', chr,'.txt.gz'))[,CHR:=gsub("chr","",Chromosome)][, Chromosome:=NULL])-> rec #need to fix this path
for(chr in 1:22){colnames(rec[[chr]])<-c('POS', 'RATE_cM_Mb', 'MAP_cM', 'CHR')}
lapply(1:22, function(chr) fread(paste0('zcat /project/mathilab/data/maps/maps_b37/maps_chr.', chr, '.gz')))-> maps #need to fix this path
for(chr in 1:22){colnames(maps[[chr]])[1]<-"POS"}
lapply(1:22, function(chr) setkey(rec[[chr]],POS))
lapply(1:22, function(chr) setkey(maps[[chr]],POS))
lapply(1:22, function(chr) maps[[chr]][rec[[chr]], nomatch=0])-> map
remove(maps)
} else if(args[2]=='CEU'){
lapply(1:22, function(chr) fread(paste0('/project/mathilab/data/maps/1000_genomes_maps/hg19/CEU/CEU_recombination_map_hapmap_format_hg19_chr_', chr, '.txt')))-> map
for(chr in 1:22){colnames(map[[chr]])[2]<-"POS"}
for(chr in 1:22){map[[chr]][,CHR:=gsub("chr", "", Chromosome)][,Chromosome:=NULL]}
for(chr in 1:22){colnames(map[[chr]])<-c('POS', 'RATE_cM_Mb', 'MAP_cM', 'CHR')}
}
lapply(1:22, function(chr) map[[chr]][,POS2:=POS])
lapply(1:22, function(chr) arrange(map[[chr]], CHR, POS))
lapply(1:22, function(chr) setDT(map[[chr]]))
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
                f[, diff:=AA.rate(f$POS2)-AA.rate(f$POS)][,CHR:=chr]
                } else if (w_map=="CEU"){
                CEU.rate<-approxfun(map[[chr]]$POS, map[[chr]]$MAP_cM, rule=2)
                e<-data.table(POS=a, POS2=b, Win=paste0(a,"|",b))
                setkey(e, POS, POS2)
                setkey(map[[chr]], POS, POS2)
                f<-foverlaps(map[[chr]], e)
                f[, diff:=CEU.rate(f$POS2)-CEU.rate(f$POS)]
		}
                chrs[[chr]]<-f
                
}
do.call(rbind, chrs)-> f #combine into one data.table with all chromosomes
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
hei<-lapply(1:22, function(chr) readRDS(paste0('~/height_prediction/', args[3],'/ukb_afr/output/hei_phys_100000_0.0005_v2.Rds'))[[chr]]) 
hei<-do.call(rbind, hei)
hei[MarkerName %in% q1$MarkerName]-> hei_q1
hei[MarkerName %in% q2$MarkerName]-> hei_q2
hei[MarkerName %in% q3$MarkerName]-> hei_q3
hei[MarkerName %in% q4$MarkerName]-> hei_q4
prs<-vector('list', 4)
names(prs)<-c("q1","q2","q3", "q4")
prs[['q1']]<-PolScore(hei2=hei_q1)
cat('q1 done\n')
prs[['q2']]<-PolScore(hei2=hei_q2)
cat('q2 done\n')
prs[['q3']]<-PolScore(hei2=hei_q3)
cat('q3 done\n')
prs[['q4']]<-PolScore(hei2=hei_q4)
cat('q4 done\n')
PRS<-prs
remove(prs)
saveRDS(PRS,file=paste0("~/height_prediction/strat_prs/output/prs_ukb_afr_", args[3], '_', rate.dist, "_", w_map, "_v2.Rds")) #store results
obj<-c(nrow(hei_q1), nrow(hei_q2), nrow(hei_q3), nrow(hei_q4))
saveRDS(obj, file=paste0("~/height_prediction/strat_prs/output/Nr_SNPs_ukb_afr_", args[3], '_',rate.dist, "_", w_map, "_v2.Rds")) #store results
saveRDS(afr, file=paste0("~/height_prediction/strat_prs/output/rec_quant_ukb_afr_",args[3], '_',rate.dist, "_", w_map, "_v2.Rds")) #store results
#Make a list
PRS2<-vector('list', length(PRS[['q1']]))
names(PRS2)<-names(PRS[['q1']])
for(J in names(PRS2)){
        PRS2[[J]]<-data.table(q1=PRS[['q1']][[J]],q2=PRS[['q2']][[J]],q3=PRS[['q3']][[J]], q4=PRS[['q4']][[J]])
        cat(J, '\r')
}
do.call(rbind,PRS2)-> PRS2 #combine into one data.table
rownames(PRS2)<-names(PRS[[1]])
saveRDS(PRS2, file=paste0('~/height_prediction/strat_prs/output/PRS2_ukb_afr_', args[3], '_', rate.dist,"_", w_map, '_v2.Rds')) #store results


#########################################################
#### Second part: calculate partial R2 for PRS. #########
#########################################################

#read in phenotype data
fread('~/height_prediction/input/ukb_afr/UKB_AFR_pheno.txt')-> Pheno_UKB_afr #need to fix this path

#fix some columns and set keys for merging data tables
Pheno_UKB_afr[,ID:=paste0(ID, "_", ID)]
setkey(Pheno_UKB_afr, ID)

#add ancestry

ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/ukb_afr/rfmix_anc_chr', X, '.txt'))))

anc_ukb_afr<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

anc_ukb_afr[,ID:=SUBJID][,SUBJID:=NULL]
as.character(anc_ukb_afr$ID)-> anc_ukb_afr$ID
anc_ukb_afr[,ID:=paste0(ID, "_", ID)]
setkey(anc_ukb_afr, ID)


#Make a list to store results
PGS2_UKB_afr<-vector('list', 4)
names(PGS2_UKB_afr)<-c("q1","q2","q3","q4")
for (I in names(PGS2_UKB_afr)){
        data.table(ID=rownames(PRS2), PGS=PRS2[,get(I)])-> PGS2_UKB_afr[[I]]
        setkey(PGS2_UKB_afr[[I]], ID)
        PGS2_UKB_afr[[I]][Pheno_UKB_afr, nomatch=0]-> PGS2_UKB_afr[[I]]
        PGS2_UKB_afr[[I]][anc_ukb_afr, nomatch=0]-> PGS2_UKB_afr[[I]]
        PGS2_UKB_afr[[I]][,age2:=Age^2]
        PGS2_UKB_afr[[I]][AFR_ANC>=0.05]-> PGS2_UKB_afr[[I]]
        PGS2_UKB_afr[[I]][-which(is.na(PGS2_UKB_afr[[I]][,Height])),]-> PGS2_UKB_afr[[I]]
}

#run linear models
lapply(PGS2_UKB_afr, function(X) lm(Height~Sex, X))-> lm0_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS, X))-> lm1_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Age, X))-> lm2_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~age2, X))-> lm3_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~EUR_ANC, X))-> lm4_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS+Age, X))-> lm5_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~PGS+age2, X))-> lm6_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Sex+Age+age2+EUR_ANC, X))-> lm7_UKB_afr
lapply(PGS2_UKB_afr, function(X) lm(Height~Sex+Age+age2+EUR_ANC+PGS, X))-> lm8_UKB_afr

#Get partial R2, i.e, the proportion of variation in height explained by the PRS
lapply(c("q1","q2","q3","q4"), function(I) boot(data=PGS2_UKB_afr[[I]], statistic=rsq.R2,R=3000, formula1=Height~Sex+Age+age2+EUR_ANC, formula2=Height~Sex+Age+age2+EUR_ANC+PGS))-> results.UKB_afr
boots.ci.UKB_afr<-lapply(results.UKB_afr, function(X) boot.ci(X, type = c("norm", 'basic', "perc")))

partial_r2_UKB_afr<-lapply(1:length(PGS2_UKB_afr), function(X) partial.R2(lm7_UKB_afr[[X]], lm8_UKB_afr[[X]])) 
names(partial_r2_UKB_afr)<- names(PGS2_UKB_afr)


saveRDS(partial_r2_UKB_afr,file=paste0('strat_prs/output/part_R2_ukb_afr_', args[3], '_', rate.dist, "_", w_map, '_v2.Rds')) #store results
saveRDS(boots.ci.UKB_afr, file=paste0('strat_prs/output/results_ukb_afr_', args[3], '_', rate.dist,"_", w_map, '_v2.Rds'))


