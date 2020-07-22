#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
#Load packages #########################
suppressPackageStartupMessages({library("optparse")
library(data.table)
library(dplyr)
library(asbio)
library(boot)
})
options(scipen=999)
########################################
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
rate.dist<-as.numeric(args[1])
w_map<-args[2]
PRS<-vector('list', 22)
whi<-do.call(rbind, readRDS(paste0('~/height_prediction/', args[3], '/WHI/output/hei_phys_100000_0.0005_v2.Rds')))[,.(CHR,POS,MarkerName,REF,ALT,Allele1,Allele2,b,p)]
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
		f[, diff:=AA.rate(f$POS2)-AA.rate(f$POS)][, CHR:=chr]
		chrs[[chr]]<-f
		} else if (w_map=="CEU"){
                CEU.rate<-approxfun(map[[chr]]$POS, map[[chr]]$MAP_cM, rule=2)
                e<-data.table(POS=a, POS2=b, Win=paste0(a,"|",b))
                setkey(e, POS, POS2)
                setkey(map[[chr]], POS, POS2)
                f<-foverlaps(map[[chr]], e)
                f[, diff:=CEU.rate(f$POS2)-CEU.rate(f$POS)][, CHR:=chr]
		chrs[[chr]]<-f
		}
}
do.call(rbind, chrs)-> f #combine into one data.table with all chromosomes
f[,Quantile:=cut(diff, breaks=quantile(diff), na.rm=T, include.lowest=T)]
whi[, POS1:=POS-(rate.dist/2)][,POS2:=POS+(rate.dist/2)]

#Stratify genome into 4 quantiles of recombination rate
if(w_map=="CEU"){
	whi[, diff:=CEU.rate(POS2)-CEU.rate(POS1), by=CHR]
        } else if (w_map=="AA"){
	 whi[, diff:=AA.rate(POS2)-AA.rate(POS1), by=CHR]
}
whi[,Quantile:=cut(diff, breaks=quantile(diff), na.rm=T, include.lowest=T)]
q1<-split(whi, by='Quantile')[[1]]
q2<-split(whi, by='Quantile')[[2]]
q3<-split(whi, by='Quantile')[[3]]
q4<-split(whi, by='Quantile')[[4]]
whi$Quantile
cat('checkpoint\n')
#calculate PRS for each of these quantiles
hei<-lapply(1:22, function(chr) readRDS(paste0('~/height_prediction/', args[3], '/WHI/output/hei_phys_100000_0.0005_v2.Rds'))[[chr]]) #need to fix this path
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
saveRDS(PRS,file=paste0("~/height_prediction/strat_prs/output/prs_WHI_", args[3],'_',  rate.dist, "_", w_map, "_v2.Rds")) #store results
obj<-c(nrow(hei_q1), nrow(hei_q2), nrow(hei_q3), nrow(hei_q4))
saveRDS(obj, file=paste0("~/height_prediction/strat_prs/output/Nr_SNPs_WHI_",args[3], '_',rate.dist, "_", w_map, "_v2.Rds")) #store results
saveRDS(whi, file=paste0("~/height_prediction/strat_prs/output/rec_quant_WHI_",args[3], '_',rate.dist, "_", w_map, "_v2.Rds")) #store results
#Make a list
PRS2<-vector('list', length(PRS[['q1']]))
names(PRS2)<-names(PRS[['q1']])
for(J in names(PRS2)){
	PRS2[[J]]<-data.table(q1=PRS[['q1']][[J]],q2=PRS[['q2']][[J]],q3=PRS[['q3']][[J]], q4=PRS[['q4']][[J]])
	cat(J, '\r')
}
do.call(rbind,PRS2)-> PRS2 #combine into one data.table
rownames(PRS2)<-names(PRS[[1]])
saveRDS(PRS2, file=paste0('~/height_prediction/strat_prs/output/PRS2_WHI_', args[3], '_', rate.dist,"_", w_map, '_v2.Rds')) #store results

#########################################################
#### Second part: calculate partial R2 for PRS. #########
#########################################################

#read in phenotype data
fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI #need to fix this path

#fix some columns and set keys for merging data tables
Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)


#add ancestry
#add ancestry
ancestry<-do.call(rbind, lapply(1:22, function(X) fread(paste0('~/height_prediction/input/WHI/rfmix_anc_chr', X, '.txt'))))

anc_WHI<-ancestry %>% group_by(SUBJID) %>% summarise(AFR_ANC=mean(AFR_ANC), EUR_ANC=1-mean(AFR_ANC)) %>% as.data.table #mean across chromosomes for each individual

setkey(anc_WHI, SUBJID)

cat('All good so far\n')

Pheno_WHI[anc_WHI, nomatch=0]-> PGS_WHI
PGS_WHI[,AGE2:=AGE^2]
sd1_f<-sd(PGS_WHI$HEIGHTX, na.rm=T)
m1_f<-mean(PGS_WHI$HEIGHTX, na.rm=T)
PGS_WHI<-PGS_WHI[HEIGHTX>=m1_f-(2*sd1_f) & HEIGHTX<=m1_f+(2*sd1_f)]
nrow(PGS_WHI)
setkey(PGS_WHI, SUBJID)
nrow(PGS_WHI)

PGS2_WHI<-vector('list', 4)
names(PGS2_WHI)<-c("q1","q2","q3","q4")

for (I in names(PGS2_WHI)){
        data.table(ID=rownames(PRS2), PGS=PRS2[,get(I)])-> PGS2_WHI[[I]]
        setkey(PGS2_WHI[[I]], ID)
        PGS2_WHI[[I]]<-PGS2_WHI[[I]][PGS_WHI, nomatch=0]
}

#run linear models
lapply(PGS2_WHI, function(X) lm(HEIGHTX~PGS, X))-> lm1_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~AGE, X))-> lm2_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~AGE2, X))-> lm3_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~EUR_ANC, X))-> lm4_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~PGS+AGE, X))-> lm5_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~PGS+AGE2, X))-> lm6_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~AGE+AGE2+EUR_ANC, X))-> lm7_WHI
lapply(PGS2_WHI, function(X) lm(HEIGHTX~AGE+AGE2+EUR_ANC+PGS, X))-> lm8_WHI

lapply(c("q1","q2","q3","q4"), function(I) boot(data=PGS2_WHI[[I]], statistic=rsq.R2,R=3000, formula1=HEIGHTX~AGE+AGE2+EUR_ANC, formula2=HEIGHTX~AGE+AGE2+EUR_ANC+PGS))-> results.WHI
boots.ci.WHI<-lapply(results.WHI, function(X) boot.ci(X, type = c("norm", 'basic', "perc")))
#Get partial R2, i.e, the proportion of variation in height explained by the PRS
partial_r2_WHI<-lapply(1:length(PGS2_WHI), function(X) partial.R2(lm7_WHI[[X]], lm8_WHI[[X]]))
names(partial_r2_WHI)<- names(PGS2_WHI)

saveRDS(partial_r2_WHI,file=paste0('~/height_prediction/strat_prs/output/part_R2_WHI_', args[3], '_', rate.dist,"_", w_map, '_v2.Rds'))  #store results.
saveRDS(boots.ci.WHI, file=paste0('~/height_prediction/strat_prs/output/results_WHI_', args[3], '_', rate.dist,"_", w_map, '_v2.Rds'))
