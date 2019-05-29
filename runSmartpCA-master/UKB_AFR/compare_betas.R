library(data.table)
library(dplyr)
source('~/height_prediction/scripts/my_manhattan.R')
plink<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/association_v3.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'
#
plink2<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/test.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","A1","TEST"," OBS_CT","BETA", "SE","T_STAT", "UNADJ")
setkey(plink, MarkerName, CHR, UNADJ)
setkey(plink2, MarkerName, CHR, UNADJ)
plink[plink2, nomatch=0]-> final_plink
final_plink$CHR<-as.numeric(final_plink$CHR)
arrange(final_plink, CHR,POS) %>% as.data.table -> final_plink
setkey(final_plink, MarkerName, CHR, POS)
#read file with UKB effect sizes and get the positions of SNPs
ukb_height<-fread('zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz')
ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL] 
ukb_height[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
ukb_height[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
ukb_height[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height
ukb_height$CHR<-as.numeric(ukb_height$CHR)
ukb_height$POS<-as.numeric(ukb_height$POS)
setkey(ukb_height, MarkerName, CHR, POS)

#combine

final_plink[ukb_height, nomatch=0]-> combo
combo[,TEST:=NULL]
#read local ancestry stuff

local_anc<-vector('list', 22)
for(I in 1:22){
	local_anc[[I]]<-fread(paste0('~/height_prediction/gwas/ukb_afr/output/AS_Beta_chr', I, 'example.txt'))
}

do.call(rbind, local_anc)-> local_anc

colnames(local_anc)[4]<-"MarkerName"
setkey(local_anc, MarkerName, CHR, POS)

combo[local_anc, nomatch=0]-> combo_local

#now restrict to onlu PRS snps
fread('/project/mathilab/bbita/gwas_admix/new_height/ukb_afr_betas_100000_0.0005.txt')-> prs_snps

combo_prs<-combo_local[MarkerName %in% prs_snps$MarkerName]

with(combo_prs, cor.test(BETA,POP1)) #POP1 is AFR in the local ancestry analysis, POP2 is EUR, ALL is estimate from all chromosomes. T



##QQ plot for plink UKB

observed <- sort(combo_prs$UNADJ)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))


observed2 <- sort(combo_prs$GC)
lobs2 <- -(log10(observed2))

expected2 <- c(1:length(observed2))
lexp2 <- -(log10(expected2 / (length(expected2)+1)))


df<-data.table(Obs=c(lobs, lobs2), Exp=c(lexp, lexp2), Class=c(rep("UNADJ", length(observed)),rep("GC", length(observed2))))

ggplot(df, aes(x=Exp, y=Obs, colour=Class)) + geom_point() +   geom_abline(intercept = 0, slope = 1, col="black") + ggtitle("Plink GWAS for UKB admixed inviduals")
ggsave('figs/qqplot_LA_UKB.pdf')
#this works



##manhattan for plink

library(qqman)
df2<-combo_prs[,.(CHR,POS,MarkerName,GC)]
df2[,BP:=POS][,POS:=NULL]
df2[,P:=GC][,GC:=NULL]
df2[,SNP:=MarkerName][,MarkerName:=NULL]
arrange(df2, CHR, BP) %>% as.data.table-> df2
