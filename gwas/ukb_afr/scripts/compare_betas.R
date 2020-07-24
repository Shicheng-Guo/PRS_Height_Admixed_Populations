library(data.table)
library(dplyr)
source('~/height_prediction/scripts/my_manhattan.R')
plink<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/association_v3.Res.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'
gc()
#
plink2<-fread('~/height_prediction/runSmartpCA-master/UKB_AFR/test3.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","A1","TEST"," OBS_CT","PLINK", "SE","T_STAT", "UNADJ")
setkey(plink, MarkerName, CHR, UNADJ)
setkey(plink2, MarkerName, CHR, UNADJ)
plink[plink2, nomatch=0]-> final_plink
final_plink$CHR<-as.numeric(final_plink$CHR)
setkey(final_plink, CHR, POS)
final_plink[order(CHR,POS)]-> final_plink
remove(plink, plink2)
gc()

#read file with UKB effect sizes and get the positions of SNPs
ukb_height<-fread('zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz')[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]
gc()
ukb_height[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
gc()
ukb_height[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
gc()
ukb_height[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height
gc()
ukb_height$CHR<-as.numeric(ukb_height$CHR)
gc()
ukb_height$POS<-as.numeric(ukb_height$POS)
gc()
setkey(ukb_height,CHR, POS)
gc()
#combine
final_plink[ukb_height, nomatch=0]-> combo
remove(ukb_height, final_plink)
gc()
combo[,TEST:=NULL]
combo[,Beta_Diff:=abs(PLINK-b)]
combo[, Beta_Diff_Chisq:=((PLINK-b)^2)/((SE)^2+(i.SE)^2)]
fwrite(combo, file="~/height_prediction/gwas/ukb_afr/output/plink_effect_sizes.txt")
#read local ancestry stuff

local_anc<-vector('list', 22)
for(I in 1:22){
	local_anc[[I]]<-fread(paste0('~/height_prediction/gwas/ukb_afr/output/AS_Beta_chr', I, 'example.txt'))
}

do.call(rbind, local_anc)-> local_anc

colnames(local_anc)[4]<-"MarkerName"
setkey(local_anc, CHR, POS)
setkey(combo, CHR, POS)

combo[local_anc, nomatch=0]-> combo_local
fwrite(combo_local, file="~/height_prediction/gwas/ukb_afr/output/plink_local_anc_effect_sizes.txt")
#now restrict to onlu PRS snps
prs_snp<-do.call(rbind, readRDS('~/height_prediction/gwas/ukb_afr/output/hei_phys_100000_0.0005.Rds'))[,1:5]
setkey(combo_local, CHR,POS)
colnames(prs_snp)[1]<-'CHR'
setkey(prs_snp, CHR, POS)
combo_prs<-combo_local[prs_snp, nomatch=0]

saveRDS(combo_prs, file='~/height_prediction/gwas/ukb_afr/output/combo_prs.Rds')
fwrite(combo_prs, file='~/height_prediction/gwas/ukb_afr/output/combo_prs.txt')

with(combo_prs, cor.test(PLINK,ALL)) #77.8%

with(combo_local, cor.test(PLINK,ALL)) #70.5%

###


##QQ plot for plink UKB

observed <- sort(combo_local$UNADJ)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

observed2 <- sort(combo_local$GC)
lobs2 <- -(log10(observed2))

expected2 <- c(1:length(observed2))
lexp2 <- -(log10(expected2 / (length(expected2)+1)))

df<-data.table(Obs=c(lobs, lobs2), Exp=c(lexp, lexp2), Class=c(rep("UNADJ", length(observed)),rep("GC", length(observed2))))

ggplot(df, aes(x=Exp, y=Obs, colour=Class)) + geom_point() +   geom_abline(intercept = 0, slope = 1, col="black") + ggtitle("Plink GWAS for UKB admixed inviduals")
ggsave('~/height_prediction/runSmartpCA-master/UKB_AFR/figs/qqplot_plink_UKB.pdf')

