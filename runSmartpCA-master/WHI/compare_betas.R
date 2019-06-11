library(data.table)
library(ggplot2)
library(dplyr)
source('~/height_prediction/scripts/my_manhattan.R')
plink<-fread('~/height_prediction/runSmartpCA-master/WHI/association_v3.Res.Height.glm.linear.adjusted', header=T, fill=T)
colnames(plink)[2]<-'MarkerName'
colnames(plink)[1]<-'CHR'
#
plink2<-fread('~/height_prediction/runSmartpCA-master/WHI/test3.txt', fill=T)
colnames(plink2)<-c("CHR","POS", "MarkerName","REF","ALT","A1","TEST"," OBS_CT","BETA", "SE","T_STAT", "UNADJ")
setkey(plink, MarkerName, CHR, UNADJ)
setkey(plink2, MarkerName, CHR, UNADJ)
plink[plink2, nomatch=0]-> final_plink
final_plink$CHR<-as.numeric(final_plink$CHR)
arrange(final_plink, CHR,POS) %>% as.data.table -> final_plink
setkey(final_plink,CHR, POS)
#read file with UKB effect sizes and get the positions of SNPs
ukb_height<-fread('zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz')
ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL] 
ukb_height[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
ukb_height[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
ukb_height[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height
ukb_height$CHR<-as.numeric(ukb_height$CHR)
ukb_height$POS<-as.numeric(ukb_height$POS)
setkey(ukb_height, CHR, POS)

#combine
gc()

final_plink[ukb_height, nomatch=0]-> combo
combo[,TEST:=NULL]
#read local ancestry stuff

local_anc<-vector('list', 22)
for(I in 1:22){
	local_anc[[I]]<-fread(paste0('~/height_prediction/gwas/WHI/output/AS_Beta_chr', I, 'example.txt'))
}

do.call(rbind, local_anc)-> local_anc

colnames(local_anc)[4]<-"MarkerName"
setkey(local_anc, MarkerName, CHR, POS)
setkey(combo, MarkerName, CHR, POS)

combo[local_anc, nomatch=0]-> combo_local

#now restrict to onlu PRS snps
fread('/project/mathilab/bbita/gwas_admix/new_height/whi_betas_100000_0.0005.txt')-> prs_snps

combo_prs<-combo_local[MarkerName %in% prs_snps$MarkerName]

with(combo_prs, cor.test(BETA,POP1)) # 0.70
with(combo_prs, cor.test(BETA,POP2)) #0.58
with(combo_prs, cor.test(BETA,ALL)) #0.977


with(combo_local, cor.test(BETA,POP1)) # 0.698922
with(combo_local, cor.test(BETA,POP2)) #0.517691
with(combo_local, cor.test(BETA,ALL)) #0.9515042

combo_prs[, Beta_Diff:=abs(BETA-ALL)] #investigate which snps have large beta diff
combo_local[,Beta_Diff:=abs(BETA-ALL)]

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

ggplot(df, aes(x=Exp, y=Obs, colour=Class)) + geom_point() +   geom_abline(intercept = 0, slope = 1, col="black") + ggtitle("Plink GWAS for WHI admixed inviduals")
ggsave('~/height_prediction/runSmartpCA-master/WHI/figs/qqplot_plink_WHI.pdf')
dev.off()

library(qqman)
combo_local[, P_all:=2*(pt(abs(Tstat_all), df=N-13, lower.tail=F))]
combo_local[, P_pop1:=2*(pt(abs(Tstat_pop1), df=N-13, lower.tail=F))]
combo_local[, P_pop2:=2*(pt(abs(Tstat_pop2), df=N-13, lower.tail=F))]
combo_local[, P_all_GC:= pchisq(qchisq(P_all, df=1)/qchisq(0.5, df=1), df=1)]
combo_local[, P_pop1_GC:= pchisq(qchisq(P_pop1, df=1)/qchisq(0.5, df=1), df=1)]
combo_local[, P_pop2_GC:= pchisq(qchisq(P_pop2, df=1)/qchisq(0.5, df=1), df=1)]

df2<-combo_local[,.(CHR,POS,MarkerName,GC)]
df2[,BP:=POS][,POS:=NULL]
df2[,P:=GC][,GC:=NULL]
df2[,SNP:=MarkerName][,MarkerName:=NULL]
arrange(df2, CHR, BP) %>% as.data.table-> df2

pdf('~/height_prediction/runSmartpCA-master/WHI/figs/manhattan_plink_WHI.pdf')
manhattan(df2, genomewideline=F, suggestiveline=F)
dev.off()

df2<-combo_local[,.(CHR,POS,MarkerName,P_all)]
df2[,P:=P_all][,P_all:=NULL]
df2[,SNP:=MarkerName][,MarkerName:=NULL]
df2[,BP:=POS][,POS:=NULL]
arrange(df2, CHR, BP) %>% as.data.table-> df2

pdf('~/height_prediction/runSmartpCA-master/WHI/figs/manhattan_LA_WHI.pdf')
manhattan(df2, genomewideline=F, suggestiveline=F)
dev.off()
#QQ plot for LA WHI

observed <- sort(combo_local$P_all)
lobs <- -(log10(observed))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

observed_GC <- sort(combo_local$P_all_GC)
lobs_GC <- -(log10(observed_GC))

expected_GC <- c(1:length(observed_GC))
lexp_GC <- -(log10(expected_GC / (length(expected_GC)+1)))

observed2 <- sort(combo_local$P_pop1)
lobs2 <- -(log10(observed2))

expected2 <- c(1:length(observed2))
lexp2 <- -(log10(expected2 / (length(expected2)+1)))

observed2_GC <- sort(combo_local$P_pop1_GC)
lobs2_GC <- -(log10(observed2_GC))


expected2_GC <- c(1:length(observed2_GC))
lexp2_GC <- -(log10(expected2_GC / (length(expected2_GC)+1)))

observed3 <- sort(combo_local$P_pop2)
lobs3 <- -(log10(observed3))

expected3 <- c(1:length(observed3))
lexp3 <- -(log10(expected3 / (length(expected3)+1)))

observed3_GC <- sort(combo_local$P_pop2_GC)
lobs3_GC <- -(log10(observed3_GC))

expected3_GC <- c(1:length(observed3_GC))
lexp3_GC <- -(log10(expected3_GC / (length(expected3_GC)+1)))

df<-data.table(Obs=c(lobs, lobs_GC, lobs2,lobs2_GC, lobs3, lobs3_GC), Exp=c(lexp, lexp_GC,lexp2,lexp2_GC, lexp3, lexp3_GC), Class=c(rep("ALL", length(observed)), rep("ALL_GC", length(observed_GC)),rep("POP1", length(observed2)), rep("POP1_GC", length(observed2_GC)), rep('POP2',length(observed3)), rep("POP2_GC", length(observed3_GC))))

gc()

pdf('~/height_prediction/runSmartpCA-master/WHI/figs/qqplot_LA_WHI.pdf')
ggplot(df, aes(x=Exp, y=Obs, colour=Class)) + geom_point(cex=0.5) +   geom_abline(intercept = 0, slope = 1, col="black") + ggtitle("Local ancestry GWAS for WHI admixed inviduals")
dev.off()
dev.off()
