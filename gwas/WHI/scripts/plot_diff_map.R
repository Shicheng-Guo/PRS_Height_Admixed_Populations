#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)

}
#Load packages
library("optparse")
library(data.table)
library(biomaRt)
library(ggplot2)

beta<-vector('list', 22)
prun<-args[1]
#prun<-"phys_100000_0.0005"
hei<-readRDS(paste0('~/height_prediction/gwas/WHI/output/hei_', prun, '_v2.Rds'))
for(chr in 1:22){
#chr<-22
cat('start chr ')
cat(chr)
cat('\n')
#rate.dist <- 20000
rate.dist<-as.numeric(args[2])
betas <-read.table(paste0('~/height_prediction/gwas/WHI/output/AS_Beta_chr', chr, 'example.txt'), as.is=T, header=T)
maps<-fread(paste0('zcat /project/mathilab/data/maps_b37/maps_chr.', chr, '.gz'))
snps <- read.table(paste0("~/height_prediction/input/WHI/WHI_b37_strand_include_kgCY_chr", chr, ".phsnp"), as.is=TRUE)
colnames(snps) <- c("ID", "CHR", "Map", "POS", "REF", "ALT")
betas <- cbind(snps, betas)
setDT(betas)
snp_list<-hei[[chr]]
#snp_list<-readRDS(paste0('/project/mathilab/bbita/gwas_admix/new_height/WHI/prunned_1kg/LD_prunned_hei_chr', chr, "_", prun, ".Rds"))
grch37_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_snp")
ID_posgrch37<-getBM(attributes = c('refsnp_id','allele', 'chr_name','chrom_start'),filters = 'snp_filter',values = snp_list$i.MarkerName,mart = grch37_snp)
setDT(ID_posgrch37)
colnames(ID_posgrch37)<-c('MarkerName','Allele','CHR','POS')
ID_posgrch37[CHR==chr]-> ID_posgrch37

cat('checkpoint \n')
AA.rate <- approxfun(maps$Physical_Pos, maps$AA_Map)
YRI.rate <- approxfun(maps$Physical_Pos, maps$YRI_LD)
CEU.rate <- approxfun(maps$Physical_Pos, maps$CEU_LD)
COMBINED.rate<-approxfun(maps$Physical_Pos,maps$COMBINED_LD)
#TODO Have to think about what to do about point that are outside the map...
#Now estimate the AA rate and difference in rates at each SNP
betas$AA.rate <- 0
betas$CEU_YRI_diff.rate <- 0
betas$CEU.rate <- 0
betas$YRI.rate <- 0
betas$COMBINED.rate<-0
for(i in 1:NROW(betas)){
    cat(i, '\r')
    AA.x <- AA.rate(betas$POS[i]+(rate.dist/2))-AA.rate(betas$POS[i]-(rate.dist/2))
    betas$AA.rate[i] <- AA.x
    CEU.x <- CEU.rate(betas$POS[i]+(rate.dist/2))-CEU.rate(betas$POS[i]-(rate.dist/2))
    betas$CEU.rate[i] <- CEU.x
    YRI.x <- YRI.rate(betas$POS[i]+(rate.dist/2))-YRI.rate(betas$POS[i]-(rate.dist/2))
    betas$YRI.rate[i]<-YRI.x
    COMBINED.x<-COMBINED.rate(betas$POS[i]+(rate.dist/2))-COMBINED.rate(betas$POS[i]-(rate.dist/2))
    betas$COMBINED.rate[i]<-COMBINED.x
    betas$CEU_YRI_diff.rate[i] <- abs(CEU.x-YRI.x)
}
betas$y <- abs(betas$POP1-betas$POP2)
#At this point you will restrict the betas to the SNPS that you are using in the PRS
betas<-betas[which(betas$POS %in% ID_posgrch37$POS),] #
setkey(betas, CHR, POS)
setkey(snp_list, CHR, POS)
betas[snp_list[,.(CHR, POS, b,SE)]]-> betas
betas[,Tstat_ukb:=b/SE]
betas[,y2:=abs(betas$b-betas$POP1)]
betas[,Beta_Diff_Chisq:=((Tstat_pop1-Tstat_ukb)^2)/2]
beta[[chr]]<-betas
cat('chr ')
cat(chr)
cat(' done\n')
}
#saveRDS(beta, 'betas.Rds')
#Plot against AA map - no significant effect!
#pdf(paste0("chr", chr, "_rate_against_AA_map.pdf"))
do.call(rbind, beta)-> beta
setDT(beta)
beta[order(CHR, POS)]-> beta
as.factor(beta$CHR)-> beta$CHR
pdf(paste0('~/height_prediction/gwas/WHI/figs/all_chr_rate_against_AA_map_', args[2], '.pdf'))
ggplot(beta, aes(x=AA.rate, y=Beta_Diff_Chisq)) + geom_point() + geom_smooth(method='lm') + facet_wrap(~CHR, nrow=5, scales='free') 
dev.off()

#Plot against AA map - no significant effect! If anything a negative effect. 
pdf(paste0("~/height_prediction/gwas/WHI/figs/all_chr_rate_against_Diff_map_", args[2], ".pdf"))
ggplot(beta, aes(x=CEU_YRI_diff.rate,y=Beta_Diff_Chisq)) + geom_point() + geom_smooth(method='lm') + ggtitle(paste0(args[2], " distance WHI"))
dev.off()

pdf(paste0("~/height_prediction/gwas/WHI/figs/all_chr_rate_against_AA_map_combined", args[2], ".pdf"))
ggplot(beta[AA.rate<=0.05], aes(x=AA.rate, y=Beta_Diff_Chisq)) + geom_point() + geom_smooth(method='lm') + ggtitle(paste0(args[2], " distance WHI"))
dev.off()

pdf(paste0("~/height_prediction/gwas/WHI/figs/all_chr_rate_against_CEU_map_combined", args[2], ".pdf"))
ggplot(beta[CEU.rate<=0.05], aes(x=CEU.rate, y=Beta_Diff_Chisq)) + geom_point() + geom_smooth(method='lm') + ggtitle(paste0(args[2], " distance WHI"))
dev.off()

saveRDS(beta, file=paste0('~/height_prediction/gwas/WHI/output/betas_', args[1], '_', args[2], '.Rds'))


