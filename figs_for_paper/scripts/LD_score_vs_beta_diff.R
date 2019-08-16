setwd('~/height_prediction/figs_for_paper/')
system('wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2')
system('tar xvfj eur_w_ld_chr.tar.bz2')
system('rm eur_w_ld_chr.tar.bz2')
system('rm eur_w_ld_chr.tar.bz2.1')

library(data.table)
library(ggplot2)

ld<-do.call(rbind, lapply(1:22, function(X) fread(paste0('zcat ~/height_prediction/figs_for_paper/eur_w_ld_chr/', X,'.l2.ldscore.gz'))))
beta<-readRDS('~/height_prediction/gwas/ukb_afr/output/betas_phys_100000_0.0005_10000.Rds')
colnames(ld)[3]<-'POS'
setkey(ld, CHR, POS)
beta$CHR<-as.integer(beta$CHR)
setkey(beta, CHR, POS)
beta[ld,nomatch=0]-> test

test[,Quantile:=cut(L2, breaks=quantile(L2, probs=seq(0,1, by=0.1), na.rm=T), include.lowest=T)]
test[, MeanBetaDiffChisq:=mean(Beta_Diff_Chisq, na.rm=T), by=Quantile]
test[, MedianL2:=median(L2, na.rm=T), by=Quantile]


ggplot(test, aes(x=L2, y=MeanBetaDiffChisq)) + geom_point() + geom_smooth(method='lm', se=F)
ggsave('~/height_prediction/figs_for_paper/figs/test_LD.pdf')



ggplot(test, aes(x=MedianL2, y=MeanBetaDiffChisq)) + geom_point() + geom_smooth(method='lm', se=F) + geom_line()
ggsave('~/height_prediction/figs_for_paper/figs/test_LD_v2.pdf')
