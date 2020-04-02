#!/usr/bin/env Rscript
library(data.table)
library(snpStats)
library(ggplot2)
library(dplyr)
library(parallel)
library(mgcViz)
library(SOAR)
library(grid)
oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
Sys.setenv(R_LOCAL_CACHE=".R_Test")
#########################
cat('checkpoint number 1\n')
##
##Arslan's function
qexp=function(df){
  #number of rows in df
  nsnps=nrow(df)
  #generate P-values from chisquare statistic
  #sort df on P-value
  df=df%>%
    mutate(P=pchisq(Beta_Diff_Chisq,df=1,lower.tail = F))%>%
    arrange(P) %>% as.data.table
  #generate P-values from uniform distribution
  df=df%>%
    mutate(EXP_P=seq(1,nsnps)/(nsnps+1),
           EXP_CHI=qchisq(EXP_P,lower.tail=F,df=1)) %>% as.data.table
  return(df)
}
#############################

beta1<-readRDS('~/height_prediction/loc_anc_analysis/output/final_plink.Rds')
setkey(beta1, CHR, POS)
Store(beta1)
gc()
fread(paste0("zcat ~/height_prediction/gwas/input/50_raw_filtered.txt.gz"))-> ukb_height #read in GWAS summary statistics for height from the UK Biobank
gc()
ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele.
gc()
cat('checkpoint number 2\n')
ukb_height[CHR %in% 1:22]-> ukb_height
gc()
ukb_height[, N:=n_complete_samples][, AC:=NULL][, b:=beta][,p:=pval][,n_complete_samples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
gc()
ukb_height[,.(Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height
gc()
ukb_height$CHR<-as.integer(ukb_height$CHR)
gc()
ukb_height$POS<-as.integer(ukb_height$POS)
gc()
setkey(beta1, CHR, POS)
setkey(ukb_height, CHR, POS)
gc()
cat('checkpoint number 3\n')
beta2<-ukb_height[beta1, nomatch=0]
summary(beta2[PLINK<3 & PLINK>-3]) #this gives b and PLINK similar distributions.
Store(beta1, ukb_height)
gc()
beta2[,Beta_Diff:=b-PLINK]
beta2[,Beta_Diff_Chisq:=(Beta_Diff/sqrt(((SE^2)+(SE_plink^2))))^2]
gc()
pdf('~/height_prediction/figs_for_paper/figs/diff_qqplot_all.pdf')
#ggplot(beta2, aes(sample=Beta_Diff_Chisq)) + stat_qq(distribution=qchisq,dparams=list(df=1)) + stat_qq_line(distribution=qchisq,dparams=list(df=1))  
test_dat_all=qexp(beta2)
ggplot(test_dat_all)+
  geom_point(aes(EXP_CHI,Beta_Diff_Chisq))+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=4)) +
  labs(x=bquote(log[10]~"Expected"~chi^2),
       y=bquote(log[10]~"Observed"~chi^2))

dev.off()

gc()
cat('checkpoint number 4\n')
whi<-readRDS('~/height_prediction/gwas/WHI/output/plink_whi.Rds')
pdf('~/height_prediction/figs_for_paper/figs/diff_qqplot_PRS.pdf')
#b_plot<-ggplot(whi, aes(sample=Beta_Diff_Chisq)) + stat_qq(distribution='qchisq',dparams=list(df=1)) + stat_qq_line(distribution=qchisq,dparams=list(df=1))
test_dat=qexp(whi)
b_plot<-ggplot(test_dat)+
  geom_point(aes(EXP_CHI,Beta_Diff_Chisq))+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=3)) +
  labs(x=bquote(log[10]~"Expected"~chi^2),
       y=bquote(log[10]~"Observed"~chi^2))

print(b_plot)
dev.off()

gc()
#beta2[,PLINK:=scale(PLINK, center=-0.0005786,scale=7)] #mean and SD of b/sd of plink /10
gc()
cat('checkpoint number 5\n')
beta2$CHR<-as.integer(beta2$CHR)
beta2$POS<-as.integer(beta2$POS)
whi$CHR<-as.integer(whi$CHR)
whi$POS<-as.integer(whi$POS)
setkey(beta2, CHR, POS)
gc()
setkey(whi, CHR, POS)
gc()
cat('checkpoint number 6\n')
beta3<-beta2[whi, nomatch=0][,i.Beta_Diff:=NULL][,i.Beta_Diff_Chisq:=NULL]
summary(beta3[,.(PLINK, b)]) #here we can see that for PRS SNPs the distributions are pretty similar. 

png('~/height_prediction/figs_for_paper/figs/beta_cor_all.png')

ggplot(beta2, aes(x=b, y=PLINK)) + geom_point() + geom_smooth(method='lm')
dev.off()

Store(beta2)
gc()

pdf('~/height_prediction/figs_for_paper/figs/beta_cor_PRS.pdf')
a_plot<-ggplot(beta3, aes(x=b, y=PLINK)) + geom_point() + geom_smooth(method='lm') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18) 
, axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))+
labs(x=bquote("UKB_eur"~beta), y=bquote("UKB_afr"~beta))
print(a_plot)
dev.off()


png('~/height_prediction/figs_for_paper/figs/panel_inset_fig3.png', res=300, width=9, height=8, units="in")
c_plot<-a_plot+  annotation_custom(
    ggplotGrob(b_plot), 
    xmin = 0.25, xmax = 1, ymin = 1, ymax = 2
  )
print(c_plot)
#vp<-viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.2)
#print(b_plot, vp=vp)
dev.off()
Store(c_plot)
saveRDS(c_plot, "~/height_prediction/output/c_plot.Rds")
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
