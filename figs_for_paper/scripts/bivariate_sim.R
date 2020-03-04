#!/usr/bin/env Rscript
library(mvtnorm)
library(ggplot2)
library(wesanderson)
library("optparse")
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(wesanderson)
library(rlist)
library(asbio)
library(GGally)
library(tidyr)
library(hexbin)
library(psychometric)
library(boot)
library(RColorBrewer)
library(cowplot)

####################################
args<-"gwas"
#OR plots
readRDS(paste0('~/height_prediction/', args[1], '/WHI/output/results.WHI.Rds'))-> results.WHI
readRDS(paste0('~/height_prediction/',args[1], '/JHS/output/results.JHS.Rds'))-> results.JHS
readRDS(paste0('~/height_prediction/',args[1], '/ukb_afr/output/results.UKB_afr.Rds'))-> results.UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/results.HRS_eur.Rds'))-> results.HRS_eur
readRDS(paste0('~/height_prediction/',args[1], '/HRS_afr/output/results.HRS_afr.Rds'))-> results.HRS_afr
readRDS(paste0('~/height_prediction/',args[1], '/WHI/output/B_WHI.Rds'))->B_WHI
readRDS(paste0('~/height_prediction/',args[1], '/JHS/output/B_JHS.Rds'))->B_JHS
readRDS(paste0('~/height_prediction/',args[1], '/ukb_afr/output/B_UKB_afr.Rds'))-> B_UKB_afr
readRDS(paste0('~/height_prediction/',args[1], '/HRS_eur/output/B_HRS_eur.Rds'))-> B_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/B_HRS_afr.Rds'))-> B_HRS_afr
readRDS(paste0('~/height_prediction/',args[1],'/ukb_afr/output/PGS3_UKB_afr.Rds'))-> PGS3_UKB_afr
readRDS(paste0('~/height_prediction/',args[1],'/WHI/output/PGS3_WHI.Rds'))-> PGS3_WHI
readRDS(paste0('~/height_prediction/',args[1],'/JHS/output/PGS3_JHS.Rds'))-> PGS3_JHS
readRDS(paste0('~/height_prediction/',args[1],'/HRS_eur/output/PGS3_HRS_eur.Rds'))-> PGS3_HRS_eur
readRDS(paste0('~/height_prediction/',args[1],'/HRS_afr/output/PGS3_HRS_afr.Rds'))-> PGS3_HRS_afr

for(I in names(B_JHS)){ #JHS lacks the LD prunning methods
        ALL<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
        ALL$Dataset<-factor(ALL$Dataset, levels=c("UKB_afr", "WHI_afr", "JHS_afr", "HRS_afr",  "HRS_eur"))
}


ALL2<-vector('list', length(names(B_JHS)))
names(ALL2)<-names(B_JHS)

for(I in names(B_JHS)){
                ALL2[[I]]<-rbind(B_JHS[[I]][1:2,][, Dataset:='JHS_afr'], B_WHI[[I]][1:4,][, Dataset:='WHI_afr'], B_UKB_afr[[I]][1:4,][,Dataset:='UKB_afr'],B_HRS_afr[[I]][1:2,][, Dataset:='HRS_afr'],  B_HRS_eur[[I]])
        #rbind(ALL2[[I]][!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL2[[I]][Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL2[[I]]
                tmp<-lm(R_sq~Med_Eur_Anc,weights=1/
        c(var(results.JHS[[I]][[1]]$t), var(results.JHS[[I]][[2]]$t), var(results.WHI[[I]][[1]]$t), var(results.WHI[[I]][[2]]$t), var(results.WHI[[I]][[3]]$t), var(results.WHI[[I]][[4]]$t), var(results.UKB_afr[[I]][[1]]$t),var(results.UKB_afr[[I]][[2]]$t), var(results.UKB_afr[[I]][[3]]$t), var(results.UKB_afr[[I]][[4]]$t), var(results.HRS_afr[[I]][[1]]$t), var(results.HRS_afr[[I]][[2]]$t),var(results.HRS_eur[[I]]$t)), data=ALL2[[I]])
#       }
        #rbind(ALL3[[I]][!(Dataset %in% c('UKB_EUR', 'pennBB_EA'))], ALL3[[I]][Dataset %in% c('UKB_EUR', 'pennBB_EA')][, Med_Eur_Anc:=1])-> ALL3[[I]]
        ALL2[[I]][,Set:=I]
        readRDS(paste0('~/height_prediction/', args[1],'/WHI/output/Nr_SNPs_WHI.Rds'))[Name==I][, Nr]->a
        readRDS(paste0('~/height_prediction/', args[1],'/ukb_afr/output/Nr_SNPs_UKB_afr.Rds'))[Name==I][, Nr]->b
        readRDS(paste0('~/height_prediction/', args[1],'/JHS/output/Nr_SNPs_JHS.Rds'))[Name==I][, Nr]->d
        readRDS(paste0('~/height_prediction/', args[1],'/HRS_eur/output/Nr_SNPs_HRS.Rds'))[Name==I][, Nr]->f
        #readRDS('Nr_SNPs_pennBB_afr.Rds')[Name==I][, Nr]->f
        #readRDS('Nr_SNPs_pennBB_eur.Rds')[Name==I][, Nr]->g
        ALL2[[I]][,Intercept:=coef(tmp)[[1]]][,Slope:=coef(tmp)[[2]]]
        ALL2[[I]][,Slope_Intercept:=sum(coef(tmp))]
        ALL2[[I]][, Nr_SNPs_WHI:=a]
        ALL2[[I]][, Nr_SNPs_UKB:=b]
        ALL2[[I]][, Nr_SNPs_JHS:=d]
        ALL2[[I]][, Nr_SNPs_HRS_eur:=f]
        ALL2[[I]][, Nr_SNPs_HRS_afr:=f]
        cat(I, ' \n')
}
do.call(rbind,ALL2)[,.(Set,Intercept,Slope_Intercept, Slope, Nr_SNPs_WHI, Nr_SNPs_UKB, R_sq, Med_Eur_Anc)]->ALL4

combo<-vector('list', length(PGS3_JHS))
names(combo)<-names(PGS3_JHS)

for (I in names(PGS3_JHS)){
        rbind(PGS3_WHI[[I]][,.(SUBJID,AGE, age2, HEIGHTX,PGS, SEX,EUR_ANC)][,SUBJ_ID:=SUBJID][, age:=AGE][, sex:='FEMALE'][, SUBJID:=NULL][, AGE:=NULL][, SEX:=NULL][,Dataset:='WHI_afr'][, Res.Height:=resid(lm(HEIGHTX~age+age2))],
        PGS3_HRS_afr[[I]][,SUBJ_ID:=ID][, age:=AGE][, age2:=AGE2][, HEIGHTX:=HEIGHT][, sex:=SEX][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='HRS_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
        PGS3_HRS_eur[[I]][,SUBJ_ID:=ID][, age:=AGE][, age2:=AGE2][, HEIGHTX:=HEIGHT][, sex:=SEX][,.(SUBJ_ID,age, age2, HEIGHTX,PGS, sex, EUR_ANC)][,Dataset:='HRS_eur'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
        PGS3_UKB_afr[[I]][, sex:=Sex][, HEIGHTX:=Height][,SUBJ_ID:=ID][,.(SUBJ_ID, Age, age2, HEIGHTX, PGS, sex, EUR_ANC)][, age:=Age][, Age:=NULL][,Dataset:='UKB_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))],
        PGS3_JHS[[I]][,.(SUBJID, age, age2, HEIGHTX, PGS, sex, EUR_ANC)][,SUBJ_ID:=SUBJID][, SUBJID:=NULL][,Dataset:='JHS_afr'][, Res.Height:=resid(lm(HEIGHTX~sex+age+age2))])[,Prun_Set:=I][,Nr_SNPs:=0]-> combo[[I]]
        readRDS(paste0('~/height_prediction/', args[1],'/WHI/output/Nr_SNPs_WHI.Rds'))[Name==I][, Nr]->a
        readRDS(paste0('~/height_prediction/', args[1],'/ukb_afr/output/Nr_SNPs_UKB_afr.Rds'))[Name==I][, Nr]->b
        readRDS(paste0('~/height_prediction/', args[1],'/JHS/output/Nr_SNPs_JHS.Rds'))[Name==I][, Nr]->d
        readRDS(paste0('~/height_prediction/', args[1],'/HRS_eur/output/Nr_SNPs_HRS.Rds'))[Name==I][, Nr]->f
        combo[[I]][Dataset=='WHI_afr']$Nr_SNPs<-a
        combo[[I]][Dataset=='JHS_afr']$Nr_SNPs<-d
        combo[[I]][Dataset=='UKB_afr']$Nr_SNPs<-b
        combo[[I]][Dataset=='HRS_afr']$Nr_SNPs<-f
        combo[[I]][Dataset=='HRS_eur']$Nr_SNPs<-f
        combo[[I]][, Std.PRS:=scale(PGS), by=Dataset]-> combo[[I]]
        cat(I, ' done\n')
}
do.call(rbind,combo)-> combo2
F2_x<-function(x, data){
#p(x|y)
y<-1-x
z<-quantile(data[,Res.Height], probs=x)[[1]]
z2<-quantile(data[,Std.PRS], probs=x)[[1]]
a<-nrow(data[Std.PRS>=z2])
b<-nrow(data[Std.PRS>=z2][Res.Height>=z])
res<-((b+0.5)/(a+0.5))/y #add 0.5 to b to avoid 0 in numerator. Haldane-Ascombe correction
#res<-log(res)
return(res)
}

F3_x<-function(x, data, data2){
#p(x|y)
y<-1-x
z<-quantile(data[,Res.Height], probs=x)[[1]]
z2<-quantile(data2[Dataset=='HRS_eur'][,Std.PRS], probs=x)[[1]]
a<-nrow(data[Std.PRS>=z2])
tp<-nrow(data[Std.PRS>=z2][Res.Height>=z])
fp<-nrow(data[Std.PRS<z2][Res.Height>=z])
tn<-nrow(data[Std.PRS<z2][Res.Height<z])
fn<-nrow(data[Std.PRS>=z2][Res.Height<z])
res<-list(F=((tp+0.5)/(a+0.5))/y, G=((fp+0.5)/(a+0.5))/y, f=tp/(tp+fn), g=fp/(fp+tn)) #add 0.5 to b to avoid 0 in numerator. Haldane-Ascombe correction
#res<-log(res)
return(res)
}

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='WHI_afr'],  data2=combo[[Y]]), Quantile=X))))-> AA
names(AA)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='WHI_afr']), Quantile=X))))-> AAA
names(AAA)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='JHS_afr'],  data2=combo[[Y]]), Quantile=X))))-> AJ
names(AJ)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='JHS_afr']), Quantile=X))))-> AAJ
names(AAJ)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='UKB_afr'],  data2=combo[[Y]]), Quantile=X))))-> AU
names(AU)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='UKB_afr']), Quantile=X))))-> AAU
names(AAU)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='HRS_eur'],  data2=combo[[Y]]), Quantile=X))))-> AP
names(AP)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='HRS_eur']), Quantile=X))))-> AAP
names(AAP)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F3_X=F3_x(x=X,data=combo[[Y]][Dataset=='HRS_afr'],  data2=combo[[Y]]), Quantile=X))))-> APE
names(APE)<-names(B_JHS)

lapply(1:length(names(B_JHS)), function(Y) do.call(rbind,lapply(seq(from=0.5,to=0.99, by=0.01), function(X) data.frame(F2_X=F2_x(x=X,data=combo[[Y]][Dataset=='HRS_afr']), Quantile=X))))-> AAPE
names(AAPE)<-names(B_JHS)

I<-names(AA)[63]

#png(paste0('~/height_prediction/figs_for_paper/figs/OR_WHI_', I,  ".png"))
plotA<-ggplot(AA[[I]], aes(x=Quantile, y=F3_X.F)) +
geom_point(size=2) + labs(title="WHI_afr", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12),legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()
#png(paste0('~/height_prediction/figs_for_paper/figs/OR_JHS_', I,  ".png"))
plotB<-ggplot(AJ[[I]], aes(x=Quantile, y=F3_X.F)) +
geom_point(size=2)  + labs(title="JHS_afr", y="OR")+geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12),legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()

#png(paste0('~/height_prediction/figs_for_paper/figs/OR_UKB_afr_', I,  ".png"))
plotC<-ggplot(AU[[I]], aes(x=Quantile, y=F3_X.F)) +
geom_point(size=2) + labs(title="UKB_afr", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12),legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()

#png(paste0('~/height_prediction/figs_for_paper/figs/OR_HRS_afr_', I,  ".png"))
plotD<-ggplot(APE[[I]], aes(x=Quantile, y=F3_X.F)) +
geom_point(size=2) + labs(title="HRS_afr", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12),legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()


plotE<-ggplot(AP[[I]], aes(x=Quantile, y=F3_X.F)) +
geom_point(size=2) + labs(title="HRS_eur", y="OR") + geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12),legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

###################################
#expectation
####################################
expected.or <- function(x, rho){
val <- qnorm(1-x)
prob <- pmvnorm(lower=c(val, val), upper=Inf, sigma=matrix(c(1,rho,rho,1),2,2))[1]
OR <- prob/x/x
return(OR)
}

x.vals <- seq(0.01,0.5, l=50)
rho125 <- sapply(x.vals, expected.or, sqrt(0.125))
rho025 <- sapply(x.vals, expected.or, sqrt(0.025))
rho038 <- sapply(x.vals, expected.or, sqrt(0.038))
rho041 <- sapply(x.vals, expected.or, sqrt(0.041))


my_col<-wes_palette("Darjeeling1",4)

#pdf('~/height_prediction/figs_for_paper/figs/OR_expected.pdf')
#plot(1-x.vals, rho12, type="l", xlim=c(0.5,1), col=my_col[1], bty="n", xlab="Quantile", ylab="Bivariate normal OR", lwd=2, cex.lab=1.5)
#lines(1-x.vals, rho041, type="l",col=my_col[2], lwd=2)
#lines(1-x.vals, rho036, type="l",col=my_col[3], lwd=2)
#lines(1-x.vals, rho034, type="l",col=my_col[4], lwd=2)
#abline(h=1, col='orange', lty=2)
#legend("topleft", c(expression(paste(R^2, "=0.12",sep="")), expression(paste(R^2,"=0.041", sep="")), expression(paste(R^2,"=0.036", sep="")), expression(paste(R^2,"=0.034", sep=""))), col=my_col, bty="n", lty=1, lwd=2)
#dev.off()

data.table(Name=c(rep("Rho_12.5",50), rep("Rho_2.5",50), rep("Rho_3.8",50) , rep("Rho_4.1", 50)), Bivariate_OR=c(rho125,rho025, rho038, rho041), Quantile=rep(1-x.vals,4))-> dt

plotF<-ggplot(dt, aes(x=Quantile, y=Bivariate_OR, colour=Name)) + geom_line() +
labs(title="Expectation", y="Expected OR") +  geom_hline(yintercept=1, linetype="dashed", color = "orange") +
theme(axis.title.y = element_text(size = 15), axis.title.x=element_text(size=15), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), legend.text=element_text(size=12),legend.position=c(.5,.75), legend.title=element_blank(),legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
coord_cartesian(ylim = c(1,8))
#ggsave('~/height_prediction/figs_for_paper/figs/OR_expected.pdf')

png('~/height_prediction/figs_for_paper/figs/OR_panel1.png', unit="in", res=300, height=12, width=7)
plot_grid(plotA, plotB,plotC, plotD, plotE, plotF,  labels = c("A", "B", "C","D", "E", "F"), nrow=3, align="v")
dev.off()
ggsave(paste0('~/height_prediction/figs_for_paper/figs/OR_panel1.pdf'))
