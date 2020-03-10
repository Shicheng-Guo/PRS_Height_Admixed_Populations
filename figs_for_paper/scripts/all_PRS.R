#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(reshape2)
library(asbio)
library(ggplot2)
library(dplyr)
########################
PRS1_all_res<-readRDS(file='~/height_prediction/loc_anc_analysis/output/PRS1_all_res.Rds')
PRS2_all_res<-readRDS(file='~/height_prediction/loc_anc_analysis/output/PRS2_all_res.Rds')
PRS1_2_all_res<-readRDS(file='~/height_prediction/loc_anc_analysis/output/PRS1_2_all_res.Rds')

PRS1_2_all_res[,Dataset2:=ifelse(Dataset=='WHI', "Women's Health Initiative", "Jackson Heart Study + \nHealth and Retirement Study")]
PRS1_2_all_res[,Dataset:=NULL]
dt<-melt(PRS1_2_all_res,id=c('Alpha', 'Dataset2', 'Test'))

ggplot(dt, aes(x=Alpha, y=value,colour=Test)) + facet_wrap(~Dataset2) +
geom_line(size=1.2) + coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.048)) + scale_color_manual(values=c("darkseagreen4", "darkslateblue", "darkorange3")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial~R^2))

ggsave('~/height_prediction/loc_anc_analysis/figs/PRS1_2.pdf')


PRS3<-readRDS('~/height_prediction/loc_anc_analysis/output/part_R2_WHI.Rds')
#rbind(LA, data.table(Alfa=LA$Alfa, part_R2=NA,Dataset='JHS+HRS'))-> LA #for now
test4[, PRS3:=c(unlist(LA[1:101])/100, rep(NA, 101))]
#test5[alpha<=0.5]-> test5

#test5[, Dataset2:=ifelse(Dataset=='WHI', "Women's Health Initiative", "Jackson Heart Study + \nHealth and Retirement Study")]


dt[,PRS1:=PRS1/test4$PRS1[1]]
dt[,PRS2:=PRS2/test4$PRS2[1]]
dt[,PRS3:=PRS3/test4$PRS3[1]]
melt(dt,id=c('Dataset', 'alpha'))->dt2
dt2[Dataset=='WHI']-> dt2
dt2[, Dataset2:=ifelse(Dataset=='WHI', "Women's Health Initiative", "Jackson Heart Study + \nHealth and Retirement Study")]

ggplot(test5, aes(x=alpha, y=value,colour=variable))  + 
#facet_wrap(~Dataset2) + 
geom_line(size=1.2) + 
#coord_cartesian(xlim=c(0,0.5), ylim=c(0.02, 0.048)) + 
scale_color_manual(values=c("darkseagreen4", "darkslateblue", "deeppink4", "gray7")) +
geom_hline(yintercept=0.0412) +
theme(strip.text.x = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial~R^2))
ggsave('~/height_prediction/loc_anc_analysis/figs/multi_prs.pdf')

ggplot(dt2, aes(x=alpha, y=value,colour=variable))  +
#facet_wrap(~Dataset2) +
geom_line(size=1.2) +
#coord_cartesian(xlim=c(0,0.5), ylim=c(0, 2)) +
scale_color_manual(values=c("darkseagreen4", "darkslateblue", "deeppink4", "gray7")) +
geom_hline(yintercept=1) +
theme(strip.text.x = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 18), axis.title.x=element_text(size=18),axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.key=element_blank(), legend.background=element_blank(),legend.title=element_blank(), legend.text=element_text(size=18)) + labs(x=expression(alpha), y=expression(Partial~R^2~fold))
ggsave('~/height_prediction/loc_anc_analysis/figs/multi_prs_fold.pdf')

txtStop()
