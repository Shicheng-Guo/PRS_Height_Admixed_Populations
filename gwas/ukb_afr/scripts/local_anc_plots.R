library(data.table)
library(reshape2)
library(ggplot2)

#a<-c("100", "200", "300", "400","500","600","700", "800", "900", "1000","1500", "2000","2500", "3000","3500", "4000","4500","5000","5500", "6000", "6500", "7000", "7500", "8000", "900", "9500", "10000","12500", "15000", "17500", "20000", "30000","40000","50000", "75000","100000","200000", "300000", "400000","500000","1000000")

#this<-c("100","200", "300", "400","500","600","700", "800", "900", "1000","1500", "2000","2500","3500", "4000","4500","5000","5500", "6000", "6500", "8000", "9000", "9500", "10000","12500", "15000", "17500", "20000", "30000","40000","50000", "75000","100000","200000", "300000", "400000","500000","1000000")

#this<-c("100","200","1000","2000", "10000", "20000", "100000", "200000")
this<-c("100", "200","1000","10000","20000", "100000", "200000")
all<-lapply(paste0("~/height_prediction/gwas/ukb_afr/outpout/betas_phys_100000_0.0005_", this, ".Rds"), function(X) setDT(readRDS(X)))
names(all)<-this

#all[[9]][,y2:=NULL]
sapply(1:length(all), function(X) all[[X]][, Win_Size:=this[X]])
for(I in 1:length(all)){
a<-coef(lm(y2~CEU_YRI_diff.rate, data=all[[I]]))[[2]];
a1<-confint(lm(y2~CEU_YRI_diff.rate, data=all[[I]]))[2,]
bla<-coef(lm(y2~AA.rate, data=all[[I]]))[[2]]
b1<-confint(lm(y2~AA.rate, data=all[[I]]))[2,]
d<-coef(lm(y2~CEU.rate, data=all[[I]]))[[2]]
d1<-confint(lm(y2~CEU.rate, data=all[[I]]))[2,]
f<- coef(lm(y2~YRI.rate, data=all[[I]]))[[2]]
f1<- confint(lm(y2~YRI.rate, data=all[[I]]))[2,]
g<-coef(lm(y2~COMBINED.rate, data=all[[I]]))[[2]]
g1<-confint(lm(y2~COMBINED.rate, data=all[[I]]))[2,]
all[[I]][,Slope_diff:=a];
all[[I]][,Slope_YRI:=f]
all[[I]][,Slope_AA:=bla];
all[[I]][, Slope_CEU:=d]
all[[I]][,Slope_Combined:=g]
all[[I]][,`:=`(c1_diff=a1[1], c2_diff=a1[2])]
all[[I]][,`:=`(c1_AA=b1[1], c2_AA=b1[2])]
all[[I]][,`:=`(c1_CEU=d1[1], c2_CEU=d1[2])]
all[[I]][,`:=`(c1_YRI=f1[1], c2_YRI=f1[2])]
all[[I]][,`:=`(c1_Comb=g1[1], c2_Comb=g1[2])]
cat(I,"\r")
}

do.call(rbind, all)-> all2
#all2<-all2[,.(AA.rate, CEU.rate, YRI.rate, CEU_YRI_diff.rate, y, Win_Size, Slope_diff, Slope_YRI, Slope_CEU, Slope_AA)]
all2[,AA.rate:c2_Comb]-> all2

melt(all2, id.vars=c("AA.rate","CEU_YRI_diff.rate","CEU.rate","YRI.rate","COMBINED.rate","y2","y", "b", "Win_Size","c1_diff", "c2_diff", "c1_AA", "c2_AA", "c1_CEU", "c2_CEU", "c1_YRI", "c2_YRI", "c1_Comb", "c2_Comb"))-> all3

all3[, Win_Size:=as.numeric(Win_Size)]

rbind(all3[variable=="Slope_diff"][,c1:=c1_diff][, c2:=c2_diff],
all3[variable=="Slope_AA"][,c1:=c1_AA][, c2:=c2_AA],
all3[variable=="Slope_YRI"][,c1:=c1_YRI][, c2:=c2_YRI],
all3[variable=="Slope_CEU"][,c1:=c1_CEU][, c2:=c2_CEU],
all3[variable=="Slope_Combined"][, c1:=c1_Comb][,c2:=c2_Comb])-> all3

#all3[Win_Size<=600000]-> all3
pdf('/height_prediction/gwas/ukb_afr/figs/local_anc.pdf')

ggplot(all3, aes(x=Win_Size,y=value, group=variable, color=variable)) + geom_point()+ geom_line() + scale_x_log10()
dev.off()


pdf('/height_prediction/gwas/ukb_afr/figs/local_anc_v2.pdf')

ggplot(all3, aes(x=Win_Size,y=value, group=variable, color=variable)) + geom_errorbar(aes(x=Win_Size, group=variable, color=variable,ymin=c1, ymax=c2), width=.2) + geom_point()+ geom_line() + scale_x_log10()
dev.off()
