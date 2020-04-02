library(gamlss)
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
#df is the dataframe containing PValues
qexp=function(df){
  #number of rows in df
  nsnps=nrow(df)
  #generate P-values from chisquare statistic
  #sort df on P-value
  df=df%>%
    mutate(P=pnorm(RES,lower.tail = F))%>%
    arrange(P) %>% as.data.table
  #generate P-values from uniform distribution
  df=df%>%
    mutate(EXP_P=seq(1,nsnps)/(nsnps+1),
           EXP_NORM=qnorm(EXP_P,lower.tail=F)) %>% as.data.table
  return(df)
}
#
data <- fread('~/height_prediction/output/WHI_JHS_UKB_HRS_summ.txt')
data <- data[!is.na(data$HEIGHTX),][Prun_Set=='phys_100000_0.0005']
data[,":="(sex=ifelse(sex==1, "Male", sex))]
data[,":="(sex=ifelse(sex==2, "Female", sex))]
data[,":="(sex=ifelse(sex=="FEMALE", "Female", sex))]
data[, EUR_ANC:=ifelse(Dataset=='HRS_eur', 1, EUR_ANC)]
#residuals
model<-lm(HEIGHTX~sex+Dataset+sex*Dataset+age+age*Dataset+age*sex+age2+age2*Dataset+age2+sex, data=data)
data$RES <- scale(model$residuals) #take height residuals
#plot qq
dt_fem<-rbind(qexp(data[Dataset=='HRS_eur'][sex=='Female']),
        qexp(data[Dataset=='HRS_afr'][sex=='Female']),
        qexp(data[Dataset=='UKB_afr'][sex=='Female']),
        qexp(data[Dataset=='JHS_afr'][sex=='Female']),
        qexp(data[Dataset=='WHI_afr']))
dt_male<-rbind(qexp(data[Dataset=='HRS_eur'][sex=='Male']),
        qexp(data[Dataset=='HRS_afr'][sex=='Male']),
        qexp(data[Dataset=='UKB_afr'][sex=='Male']),
        qexp(data[Dataset=='JHS_afr'][sex=='Male']))
dt<-rbind(qexp(data[Dataset=='HRS_eur']),
        qexp(data[Dataset=='HRS_afr']),
        qexp(data[Dataset=='UKB_afr']),
        qexp(data[Dataset=='JHS_afr']),
        qexp(data[Dataset=='WHI_afr']))
as.factor(dt$Dataset)-> dt$Dataset
as.factor(dt_fem$Dataset)-> dt_fem$Dataset
as.factor(dt_male$Dataset)-> dt_male$Dataset
dt[, Cat:='ALL']
dt_fem[, Cat:='F']
dt_male[, Cat:='M']
dt2<-rbind(dt, dt_fem, dt_male)
ggplot(dt2,aes(x=EXP_NORM, y=RES, colour=Dataset))+ facet_wrap(~Cat) +
  geom_point()+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_bw()+
  labs(x=bquote(log[10]~"Expected"),
       y=bquote(log[10]~"Observed"))
ggsave('~/height_prediction/figs_for_paper/figs/qq_res_post_filter.pdf')
#prune out extreme values
model<-lm(HEIGHTX~sex+Dataset+age+age2+sex*Dataset+sex*age+sex*age2+Dataset*age+Dataset*age2, data=data)
data$RES <- scale(model$residuals) #take height residuals
ggplot(data, aes(y=RES,x=EUR_ANC, colour=Dataset)) + geom_point(alpha=0.5)
ggsave('~/height_prediction/figs_for_paper/figs/all_dtsets_res.pdf')
check<-data %>% group_by(Dataset) %>% summarise(Age=mean(age), Res.Height.Var=var(RES))
gamlss.mod1 <- gamlss(RES~EUR_ANC, family=NO2, data=data) #model residuals as a function of Eur Anc.Normal Distribution (With Variance As Sigma Parameter) For Fitting A GAMLSS
gamlss.mod2 <- gamlss(RES~EUR_ANC, sigma.formula=~EUR_ANC, family=NO2(sigma.link=identity), data=data)
summary(gamlss.mod2)

c1 <-  gamlss.mod2$mu.coefficients
c2 <- gamlss.mod2$sigma.coefficients
pd <- c(c2[1], -0.24*c2[1])


## abline(c1[1], c1[2], col="red", lwd=3)
x <- seq(0, 1, l=101)
sd.diff <- sqrt(c2[1]+x*c2[2])
pd.diff <- sqrt(pd[1]+x*pd[2])

mn <- c1[1] + x*c1[2]
#
data[, Group:=ifelse(Dataset=='HRS_eur', 'EUR', 'AFR')]
data[Dataset=='HRS_eur']$EUR_ANC<-jitter(data[Dataset=='HRS_eur']$EUR_ANC)

round(nrow(data[RES>-2.5][RES>=-2.5])/nrow(data),3) #percentage of points between 2.5 and -2.5 RES: 99.6%

ggplot(data, aes(x=EUR_ANC, y=RES)) + 
geom_point(cex=0.3, alpha=0.1, col='gray') + 
geom_segment(aes(x=x[1],y=mn[1], xend=x[101], yend=mn[101], colour="Mean Height")) + 
geom_segment(aes(x=x[1],y=(mn+1)[1], xend=x[101], yend=(mn+1)[101], colour="Mean Height"), lty=2) +
geom_segment(aes(x=x[1],y=(mn-1)[1], xend=x[101], yend=(mn-1)[101], colour="Mean Height"), lty=2) +
geom_segment(aes(x=x[1], y = (mn+sd.diff)[1], xend=x[101], yend=(mn+sd.diff)[101], colour='Fitted'), lty=2) + 
geom_segment(aes(x=x[1], y = (mn-sd.diff)[1], xend=x[101], yend=(mn-sd.diff)[101], colour='Fitted'), lty=2) + 
geom_segment(aes(x=x[1], y = (mn+pd.diff)[1], xend=x[101], yend=(mn+pd.diff)[101], colour='Expected'), lty=2) + 
geom_segment(aes(x=x[1], y = (mn-pd.diff)[1], xend=x[101], yend=(mn-pd.diff)[101], colour='Expected'), lty=2) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='bottom', legend.title=element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0.1, 0)) + labs(x="European Ancestry", y="Standardized residual height") + scale_color_manual(values=c("orange", "lightblue", "gray")) + coord_cartesian(xlim=c(0,1), ylim =c(-2.5, 2.5))
ggsave("~/height_prediction/figs_for_paper/figs/FigureS13.png")


