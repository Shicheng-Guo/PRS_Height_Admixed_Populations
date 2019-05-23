library(readr)
library(data.table)
library(reshape2)
library(asbio)
library(ggplot2)
#########################
PolScore<- function(beta='all'){
        samps<-colnames(prun2)[9:(ncol(prun2)-11)]
        setDT(prun2)
        #hei2<-hei[,c(1:5,10:7301)]
        prun2[ALT==Allele1]-> temp1
        prun2[REF==Allele1]-> temp2 #im ignoring the other two rows for now
        vector('list', length(samps))-> temp_list
        names(temp_list)<-samps
	cat('checkpoint \n')
	if(beta=='all'){
        if(nrow(temp1)>0 & nrow(temp2)>0){
               for(i  in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b3]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b3]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b3]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b3]*2)-> temp_list[[i]]
                        temp_list[[i]] + sum(temp2[which(temp2[,i, with=F]=="0/0"),b3]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b3]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b3]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b3]*0)-> temp_list[[i]]
                }
        } else if (nrow(temp1)>0){
                for(i in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b3]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b3]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b3]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b3]*2)-> temp_list[[i]]
                }
        } else if (nrow(temp2)>0){
                for(i in samps){
                        sum(temp2[which(temp2[,i, with=F]=="0/0"),b3]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b3]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b3]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b3]*0)-> temp_list[[i]]
                }
        }
	  cat('checkpoint 2\n')
	}
	if(beta==1){
	        if(nrow(temp1)>0 & nrow(temp2)>0){
               for(i  in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b1]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b1]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b1]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b1]*2)-> temp_list[[i]]
                        temp_list[[i]] + sum(temp2[which(temp2[,i, with=F]=="0/0"),b1]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b1]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b1]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b1]*0)-> temp_list[[i]]
                }
        } else if (nrow(temp1)>0){
                for(i in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b1]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b1]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b1]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b1]*2)-> temp_list[[i]]
                }
        } else if (nrow(temp2)>0){
                for(i in samps){
                        sum(temp2[which(temp2[,i, with=F]=="0/0"),b1]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b1]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b1]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b1]*0)-> temp_list[[i]]
                }
        }
	cat('checkpoint 3\n')
	}
	if(beta==2){
                if(nrow(temp1)>0 & nrow(temp2)>0){
               for(i  in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b2]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b2]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b2]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b2]*2)-> temp_list[[i]]
                        temp_list[[i]] + sum(temp2[which(temp2[,i, with=F]=="0/0"),b2]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b2]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b2]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b2]*0)-> temp_list[[i]]
                }
        } else if (nrow(temp1)>0){
                for(i in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b2]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b2]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b2]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b2]*2)-> temp_list[[i]]
                }
        } else if (nrow(temp2)>0){
                for(i in samps){
                        sum(temp2[which(temp2[,i, with=F]=="0/0"),b2]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b2]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b2]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b2]*0)-> temp_list[[i]]
                }
        }
	}
#acollect sample names for this population from 1000G data.
        return(temp_list)
}

#########################

###########
rsq.R2 <- function(formula1, formula2, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula1, data=d)
  fit2<- lm(formula2, data=d)
  res<-partial.R2(fit, fit2)
  return(res)
}
##############

res<-vector('list', 22)

for(args in 1:22){

#args<-22
	what <- paste0("~/height_prediction/input/ukb_afr/UKB_kgCY_chr", args)

	betas<-fread(paste0('~/height_prediction/gwas/ukb_afr/output/AS_Beta_chr', args,'example.txt'))

	prun<-setDT(readRDS('~/height_prediction/gwas/WHI/output/hei_phys_100000_0.0005_v2.Rds')[[args]])

	fread(paste0('~/height_prediction/input/ukb_afr/UKB_kgCY_chr',args, '.phsnp'))-> snps
	colnames(snps)<-c("ID", "CHR", "V3", "POS","REF","ALT")

	tag='phys_100000_0.0005'

	cbind(betas, snps)-> betas2
	data<-betas2[which(betas2[, POS] %in% prun[, POS]),]
	prun[POS %in% betas2[, POS]]-> prun2
	#data$POS==prun2$POS
	prun2[, b1:=data$POP1]
	prun2[, b2:=data$POP2]
	prun2[, b3:=data$ALL]

	res[[args]]<-vector('list', 3)

	test<-PolScore(beta=1)
	test2<-PolScore(beta=2)
	test3<-PolScore(beta='all')
	res[[args]][[1]]<-test
	res[[args]][[2]]<-test2
	res[[args]][[3]]<-test3
	cat('chr ', args, ' done\n')
}

saveRDS(res, file="~/height_prediction/gwas/WHI/output/all_prs.Rds")

data.table(SUBJID=names(test), 
PRS_POP1=(unlist(res[[1]][[1]])+unlist(res[[2]][[1]])+unlist(res[[3]][[1]])+unlist(res[[4]][[1]])+unlist(res[[5]][[1]])+unlist(res[[6]][[1]])+unlist(res[[7]][[1]])+unlist(res[[8]][[1]])+unlist(res[[9]][[1]])+unlist(res[[10]][[1]])+unlist(res[[11]][[1]])+unlist(res[[12]][[1]])+unlist(res[[13]][[1]])+unlist(res[[14]][[1]])+unlist(res[[15]][[1]])+unlist(res[[16]][[1]])+unlist(res[[17]][[1]])+ unlist(res[[18]][[1]])+unlist(res[[19]][[1]])+unlist(res[[20]][[1]])+unlist(res[[21]][[1]])+unlist(res[[22]][[1]])), 
PRS_POP2=(unlist(res[[1]][[2]])+unlist(res[[2]][[2]])+unlist(res[[3]][[2]])+unlist(res[[4]][[2]])+unlist(res[[5]][[2]])+unlist(res[[6]][[2]])+unlist(res[[7]][[2]])+unlist(res[[8]][[2]])+unlist(res[[9]][[2]])+unlist(res[[10]][[2]])+unlist(res[[11]][[2]])+unlist(res[[12]][[2]])+unlist(res[[13]][[2]])+unlist(res[[14]][[2]])+unlist(res[[15]][[2]])+unlist(res[[16]][[2]])+unlist(res[[17]][[2]])+ unlist(res[[18]][[2]])+unlist(res[[19]][[2]])+unlist(res[[20]][[2]])+unlist(res[[21]][[2]])+unlist(res[[22]][[2]])), 
PRS_all=(unlist(res[[1]][[3]])+unlist(res[[2]][[3]])+unlist(res[[3]][[3]])+unlist(res[[4]][[3]])+unlist(res[[5]][[3]])+unlist(res[[6]][[3]])+unlist(res[[7]][[3]])+unlist(res[[8]][[3]])+unlist(res[[9]][[3]])+unlist(res[[10]][[3]])+unlist(res[[11]][[3]])+unlist(res[[12]][[3]])+unlist(res[[13]][[3]])+unlist(res[[14]][[3]])+unlist(res[[15]][[3]])+unlist(res[[16]][[3]])+unlist(res[[17]][[3]])+ unlist(res[[18]][[3]])+unlist(res[[19]][[3]])+unlist(res[[20]][[3]])+unlist(res[[21]][[3]])+unlist(res[[22]][[3]])), 
PRS_EUR=unlist(readRDS('~/height_prediction/gwas/WHI/output/PGS_WHI_phys_100000_0.0005.Rds')))-> a


fread('~/height_prediction/input/WHI/WHI_phenotypes.txt')-> Pheno_WHI

Pheno_WHI[, SUBJID:=paste0("0_", as.character(Pheno_WHI[, SUBJID]))]
setkey(Pheno_WHI, SUBJID)
#add ancestry

anc_WHI<-cbind(fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.2.Q'), fread('~/height_prediction/input/WHI/WHI_b37_strand_prune_include.fam')[,V2])
colnames(anc_WHI)<-c("AFR_ANC","EUR_ANC","SUBJID")
anc_WHI[,SUBJID:=paste0("0_", SUBJID)]
setkey(anc_WHI, SUBJID)
setkey(a, SUBJID)
setkey(Pheno_WHI, SUBJID)

a[Pheno_WHI][anc_WHI]-> final

final[,c("SUBJID","PRS_POP1","PRS_POP2","PRS_all", "PRS_EUR", "SEX", "HEIGHTX", "EUR_ANC","AFR_ANC","AGE")]-> final2

#ok, so for ukb_afr and WHI (both), POP1 is AFR
cat('another checkpoint\n')

final2[,AGE2:=AGE^2][, PRS_eur_afr:=((mean(EUR_ANC)*scale(PRS_EUR))+(mean(AFR_ANC)*scale(PRS_POP1)))]


partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_EUR, data=final2)) #4.1%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_all, data=final2)) #0.005676512%
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP1, data=final2)) #0.003577538 %
partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_POP2, data=final2)) #0.03868322%

partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(HEIGHTX~AGE+AGE2+EUR_ANC+PRS_eur_afr, data=final2)) #0.65%


a_vec<-seq(0,1,0.1)
wei_PRS<-vector('list', length(a_vec))
names(wei_PRS)<-a_vec
part_r2<-c()
for(i in 1:length(a_vec)){

cat(i, '\n')
wei_PRS[[i]]<-(a_vec[i]*scale(final2$PRS_EUR))+((1-a_vec[i])*scale(final2$PRS_POP1))
part_r2[i]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+wei_PRS[[i]]))
}

data.table(part_R2=part_r2, alfa=a_vec)-> wanna_plot

ggplot(wanna_plot, aes(x=alfa, y=part_R2)) + geom_point(size=0.8)+ geom_line()

ggsave('~/height_prediction/gwas/WHI/figs/alfa_part_R2_WHI.pdf')


a_vec2<-seq(from=-20, to=20, by=0.1)

wei_PRS_2<-vector('list', length(a_vec2))
names(wei_PRS)<-a_vec
part_r2_2<-c()


for(i in 1:length(a_vec2)){

cat(i, '\n')
wei_PRS_2[[i]]<-scale(final2$PRS_EUR)+((final2$AFR_ANC)*(a_vec2[i])*scale(final2$PRS_POP1))
part_r2_2[i]<-partial.R2(lm(HEIGHTX~AGE+AGE2+EUR_ANC, data=final2), lm(final2$HEIGHTX~final2$AGE+final2$AGE2+final2$EUR_ANC+wei_PRS_2[[i]]))
}

data.table(part_R2=part_r2_2, kappa=a_vec2)-> wanna_plot2

ggplot(wanna_plot2, aes(x=kappa, y=part_R2)) + geom_point(size=0.8)+ geom_line()

ggsave('~/height_prediction/gwas/WHI/figs/kappa_alfa_part_R2_WHI.pdf')
