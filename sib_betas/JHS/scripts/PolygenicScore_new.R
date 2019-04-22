############################################
# A function to calculate polygenic scores##
############################################

#Input data: 
#1.a set of SNPs and their beta scores and p-values (and CHR, POS, perhaps)
#2.a population name (from the Phase 3 1000 Genomes Project)
#3.a method to calculate PGS.
#4. VCF file from 1000G

##Note:
#1. For now I am using the set of 697 SNPs defined based on p<5*10^-8 in the Wood et al (2014) paper. This set has already been LD-prunned in that paper. In the future, we may incorporate prunning methods here. Or, most likely, integration with LDpred or PLINK within this R function.
#2. Let's start with one of European ancestry, because the Wood et. al height GWAS was done in populations of European ancestry.
#3. For now I will use the only one I know how to implement (unsure if there are others or if differences rely on prunning -yes/no, clumping/prunning- and p-value thresholds only). 
#4. In this case we already have the vcf but we could add an option where vcf file from 1000G is read in, but that's too much memory so I filtered beforehand.

#************
#NOTE: WHI data and ubnphased genotypes
currentdir<-paste0(getwd(), "")
parentdir<-gsub("scripts/", "", parentdir)
PolScore4<- function(panel='WHI', tag='phys_100000_5e-08', CHR=22){
	if(panel=='JHS'){
	readRDS(paste0(parentdir, 'output/hei_', tag, '_v2.Rds'))-> hei	
	hei[[CHR]]-> hei2
        samps<-colnames(hei2)[9:(ncol(hei2)-6)]
	}
        hei2[ALT==Allele1]-> temp1
        hei2[REF==Allele1]-> temp2 #im ignoring the other two rows for now
        vector('list', length(samps))-> temp_list
        names(temp_list)<-samps
        if(nrow(temp1)>0 & nrow(temp2)>0){
               for(i  in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b]*2)-> temp_list[[i]]
                        temp_list[[i]] + sum(temp2[which(temp2[,i, with=F]=="0/0"),b]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b]*0)-> temp_list[[i]]
                }
        } else if (nrow(temp1)>0){
                for(i in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b]*2)-> temp_list[[i]]
                }
        } else if (nrow(temp2)>0){
                for(i in samps){
                        sum(temp2[which(temp2[,i, with=F]=="0/0"),b]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b]*0)-> temp_list[[i]]
                }
        }
#acollect sample names for this population from 1000G data.
        return(temp_list)
}

#******************
# Version 5
#******************

PolScore5<- function(panel='UKB', tag='phys_100000_5e-08', CHR=22){
        readRDS(paste0('/project/mathilab/bbita/gwas_admix/height/ukbiobank/ukb_data_run/hei_', tag, '.Rds'))-> hei
        hei[[CHR]][['afr']]-> hei2_afr
	hei[[CHR]][['eur']]-> hei2_eur
        samps_afr<-colnames(hei2_afr)[10:(ncol(hei2_afr)-10)]
	samps_eur<-colnames(hei2_eur)[10:(ncol(hei2_eur)-10)]
        #hei2<-hei[,c(1:5,10:7301)]
        hei2_afr[ALT==Allele1]-> temp1_afr
        hei2_afr[REF==Allele1]-> temp2_afr #im ignoring the other two rows for now
	hei2_eur[ALT==Allele1]-> temp1_eur
        hei2_eur[REF==Allele1]-> temp2_eur #im ignoring the other two rows for now
        vector('list', length(samps_afr))-> temp_list_afr
        names(temp_list_afr)<-samps_afr
	vector('list', length(samps_eur))-> temp_list_eur
        names(temp_list_eur)<-samps_eur
        if(nrow(temp1_afr)>0 & nrow(temp2_afr)>0){
               for(i  in samps_afr){
                        sum(temp1_afr[which(temp1_afr[,i, with=F]=="0/0"),b]*0) + sum(temp1_afr[which(temp1_afr[,i,with=F]=="0/1"),b]*1) + sum(temp1_afr[which(temp1_afr[,i,with=F]=="1/0"),b]*1) + sum(temp1_afr[which(temp1_afr[,i,with=F]=="1/1"),b]*2)-> temp_list_afr[[i]]
                        temp_list_afr[[i]] + sum(temp2_afr[which(temp2_afr[,i, with=F]=="0/0"),b]*2) + sum(temp2_afr[which(temp2_afr[,i,with=F]=="0/1"),b]*1) + sum(temp2_afr[which(temp2_afr[,i,with=F]=="1/0"),b]*1) + sum(temp2_afr[which(temp2_afr[,i,with=F]=="1/1"),b]*0)-> temp_list_afr[[i]]
                }
        } else if (nrow(temp1_afr)>0){
                for(i in samps_afr){
                        sum(temp1_afr[which(temp1_afr[,i, with=F]=="0/0"),b]*0) + sum(temp1_afr[which(temp1_afr[,i,with=F]=="0/1"),b]*1) + sum(temp1_afr[which(temp1_afr[,i,with=F]=="1/0"),b]*1) + sum(temp1_afr[which(temp1_afr[,i,with=F]=="1/1"),b]*2)-> temp_list_afr[[i]]
                }
	
        } else if (nrow(temp2_afr)>0){
                for(i in samps_afr){
                        sum(temp2_afr[which(temp2_afr[,i, with=F]=="0/0"),b]*2) + sum(temp2_afr[which(temp2_afr[,i,with=F]=="0/1"),b]*1) + sum(temp2_afr[which(temp2_afr[,i,with=F]=="1/0"),b]*1) + sum(temp2_afr[which(temp2_afr[,i,with=F]=="1/1"),b]*0)-> temp_list_afr[[i]]
                }
        }
       if(nrow(temp1_eur)>0 & nrow(temp2_eur)>0){
               for(i  in samps_eur){
                        sum(temp1_eur[which(temp1_eur[,i, with=F]=="0/0"),b]*0) + sum(temp1_eur[which(temp1_eur[,i,with=F]=="0/1"),b]*1) + sum(temp1_eur[which(temp1_eur[,i,with=F]=="1/0"),b]*1) + sum(temp1_eur[which(temp1_eur[,i,with=F]=="1/1"),b]*2)-> temp_list_eur[[i]]
                        temp_list_eur[[i]] + sum(temp2_eur[which(temp2_eur[,i, with=F]=="0/0"),b]*2) + sum(temp2_eur[which(temp2_eur[,i,with=F]=="0/1"),b]*1) + sum(temp2_eur[which(temp2_eur[,i,with=F]=="1/0"),b]*1) + sum(temp2_eur[which(temp2_eur[,i,with=F]=="1/1"),b]*0)-> temp_list_eur[[i]]
                }
	} else if (nrow(temp1_eur)>0){
                for(i in samps_eur){
                        sum(temp1_eur[which(temp1_eur[,i, with=F]=="0/0"),b]*0) + sum(temp1_eur[which(temp1_eur[,i,with=F]=="0/1"),b]*1) + sum(temp1_eur[which(temp1_eur[,i,with=F]=="1/0"),b]*1) + sum(temp1_eur[which(temp1_eur[,i,with=F]=="1/1"),b]*2)-> temp_list_eur[[i]]
                }

        } else if (nrow(temp2_eur)>0){
                for(i in samps_eur){
                        sum(temp2_eur[which(temp2_eur[,i, with=F]=="0/0"),b]*2) + sum(temp2_eur[which(temp2_eur[,i,with=F]=="0/1"),b]*1) + sum(temp2_eur[which(temp2_eur[,i,with=F]=="1/0"),b]*1) + sum(temp2_eur[which(temp2_eur[,i,with=F]=="1/1"),b]*0)-> temp_list_eur[[i]]
                }
        }
#acollect sample names for this population from 1000G data.
        return(list(afr=temp_list_afr, eur=temp_list_eur))
}

#*******
#* END *
#*******
