#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#**************************************
#*      CALCULATE POLYGENIC SCORES   **
#**************************************
source('~/height_prediction/scripts/PolygenicScore_v2.R')
library("optparse")
library(data.table)
library(dplyr)

library(biomaRt)
library(parallel)
#args<-c("LD_prun","LD_250000_0.01_0.5")
options(scipen=10)
options(digits=10)
home='~/height_prediction/'

if(args[1]=='sib_betas'){
         #read in betas from sibling analyses
	fread(paste0(home, args[1], "/input/sibestimates_50.tsv"))->ukb_height_sib
        ukb_height_sib[,SE:= -abs(beta)/qnorm(pval/2)]
        ukb_height_sib[, MarkerName:=SNP]
        ukb_height_sib[,CHR:=as.integer(CHR)]
        ukb_height_sib[,POS:=as.integer(POS)]
        ukb_height_sib[order(CHR,POS)]-> ukb_height_sib
        fread(paste0('zcat ', home,"gwas/input/50_raw_filtered.txt.gz"))-> ukb_height #read in GWAS summary statistics for height from the Uk Biobank
        ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)] #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
	na.omit(ukb_height)-> ukb_height
	ukb_height[CHR %in% 1:22]-> ukb_height
	gc()
	ukb_height[, N:=n_complete_samples][, AC:=NULL][, b:=beta][,p:=pval]
        ukb_height[,n_complete_samples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
        ukb_height[,.(Allele1,Allele2,SE, p, N, CHR, POS)]-> ukb_height
	ukb_height[,POS:=as.integer(POS)]
	ukb_height[,CHR:=as.integer(CHR)]
        ukb_height[order(CHR,POS)]-> ukb_height
	setkey(ukb_height,CHR, POS)
	setkey(ukb_height_sib,CHR, POS)
        ukb_height[ukb_height_sib, nomatch=0][,.(Allele1,Allele2,beta, SE, p, N, CHR, POS)][, b:=beta][,beta:=NULL]-> ukb_height_sib
	ukb_height_sib[!(POS %in% ukb_height_sib[which(duplicated(ukb_height_sib, by='POS')), POS])]-> ukb_height
} else if (args[1]=='gwas'){
      	fread(paste0('zcat ', home,"gwas/input/50_raw_filtered.txt.gz"))-> ukb_height #read in GWAS summary statistics for height from the Uk Biobank
        ukb_height[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
        na.omit(ukb_height)-> ukb_height
        ukb_height[, N:=n_complete_samples][, AC:=NULL][, b:=beta][,p:=pval]
        ukb_height[,n_complete_samples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
        ukb_height[,.(Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb_height
	ukb_height[CHR %in% seq(1:22)]-> ukb_height
        ukb_height[,POS:=as.integer(POS)]
	ukb_height[,CHR:=as.integer(CHR)]
        ukb_height[order(CHR,POS)]-> ukb_height
	ukb_height[!(POS %in% ukb_height[which(duplicated(ukb_height, by='POS')), POS])]-> ukb_height
}

cat('checkpoint\n')

if(!(args[3] %in% c("LD_50000_0.01_0.5", "LD_100000_0.01_0.5", "LD_250000_0.01_0.5"))){
	readRDS(paste0(home, args[1], '/', args[2], '/output/hei_',args[3], '.Rds'))-> snp_list #prunend based on UKB
	lapply(snp_list, function(X) X[, CHR:=CHROM][, CHROM:=NULL])-> snp_list
	lapply(snp_list, function(X) setkey(X, CHR, POS))
	setkey(ukb_height, CHR, POS)
	lapply(snp_list, function(X) X[ukb_height, nomatch=0])-> snp_list
} else{
	cat('SOMETHING\n')
	snp_list<- lapply(1:22, function(X) data.table(MarkerName=readRDS(paste0(home, args[1], '/', args[2], '/prunned_1kg/LD_prunned_hei_chr', X, '_',args[3], '.Rds'))[['keep']]))
	betas<-fread(paste0(home, args[1], '/', args[2],"/plink2/output/betas_for_plink.txt"))
	colnames(betas)[1]<-'MarkerName'
	setkey(betas, MarkerName)
	#grch37_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_snp")
	#ID_posgrch37<-lapply(1:22, function(X) getBM(attributes = c('refsnp_id','allele', 'chr_name','chrom_start'),filters = 'snp_filter',values = snp_list[[X]]$MarkerName,mart = grch37_snp))
	#lapply(ID_posgrch37, function(X) setDT(X))
	for(CR in 22:1){
		#colnames(ID_posgrch37[[CR]])<-c('MarkerName','Allele','CHR','POS')
		#setkey(ID_posgrch37[[CR]], MarkerName, CHR, POS)
		#ID_posgrch37[[CR]][,CHR:=as.integer(CHR)]
		#ID_posgrch37[[CR]][,POS:=as.integer(POS)]
		#setkey(ukb_height, CHR, POS)
		#setkey(ID_posgrch37[[CR]], CHR, POS)
		betas_chr<-merge(betas, snp_list[[CR]])
		setkey(ukb_height, CHR, POS)
		setkey(betas_chr, CHR, POS)
		#ukb_height[betas_chr, nomatch=0][, i.MarkerName:=NULL]-> ukb_height2
		betas_chr[order(CHR,POS)]-> snp_list[[CR]]
		remove(betas_chr)
		gc()
		cat(CR, '\n')
	}
	remove(betas);gc();
	names(snp_list)<-seq(1:22)
	saveRDS(snp_list, file=paste0(home, args[1], '/', args[2], '/output/hei_', args[3], '.Rds'))
}
cat('Got SNP IDs, CHR and POS.\n')
#
if((args[3] %in% c("LD_50000_0.01_0.5", "LD_100000_0.01_0.5", "LD_250000_0.01_0.5"))){
       hei<-vector('list', 22)
       names(hei)<-seq(1:22)
       for(Z in 22:1){
             	fread(paste0('zcat ',home, args[1], '/', args[2],'/output/hei_SNPs_chr', Z, '.vcf.gz'), header=T)-> hei[[Z]]
		gc()
                colnames(hei[[Z]])<-unlist(fread(paste0(home, 'input/', args[2], '/header_', args[2], '.txt'), header=F, sep="\t"))
                colnames(hei[[Z]])[1]<-'CHR'
		colnames(hei[[Z]])[3]<-'MarkerName'
		unique(hei[[Z]], by=c('CHR','POS'))-> hei[[Z]]
		gc()
                hei[[Z]][REF %in% c("A","C","G","T")]-> hei[[Z]]
		gc()
                setkey(snp_list[[Z]],CHR, POS)
                setkey(hei[[Z]], CHR, POS)
		gc()
                hei[[Z]][snp_list[[Z]], nomatch=0]-> hei[[Z]]
		hei[[Z]][,i.MarkerName:=NULL]
		hei[[Z]][,b:=BETA][,BETA:=NULL]
		cat(Z, ' done \n')
		gc()
	}
        saveRDS(hei, file=paste0(home, args[1], '/', args[2], '/output/hei_', args[3], '_v2.Rds'))
	cat('another checkpoint\n')
} else{
        saveRDS(snp_list, file=paste0(home, args[1], '/', args[2],'/output/hei_', args[3], '_v2.Rds'))
	snp_list-> hei
	remove(snp_list)
	gc()
}
cat('checkpoint 3\n')
PGS<-vector('list',22)
names(PGS)<-c(1:22)
#chose those that have no snps:
skip_me<-which(unlist(lapply(1:22, function(X) nrow(hei[[X]])))==0)
if(length(skip_me)>=1){
for (CR in seq(1:22)[-skip_me]){
	print(paste0("Chromosome is ", CR))
        PolScore2(CHR=CR, panel=args[1], panel2=args[2], tag=args[3])-> PGS[[CR]]
        cat(paste0(CR, '  done\n'))
}
}else if(length(skip_me)==0){
	for (CR in seq(1:22)){
        print(paste0("Chromosome is ", CR))
        PolScore2(CHR=CR, panel=args[1], panel2=args[2], tag=args[3])-> PGS[[CR]]
        cat(paste0(CR, '  done\n'))
}
}
samps<-names(PGS[[1]]) #sum PGS across chromosomes.
PGS2<-vector('list', length(samps))
names(PGS2)<-samps
if(length(skip_me)==1){
        PGS[[skip_me]]<-rep(NA, length(samps))
        names(PGS[[skip_me]])<-samps
}
if(length(skip_me)>1){
        for(I in skip_me){
                PGS[[skip_me[I]]]<-rep(NA, length(samps))
                names(PGS[[skip_me[I]]])<-samps
}
}
saveRDS(PGS,file=paste0(home, args[1], '/', args[2], '/output/PGS_chrs_', args[3], '.Rds'))

for (S in samps){
	sum(PGS[[1]][[S]],PGS[[2]][[S]],PGS[[3]][[S]],PGS[[4]][[S]],PGS[[5]][[S]],PGS[[6]][[S]],PGS[[7]][[S]],PGS[[8]][[S]],PGS[[9]][[S]],PGS[[10]][[S]],PGS[[11]][[S]],PGS[[12]][[S]],PGS[[13]][[S]],PGS[[14]][[S]],PGS[[15]][[S]],PGS[[16]][[S]],PGS[[17]][[S]],PGS[[18]][[S]],PGS[[19]][[S]],PGS[[20]][[S]],PGS[[21]][[S]],PGS[[22]][[S]], na.rm=T)->PGS2[[S]]
        cat(paste0(S, ' done\n'))
        }
saveRDS(PGS2, file=paste0(home, args[1], '/', args[2],'/output/PGS_', args[2], '_', args[3], '.Rds'))

#TheEnd

