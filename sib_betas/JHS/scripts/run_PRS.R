#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#**************************************
#*      CALCULATE POLYGENIC SCORES   **
#**************************************
source('PolygenicScore_new.R')
library("optparse")
library(data.table)
library(dplyr)
library(biomaRt)
library(parallel)
#args<-c("LD_prun","LD_250000_0.01_0.5")

currentdir<-paste0(getwd(),"/")
parentdir<-gsub("scripts/", "", parentdir)

cat('Tag is ', args[2], '\n')

fread('/project/mathilab/bbita/gwas_admix/height_prediction/sib_betas/input/sib_beta_gwas_p.txt')-> ukb_height #read in betas from sib analyses
ukb_height[,CHR:=as.integer(CHR)]
ukb_height[,POS:=as.integer(POS)]
ukb_height[order(CHR,POS)]-> ukb_height
cat('checkpoint\n')

if(!(args[2] %in% c("LD_50000_0.01_0.5", "LD_100000_0.01_0.5", "LD_250000_0.01_0.5"))){
	readRDS(paste0(parentdir, 'output/hei_',args[2], '.Rds'))-> snp_list #prunend based on 1KG
	lapply(snp_list, function(X) X[, CHR:=CHROM][, CHROM:=NULL])-> snp_list
	lapply(snp_list, function(X) setkey(X, CHR, POS))
	setkey(ukb_height, CHR, POS)
	lapply(snp_list, function(X) X[ukb_height, nomatch=0])-> snp_list
} else{
	cat('SOMETHING\n')
	snp_list<- lapply(1:22, function(X) data.table(MarkerName=readRDS(paste0(parentdir, 'prunned_1kg/LD_prunned_hei_chr', X, '_',args[2], '.Rds'))[['keep']]))
	grch37_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_snp")
	ID_posgrch37<-lapply(1:22, function(X) getBM(attributes = c('refsnp_id','allele', 'chr_name','chrom_start'),filters = 'snp_filter',values = snp_list[[X]]$MarkerName,mart = grch37_snp))
	lapply(ID_posgrch37, function(X) setDT(X))
	for(CR in 1:22){
		colnames(ID_posgrch37[[CR]])<-c('MarkerName','Allele','CHR','POS')
		setkey(ID_posgrch37[[CR]], MarkerName, CHR, POS)
		ID_posgrch37[[CR]][,CHR:=as.integer(CHR)]
		ID_posgrch37[[CR]][,POS:=as.integer(POS)]
		setkey(ukb_height, MarkerName, CHR, POS)
		setkey(ID_posgrch37[[CR]], MarkerName, CHR, POS)
		ukb_height[ID_posgrch37[[CR]], nomatch=0]-> ukb_height2
		ukb_height2[order(CHR,POS)]-> snp_list[[CR]]
		cat(CR, '\n')
	}
	names(snp_list)<-seq(1:22)
	saveRDS(snp_list, file=paste0(parentdir, 'output/hei_', args[2], '.Rds'))
}
cat('Got SNP IDs, CHR and POS.\n')
#
if((args[2] %in% c("LD_50000_0.01_0.5", "LD_100000_0.01_0.5", "LD_250000_0.01_0.5"))){
       hei<-vector('list', 22)
       for(Z in 1:22){
             	fread(paste0('zcat ', parentdir, 'output/hei_SNPs_ALL_chr', Z, '.vcf.gz'), header=T)-> hei[[Z]]
                colnames(hei[[Z]])<-unlist(fread('/project/mathilab/bbita/gwas_admix/new_height/JHS/header.txt', header=F, sep="\t"))
                colnames(hei[[Z]])[1]<-'CHR'
		colnames(hei[[Z]])[3]<-'MarkerName'
		unique(hei[[Z]], by=c('CHR','POS'))-> hei[[Z]]
                hei[[Z]][REF %in% c("A","C","G","T")]-> hei[[Z]]
                setkey(snp_list[[Z]],CHR, POS)
                setkey(hei[[Z]], CHR, POS)
                hei[[Z]][snp_list[[Z]], nomatch=0]-> hei[[Z]]
		hei[[Z]][,CHR:=NULL]
		cat(Z, ' done \n')
	}

        names(hei)<-seq(1:22)
        saveRDS(hei, file=paste0(parentdir, 'output/hei_', args[2], '_v2.Rds'))
	cat('another checkpoint\n')
} else{
        saveRDS(snp_list, file=paste0(parentdir, 'output/hei_', args[2], '_v2.Rds'))
}
cat('checkpoint 3\n')
PGS<-vector('list',22)
names(PGS)<-c(1:22)
for (CR in 1:22){
	print(paste0("Chromosome is ", CR))
        try(PolScore4(CHR=CR, panel='JHS', tag=args[2]))-> PGS[[CR]]
        cat(paste0(CR, '  done\n'))
}

samps<-names(PGS[[1]]) #sum PGS across chromosomes.
PGS2<-vector('list', length(samps))
names(PGS2)<-samps
for (S in samps){
	sum(PGS[[1]][[S]],PGS[[2]][[S]],PGS[[3]][[S]],PGS[[4]][[S]],PGS[[5]][[S]],PGS[[6]][[S]],PGS[[7]][[S]],PGS[[8]][[S]],PGS[[9]][[S]],PGS[[10]][[S]],PGS[[11]][[S]],PGS[[12]][[S]],PGS[[13]][[S]],PGS[[14]][[S]],PGS[[15]][[S]],PGS[[16]][[S]],PGS[[17]][[S]],PGS[[18]][[S]],PGS[[19]][[S]],PGS[[20]][[S]],PGS[[21]][[S]],PGS[[22]][[S]])->PGS2[[S]]
        cat(paste0(S, ' done\n'))
        }
saveRDS(PGS2, file=paste0(parentdir, 'output/PGS_JHS_', args[2], '.Rds'))
#TheEnd

