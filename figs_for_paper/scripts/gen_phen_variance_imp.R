#!/usr/bin/env Rscript
############################
library(data.table)
library(dplyr)
##read in a pruned set of SNPs and retain CHR and POS and write it into a file
dtsets<-vector("list", 4)
names(dtsets)<-c("HRS_afr", "HRS_eur", "UKB_afr", "UKB_eur")
G<-'100000_0.0005'
for(I in names(dtsets)){
	if(I %in% c('HRS_afr', 'HRS_eur')){
	do.call(rbind, readRDS(paste0('~/height_prediction/imputed/output/vec_all_100000_0.0005.Rds')))-> betas
	} else{	
	do.call(rbind, readRDS(paste0('~/height_prediction/imputed/output/UKB_vec_all_100000_0.0005.Rds')))-> betas
	}
	colnames(betas)[5]<-"SNP"
	res2<-vector("list",22)
	for(CR in 22:1){
		betas[,.(CHR,POS)][order(CHR,POS)]-> fwr
		fwr[CHR==CR][order(CHR,POS)]-> tp
		fwrite(tp, file=paste0("~/height_prediction/imputed/tmp/tmp",CR,"_", I, ".txt"), sep="\t", col.names=F)
		system(paste0('cp ~/height_prediction/imputed/header_', I, '.txt ~/height_prediction/imputed/tmp/res_chr',CR, "_", I, '.vcf'))
		com<-paste0('bcftools view -H -T  ~/height_prediction/imputed/tmp/tmp', CR, '_', I, '.txt ~/height_prediction/imputed/output/',I, '_', G, '_chr_', CR,'.vcf.gz >> ~/height_prediction/imputed/tmp/res_chr',CR, "_", I, '.vcf')	
		system(com);cat('sleep\n')
		Sys.sleep(10)
		comm<-paste0("plink --vcf ~/height_prediction/imputed/tmp/res_chr",CR, "_", I, ".vcf --double-id --freq --out ~/height_prediction/imputed/tmp/out_chr", CR, "_", I)
		system(comm)
		res<-fread(paste0('~/height_prediction/imputed/tmp/out_chr',CR,"_",I, '.frq'), header=T)
		res[,pq:=MAF*(1-MAF)]
		res1<-fread(paste0('~/height_prediction/imputed/tmp/res_chr', CR,"_",I, ".vcf" ))[,1:3]
		colnames(res1)<-c("CHR", "POS", "SNP")
		setkey(res1, CHR, SNP);setkey(res,CHR, SNP)
		res1[res]-> res2[[CR]]
		cat(I, " ", CR, " done\n")
		}
	do.call(rbind, res2)-> res3
	setkey(betas, CHR, POS, SNP)
	setkey(res3, CHR, POS, SNP)
	betas[res3]-> dtsets[[I]]
	dtsets[[I]][, value:=(b^2)*pq] 
	saveRDS(dtsets[[I]], file=paste0("~/height_prediction/imputed/output/tmp_dtsets_",I,".Rds"))
	system("rm ~/height_prediction/imputed/tmp/*") #clear everything
	cat(I, " done\n")
}
lapply(1:4, function(X) dtsets[[X]][, value2:=2*value])
saveRDS(dtsets, file="~/height_prediction/imputed/output/dtsets.Rds")

unlist(lapply(dtsets, function(X) sum(X$value, na.rm=T)))
# HRS_afr  HRS_eur  UKB_afr  UKB_eur
#16.86057 19.54020 17.21717 21.41683
