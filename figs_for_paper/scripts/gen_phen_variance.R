#!/usr/bin/env Rscript
############################
library(data.table)

##read in a pruned set of SNPs and retain CHR and POS and write it into a file
library(TeachingDemos)
txtStart(paste0("~/height_prediction/figs_for_paper/gen_phen_varianct.txt"))
dtsets<-vector("list", 6)
names(dtsets)<-c("WHI","JHS", "ukb_afr", "HRS_afr", "ukb_eur", "HRS_eur")

for(I in names(dtsets)){
	readRDS(paste0('~/height_prediction/gwas/', I, '/output/hei_phys_100000_0.0005_v2.Rds'))-> betas
	betas<-lapply(betas, function(X) X[,.(CHR,POS, MarkerName, REF, ALT, Allele1, Allele2, b, SE, p, N)])
	betas<-do.call(rbind, betas)
	colnames(betas)[3]<-"SNP"
	res2<-vector("list",22)
	for(CR in 22:1){
		betas[,.(CHR,POS)][order(CHR,POS)]-> fwr
		fwr[CHR==CR][order(CHR,POS)]-> tp
		fwrite(tp, file=paste0("~/height_prediction/tmp/tmp",CR,"_", I, ".txt"), sep="\t", col.names=F)
		if (I=="ukb_afr"){
			com<-paste0('bcftools view -H -T ~/height_prediction/tmp/tmp', CR, '_', I, '.txt ~/height_prediction/input/',I, '/chr', CR,'_bas_afr.vcf.gz > ~/height_prediction/tmp/res_chr',CR, "_", I, '.txt')	
		} else if (I=='ukb_eur'){
			com<-paste0('bcftools view -H -T ~/height_prediction/tmp/tmp', CR, '_', I, '.txt ~/height_prediction/input/',I, '/chr', CR,'_bas_eur.vcf.gz > ~/height_prediction/tmp/res_chr',CR, "_", I, '.txt') 
		} else {
			com<-paste0('bcftools view -H -T ~/height_prediction/tmp/tmp', CR, '_', I, '.txt ~/height_prediction/input/',I, '/chr', CR,'_bas.vcf.gz > ~/height_prediction/tmp/res_chr',CR, "_", I, '.txt')
		} 
		system(com);cat('sleep\n')
		Sys.sleep(30)
		system(paste0("cat ~/height_prediction/input/", I, "/header_", I, ".txt > ~/height_prediction/tmp/temp2_chr", CR,"_",I, ".txt"))
		system(paste0("cat ~/height_prediction/tmp/res_chr",CR,"_", I,".txt >> ~/height_prediction/tmp/temp2_chr", CR, "_",I,".txt"))
		comm<-paste0("plink --vcf ~/height_prediction/tmp/temp2_chr", CR,"_",I, ".txt --double-id --freq --out ~/height_prediction/tmp/out_chr", CR, "_", I) ## calculate allele frequency and then heterozygosity : p*q
		system(comm)
		Sys.sleep(10)
		res<-fread(paste0('~/height_prediction/tmp/out_chr',CR,"_",I, '.frq'), header=T)
		res[,pq:=MAF*(1-MAF)]
		res1<-fread(paste0('~/height_prediction/tmp/res_chr', CR,"_",I, ".txt" ))[,1:3]
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
	saveRDS(dtsets[[I]], file=paste0("~/height_prediction/output/tmp_dtsets_",I,".Rds"))
	system("rm ~/height_prediction/tmp/*") #clear everything
	cat(I, " done\n")
}

saveRDS(dtsets, file="~/height_prediction/output/dtsets.Rds")


lapply(dtsets, function(X) sum(X$value, na.rm=T))

fread('/project/mathilab/data/WHI/data/1kg/1kg_affy6_all_CEU.afreq')-> ceu
colnames(ceu)[1]<-'CHR'
colnames(ceu)[2]<-'i.MarkerName'
readRDS('~/height_prediction/output/tmp_dtsets_WHI.Rds')-> whi
setkey(whi, CHR,i.MarkerName)
setkey(ceu, CHR, i.MarkerName)
whi[ceu, nomatch=0]-> final
final[,pq:=ALT_FREQS* (1-ALT_FREQS)]
final[, value2:=(b^2)*pq]

readRDS('~/height_prediction/output/tmp_dtsets_JHS.Rds')-> jhs
setkey(jhs, CHR, i.MarkerName)
setkey(ceu, CHR, i.MarkerName)
jhs[ceu, nomatch=0]-> final2
final2[,pq:=ALT_FREQS* (1-ALT_FREQS)]
final2[, value2:=(b^2)*pq]

lapply(1:6, function(X) sum(dtsets[[X]]$value, na.rm=T))
unlist(lapply(1:2, function(X) sum(dtsets[[X]]$value, na.rm=T)))/c(sum(final$value2), sum(final2$value2))
txtStop()
