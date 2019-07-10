library(data.table)
library(rbgen)
library(dplyr)

#read in summ stats file
ukb<-fread('zcat ~/height_prediction/gwas/input/50.assoc.tsv.gz', fill=T)
ukb[,c("CHR", "POS","Allele2","Allele1") := tstrsplit(variant, ":", fixed=TRUE)][,variant:=NULL]  #fix columns. In the UKB, variants are listed as CHR:POS: REF:ALT, where ALT is the effect allele. So, in order to be compatilble with my scripts (where allele 1 is effect all
ukb[, MarkerName:=rsid][, N:=nCompleteSamples][, AC:=NULL][, b:=beta][,p:=pval]
ukb[,rsid:=NULL][,nCompleteSamples:=NULL][, beta:=NULL][, pval:=NULL][, SE:=se]
ukb[,.(MarkerName,Allele1,Allele2, b, SE, p, N, CHR, POS)]-> ukb

#prune
p_thresh<-0.0005
W<-100000
vec_all<-vector('list',22)
for(chr in 22:1){
#chr<-21
	tmp<-ukb[CHR==chr]
	tmp$POS<-as.numeric(tmp$POS)
	tmp[order(p)]->tmp
	tmp[p<=p_thresh]-> tmp
	tmp[order(p)]-> tmp
	tmp[1,POS]->p1
	tmp[, Dist:=abs(POS-p1)]
	vec2<-c();
	while(nrow(tmp[p<=p_thresh])>1){
		cat('Starting another round\n')
       		tmp[order(p)]-> tmp
        	ind<-tmp[1,]
        	append(vec2, ind$MarkerName)-> vec2
        	cat(paste0('Index SNP is ', ind$MarkerName, "\n"))
        	tmp[order(p)]-> mp
        	tmp[, Dist:=abs(POS-ind[,POS])]
        #physcial window size
        	cat(paste0('Physical window of ', W, 'kb around ', ind$MarkerName, ' removed\n'))
        	tmp[Dist>abs(W)]-> tmp
        	cat(paste0(nrow(tmp), ' SNPs left to check\n'))
        	}
	cat('Chr ', chr, ' done\n')
	ukb[CHR!=chr]-> ukb
	fwrite(as.data.frame(vec2),paste0('~/height_prediction/imputed/output/chr_', chr, '.txt'), sep="\t",col.names=F)
	vec_all[[chr]]<-vec2
	remove(vec2)
}
saveRDS(vec_all,file='~/height_prediction/imputed/output/vec_all.Rds')

#
for(chr in 22:1){
	system(paste0('plink2 --bgen /project/mathilab/data/UKB/imputed/ukb_imp_chr',chr,'_afr.bgen --extract ~/height_prediction/imputed/output/chr_', chr,'.txt --recode vcf --out ~/height_prediction/imputed/output/chr_', chr))
	system(paste0('bgzip ~/height_prediction/imputed/output/chr_', chr, '.vcf'))
	cat(chr, 'done\n')
}




#path='/project/mathilab/data/UKB/imputed/'




#data = bgen.load(paste0(path,'ukb_imp_chr21_afr.bgen'), rsids=vec_all[[21]] )
#do pruning phys 100000 0.0005

#save

