#!/usr/bin/env Rscript
############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}

library(data.table)
#args<-'HRS_eur'
#pruned<-readRDS(paste0('~/height_prediction/gwas/', args[1], '/output/hei_phys_100000_0.0005_v2.Rds'))

dt<-vector('list',22)
names(dt)<-seq(1:22)
for (chr in 22:1){
	pruned<-readRDS(paste0('~/height_prediction/gwas/', args[1], '/output/hei_chr', chr, '.Rds'))
	temp1<-pruned[ALT==Allele1]
        temp2<-pruned[REF==Allele1]
	samps<-colnames(pruned)[grep("[0-9]+_",colnames(pruned))]

	        if(nrow(temp1)>0 & nrow(temp2)>0){
                matrix(nrow=nrow(temp1)+nrow(temp2), ncol=length(samps))-> my_matrix
                colnames(my_matrix)<-samps
                rownames(my_matrix)<-c(temp1[,MarkerName],temp2[,MarkerName])
                counter<-0
                for(i  in samps){
                        my_matrix[which(temp1[,i, with=F]=="0/0"),i]<-0
                        my_matrix[which(temp1[,i, with=F]=="1/1"),i]<-2
                        my_matrix[which(temp1[,i, with=F]=="1/0"),i]<-1
                        my_matrix[which(temp1[,i, with=F]=="0/1"),i]<-1

                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="0/0"),i]<-2
                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="1/1"),i]<-0
                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="1/0"),i]<-1
                        my_matrix[nrow(temp1)+which(temp2[,i, with=F]=="0/1"),i]<-1
                        counter<-counter+1
                        cat(counter,'\r')
                }
        }	        
	lapply(1:nrow(my_matrix), function(X) table(my_matrix[X,]))-> res
	freq<-lapply(1:length(res), function(X) sum(res[[X]]['2'],(res[[X]]['1']/2), na.rm=T)/length(samps))
	dt[[chr]]<-data.table(MarkerName=c(temp1[,MarkerName], temp2[,MarkerName]), Freq_effect_allele=freq, CHR=chr, POS=c(temp1[,POS], temp2[,POS]))
	cat(chr, ' done\n')
	remove(pruned)
	gc()
}

dt2<-do.call(rbind, dt)
saveRDS(dt2,paste0('~/height_prediction/epistasis/output/table_', args[1], '.Rds'))

