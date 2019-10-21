#!/usr/bin/env Rscript

############################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
old <- Sys.time()
chr<-args[1]
plink_prun<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/plink_prun_chr', chr, '.Rds'))
plink_prun2<-readRDS(paste0('~/height_prediction/loc_anc_analysis/output/plink_prun2_chr', chr, '.Rds'))
input_MATRIX<-readRDS(file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr, '_temp.Rds'))
#write a function to calculate the ancestry-specific PRS
input<-do.call(cbind, input_MATRIX)
input2 <- sapply(input, as.numeric)
input2<-as.data.table(as.data.frame(input2))
LA_PRS<- function(X=plink_prun, X2="0_710880_A"){
        dataA<-cbind(select(X, CHR, MarkerName, i.MarkerName, POS, REF, ALT, Effect_Allele_plink, Allele1, Allele2, b, p,PLINK),  Set=unlist(input2[,get(X2)]))
        anc<-select(ancestry, contains(X2))
        afrA<-which(anc[,1]==1)
        eurA<-which(anc[,1]==2)
        dataAF<-dataA[afrA,]
        dataAF[,Set:=as.numeric(dataAF[,Set])]
        dataAF[, PRS_part:=ifelse(Effect_Allele_plink==ALT, Set*PLINK, ifelse(Allele1==REF, abs(1-Set)*PLINK, NA))]
        dataEU<-dataA[eurA,]
        dataEU[,Set:=as.numeric(dataEU[,Set])]
        dataEU[, PRS_part:=ifelse(Allele1==ALT, Set*b, ifelse(Allele1==REF, abs(1-Set)*b, NA))]
        res<-sum(rbind(dataAF,dataEU)$PRS_part, na.rm=T)
        return(res)
}
LA_PRS_v2<- function(X=plink_prun, X2="0_710880_A", alpha=as.numeric(args[3])){
        dataA<-cbind(select(X, CHR,POS, REF, ALT, Effect_Allele_plink, Allele1, b,PLINK))[ ,Set:=as.numeric(unlist(input2[,get(X2)]))]
        anc<-select(ancestry, X2)
        dataAF<-dataA[which(anc[,1]==1)][,AlleleMatch:=ifelse(Effect_Allele_plink==ALT, "TRUE", ifelse(Allele1==REF, "FALSE", NA))]
        na.omit(dataAF)-> dataAF
        dataEU<-dataA[which(anc[,1]==2)][,AlleleMatch:=ifelse(Allele1==ALT, "TRUE", ifelse(Allele1==REF, "FALSE", NA))]
        na.omit(dataEU)-> dataEU
        dataAF[, PRS_part:=ifelse(AlleleMatch=='TRUE', (1-alpha)*(Set*b)+(alpha*Set*PLINK), alpha*abs(1-Set)*PLINK+ (1-alpha)*abs(1-Set)*b)]
        dataEU[, PRS_part:=ifelse(AlleleMatch=='TRUE', Set*b, abs(1-Set)*b)]
        res<-sum(rbind(dataAF,dataEU)$PRS_part, na.rm=T)
        return(res)
}

system.time(LA_PRS()) #test
system.time(LA_PRS_v2())
colnames(input2)-> I
system.time(mclapply2(1:length(I), function(Z) LA_PRS(X2=I[Z]))-> test) #
saveRDS(test, file=paste0('~/height_prediction/loc_anc_analysis/output/chr', chr,'_',args[3], '_prs.Rds'))
