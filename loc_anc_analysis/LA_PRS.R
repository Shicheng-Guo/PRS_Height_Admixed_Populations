LA_PRS<- function(X=geno2, X2=1, Y=anc2, alpha=as.numeric(args[3])){
dataA<-cbind(data.table(Anc=unlist(Y[,X2]), plink_prun))[, Geno:=ifelse(AlM=='YES', (1-unlist(X[,X2])), unlist(X[,X2]))] ##if YES, 0=REF and 1=ALT. If NO, 1=REF, 0=ALT, thus 1-Geno
afrA<-dataA[Anc==1]
eurA<-dataA[Anc==2]
if(nrow(afrA)>=1){
afrA[,PRS_afr:=PLINK*Geno]
afrA[,PRS_eur:=b*Geno]
afrA[,PRS_part:=alpha*PRS_afr+(1-alpha)*PRS_eur]
}
if(nrow(eurA)>=1){
eurA[,PRS_afr:=0]
eurA[,PRS_eur:=b*Geno]
eurA[,PRS_part:=PRS_afr+PRS_eur]
}
res<-sum(bind_rows(afrA,eurA)$PRS_part, na.rm=T)
return(res)
}
