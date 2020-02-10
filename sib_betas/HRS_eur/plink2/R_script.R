library(data.table)
for(I in 1:11){
readRDS(paste0('../output/hei_chr', I, '.Rds'))-> A
setDT(A)
A[,MarkerName:=i.MarkerName][, i.MarkerName:=NULL][,Allele1:=NULL][, Allele2:=NULL][, b:=NULL][, SE:=NULL][, p:=NULL][, N:=NULL]
colnames(A)[1]<-"#CHROM"
colnames(A)[3]<-"ID"
unique(A, by="ID")-> A
fwrite(A, file=paste0("output/chr", I, ".vcf"), sep="\t")
a<-paste0("bgzip output/chr", I, ".vcf")
system(a)
cat(I,"\n")
}
