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
PolScore2<- function(hei2=hei_q1){
	#readRDS(paste0(home, panel, "/", panel2,  '/output/hei_', tag, '_v2.Rds'))-> hei
	samps<-colnames(hei2)[9:(ncol(hei2)-8)]
        setDT(hei2)
        #hei2<-hei[,c(1:5,10:7301)]
        hei2[ALT==Allele1 & REF==Allele2]-> temp1
        hei2[REF==Allele1 & ALT==Allele2]-> temp2 #im ignoring the other two rows for now
	if(panel=='sib_betas'){
        	samps<-colnames(hei2)[9:(ncol(hei2)-6)]
		} else if (panel=='gwas'){
	        samps<-colnames(hei2)[9:(ncol(hei2)-8)] 
	}
        hei2[ALT==Allele1]-> temp1
        hei2[REF==Allele1]-> temp2 #im ignoring the other two rows for now
        vector('list', length(samps))-> temp_list
        names(temp_list)<-samps
	cat('Number of samples is', length(samps), '\n')
        if(nrow(temp1)>0 & nrow(temp2)>0){
		matrix(nrow=nrow(temp1)+nrow(temp2), ncol=length(samps))-> my_matrix
        	colnames(my_matrix)<-samps
        	rownames(my_matrix)<-c(temp1[,MarkerName],temp2[,MarkerName])
        	b1<-temp1[,b]
       		b2<-temp2[,b]
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
		 apply(my_matrix*c(b1,b2), 2, sum, na.rm=T)-> res
		cat('Finished first for loop\n')
       		} else if (nrow(temp1)>0){
                matrix(nrow=nrow(temp1), ncol=length(samps))-> my_matrix
                colnames(my_matrix)<-samps
                rownames(my_matrix)<-temp1[,MarkerName]
                b1<-temp1[,b]
                counter<-0
		for(i in samps){
		        my_matrix[which(temp1[,i, with=F]=="0/0"),i]<-0
                        my_matrix[which(temp1[,i, with=F]=="1/1"),i]<-2
                        my_matrix[which(temp1[,i, with=F]=="1/0"),i]<-1
                        my_matrix[which(temp1[,i, with=F]=="0/1"),i]<-1
                }
		 apply(my_matrix*b1, 2, sum, na.rm=T)-> res
	cat('Finished second  for loop\n')
        } else if (nrow(temp2)>0){
	        matrix(nrow=nrow(temp2), ncol=length(samps))-> my_matrix
                colnames(my_matrix)<-samps
                rownames(my_matrix)<-temp2[,MarkerName]
                b2<-temp2[,b]
                for(i in samps){
        	
                        my_matrix[which(temp2[,i, with=F]=="0/0"),i]<-2
                        my_matrix[which(temp2[,i, with=F]=="1/1"),i]<-0
                        my_matrix[which(temp2[,i, with=F]=="1/0"),i]<-1
                        my_matrix[which(temp2[,i, with=F]=="0/1"),i]<-1        
	}
	 apply(my_matrix*b2, 2, sum, na.rm=T)-> res
        }
	cat('Finished third  for loop\n')
#acollect sample names for this population from 1000G data.
        return(res)
}


#*******
#* END *
#*******
