PolScore<- function(hei2=hei_q1){
        samps<-colnames(hei2)[9:(ncol(hei2)-8)]
        setDT(hei2)
        #hei2<-hei[,c(1:5,10:7301)]
        hei2[ALT==Allele1 & REF==Allele2]-> temp1
        hei2[REF==Allele1 & ALT==Allele2]-> temp2 #im ignoring the other two rows for now
        vector('list', length(samps))-> temp_list
        names(temp_list)<-samps
        if(nrow(temp1)>0 & nrow(temp2)>0){
               for(i  in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b]*2)-> temp_list[[i]]
                        temp_list[[i]] + sum(temp2[which(temp2[,i, with=F]=="0/0"),b]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b]*0)-> temp_list[[i]]
                }
        } else if (nrow(temp1)>0){
                for(i in samps){
                        sum(temp1[which(temp1[,i, with=F]=="0/0"),b]*0) + sum(temp1[which(temp1[,i,with=F]=="0/1"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/0"),b]*1) + sum(temp1[which(temp1[,i,with=F]=="1/1"),b]*2)-> temp_list[[i]]
                }
        } else if (nrow(temp2)>0){
                for(i in samps){
                        sum(temp2[which(temp2[,i, with=F]=="0/0"),b]*2) + sum(temp2[which(temp2[,i,with=F]=="0/1"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/0"),b]*1) + sum(temp2[which(temp2[,i,with=F]=="1/1"),b]*0)-> temp_list[[i]]
                }
        }
        return(temp_list)
}
