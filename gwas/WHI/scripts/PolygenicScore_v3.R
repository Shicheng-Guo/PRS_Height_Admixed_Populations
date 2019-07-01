PolScore3<- function(beta='ALL', data=prun2){
        samps<-colnames(data)[9:(ncol(data)-12)]
        setDT(data)
        #hei2<-hei[,c(1:5,10:7301)]
        cat('checkpoint \n')
	data[ALT==A1]-> temp1
        data[REF==A1]-> temp2 #im ignoring the other two rows for now
        if(nrow(temp1)>0 & nrow(temp2)>0){
	 	matrix(nrow=nrow(temp1)+nrow(temp2), ncol=length(samps))-> my_matrix
                colnames(my_matrix)<-samps
                rownames(my_matrix)<-c(temp1[,MarkerName],temp2[,MarkerName])
                b1<-temp1[,get(beta)]
                b2<-temp2[,get(beta)]
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
               	 	b1<-temp1[,get(beta)]
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
		        b2<-temp2[,get(beta)]
                	for(i in samps){
	                        my_matrix[which(temp2[,i, with=F]=="0/0"),i]<-2
        	                my_matrix[which(temp2[,i, with=F]=="1/1"),i]<-0
                	        my_matrix[which(temp2[,i, with=F]=="1/0"),i]<-1
                        	my_matrix[which(temp2[,i, with=F]=="0/1"),i]<-1
       			 }
         		apply(my_matrix*b2, 2, sum, na.rm=T)-> res
       		 }
        	cat('Finished third  for loop\n')
        return(res)
}

