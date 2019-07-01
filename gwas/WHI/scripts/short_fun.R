
short_fun<-function(args=22, data=prun2){
        what <- paste0("~/height_prediction/input/ukb_afr/UKB_kgCY_chr", args)
        betas<-fread(paste0('~/height_prediction/gwas/ukb_afr/output/AS_Beta_chr', args,'example.txt'))
        prun<-setDT(readRDS('~/height_prediction/gwas/WHI/output/hei_phys_100000_0.0005_v2.Rds')[[args]])
        fread(paste0(what, '.phsnp'))-> snps
        colnames(snps)<-c("ID", "CHR", "V3", "POS","REF","ALT")
        tag='phys_100000_0.0005'
        cbind(betas, snps)-> betas2
        data<-betas2[which(betas2[, POS] %in% prun[, POS]),]
        prun[POS %in% betas2[, POS]]-> prun2
        prun[POS %in% betas2[abs(Tstat_all)>1][, POS]]-> prun_tstat_all_1
        #prun[POS %in% betas2[abs(Tstat_all)>2][, POS]]-> prun_tstat_all_2
        #data$POS==prun2$POS
        prun2[, POP1:=data$POP1]
        prun2[, POP2:=data$POP2]
        prun2[, ALL:=data$ALL]
        prun_tstat_all_1[, POP1:=data$POP1]
	        prun_tstat_all_1[, POP2:=data$POP2]
        prun2_tstat_all_1[, ALL:=data$ALL]

        final_plink[CHR==args]-> plink
        setkey(plink, CHR, POS)
        setkey(prun2, CHR,POS)
        setkey(prun_tstat_all_1, CHR,POS)
        prun2[plink, nomatch=0]-> prun2
        prun_tstat_all_1[plink, nomatch=0]-> prun_tstat_all_1
        prun2[plink, nomatch=0][abs(T_STAT)>1]-> prun2_TSTAT_plink_1
        #prun2[plink, nomatch=0][abs(T_STAT)>2]-> prun2_TSTAT_plink_2
        prun2[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        prun2_TSTAT_plink_1[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        prun_tstat_all_1[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        test<-PolScore3(beta="POP1",data=prun2)
        test2<-PolScore3(beta="POP2",data=prun2)
        test3<-PolScore3(beta='ALL',data=prun2)
        test3b<-PolScore3(beta='ALL',data=prun_tstat_all_1)
       test4<-PolScore3(beta='PLINK',data=prun2) #plink
        test4b<-PolScore3(beta='PLINK',data=prun2_TSTAT_plink_1)
        return(list(test,test2,test3,test3b,test4, test4b))
}

