
short_fun_imputed<-function(args=22, data=prun2){
        prun2<-setDT(readRDS('~/height_prediction/gwas/WHI/output/hei_phys_100000_0.0005_v2.Rds')[[args]])
        final_plink[CHR==args]-> plink
        setkey(plink, CHR, POS)
        setkey(prun2, CHR,POS)
        prun2[plink, nomatch=0]-> prun2
        prun2[abs(T_STAT)>1]-> prun2_TSTAT_plink_1
        prun2[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        prun2_TSTAT_plink_1[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
	prun2_TSTAT_plink_1[,i.MarkerName.2:=NULL][,i.UNADJ:=NULL][,i.GC:=NULL][,i.BONF:=NULL][,i.HOLM:=NULL][,i.SIDAK_SS:=NULL][,i.SIDAK_SD:=NULL][,i.FDR_BH:=NULL][,i.FDR_BY:=NULL][,i.REF.1:=NULL][,i.ALT.1:=NULL][,i.TEST:=NULL][,i.OBS_CT:=NULL][,i.SE.1:=NULL][,i.T_STAT:=NULL][,i.A1:=NULL][,i.PLINK:=NULL]
	cat('start PRS\n')
        test4<-PolScore3_imputed(beta='PLINK',data=prun2) #plink
        test4b<-PolScore3_imputed(beta='PLINK',data=prun2_TSTAT_plink_1)
        return(list(test4, test4b))
}

short_fun_imputed_v2<-function(args=22, data=prun2){
        prun2<-setDT(readRDS('~/height_prediction/gwas/HRS_afr/output/hei_phys_100000_0.0005_v2.Rds')[[args]])
        final_plink[CHR==args]-> plink
        setkey(plink, CHR, POS)
        setkey(prun2, CHR,POS)
        prun2[plink, nomatch=0]-> prun2
        prun2[abs(T_STAT)>1]-> prun2_TSTAT_plink_1
	prun2[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        prun2_TSTAT_plink_1[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        prun2_TSTAT_plink_1[,i.MarkerName.2:=NULL][,i.UNADJ:=NULL][,i.GC:=NULL][,i.BONF:=NULL][,i.HOLM:=NULL][,i.SIDAK_SS:=NULL][,i.SIDAK_SD:=NULL][,i.FDR_BH:=NULL][,i.FDR_BY:=NULL][,i.REF.1:=NULL][,i.ALT.1:=NULL][,i.TEST:=NULL][,i.OBS_CT:=NULL][,i.SE.1:=NULL][,i.T_STAT:=NULL][,i.A1:=NULL][,i.PLINK:=NULL]
	cat('start PRS\n')
        test4<-PolScore3_imputed(beta='PLINK',data=prun2) #plink
        test4b<-PolScore3_imputed(beta='PLINK',data=prun2_TSTAT_plink_1)
        return(list(test4, test4b))
}


short_fun_imputed_v3<-function(args=22, data=prun2){
        prun2<-setDT(readRDS('~/height_prediction/gwas/JHS/output/hei_phys_100000_0.0005_v2.Rds')[[args]])
        final_plink[CHR==args]-> plink
        setkey(plink, CHR, POS)
        setkey(prun2, CHR,POS)
        prun2[plink, nomatch=0]-> prun2
	prun2[abs(T_STAT)>1]-> prun2_TSTAT_plink_1
        prun2[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        prun2_TSTAT_plink_1[,i.MarkerName.1:=NULL][,UNADJ:=NULL][,GC:=NULL][,BONF:=NULL][,HOLM:=NULL][,SIDAK_SS:=NULL][,SIDAK_SD:=NULL][,FDR_BH:=NULL][,FDR_BY:=NULL][,i.REF:=NULL][,i.ALT:=NULL][,TEST:=NULL][,OBS_CT:=NULL][,i.SE:=NULL][,T_STAT:=NULL][,i.MarkerName:=NULL]
        prun2_TSTAT_plink_1[,i.MarkerName.2:=NULL][,i.UNADJ:=NULL][,i.GC:=NULL][,i.BONF:=NULL][,i.HOLM:=NULL][,i.SIDAK_SS:=NULL][,i.SIDAK_SD:=NULL][,i.FDR_BH:=NULL][,i.FDR_BY:=NULL][,i.REF.1:=NULL][,i.ALT.1:=NULL][,i.TEST:=NULL][,i.OBS_CT:=NULL][,i.SE.1:=NULL][,i.T_STAT:=NULL][,i.A1:=NULL][,i.PLINK:=NULL]
	cat('start PRS\n')
	test4<-PolScore3_imputed(beta='PLINK',data=prun2) #plink
        test4b<-PolScore3_imputed(beta='PLINK',data=prun2_TSTAT_plink_1)
        return(list(test4, test4b))
}
