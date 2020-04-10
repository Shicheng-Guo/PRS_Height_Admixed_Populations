suppressMessages(library(data.table))
suppressMessages(library(dplyr))

currentdir<-paste0(getwd(),"/")
parentdir<-gsub("scripts/", "", parentdir)

LD_prun<-function(W=NULL,p_thresh=NULL, dis=opt$method, CHR=opt$chromosome, r2=NULL, dataset1=opt$dataset1, dataset2=opt$dataset2){
	as.numeric(W)-> W;as.numeric(p_thresh)-> p_thresh
	readRDS(paste0(parentdir, opt$dataset1, "/", opt$dataset2, "/output/hei_chr", CHR, '.Rds'))-> hei
	cat('checkpoint 1\n')
	hei[order(p)]-> tp
	cat('checkpoint 2\n')
	remove(hei)
	cat('Method is ', dis, "\n")
	cat('Chromosome is ', CHR, "\n")
	if (dis =='phys'){   
		tp[p<=p_thresh]-> tp
		tp %>% dplyr::arrange(p) %>% as.data.table-> tp
		tp[1,POS]->p1
		tp[, Dist:=abs(POS-p1)]
	} else if (dis == 'genet'){
		tp[p<=p_thresh]-> tp
		fread(paste0('zcat /project/mathilab/data/maps/hm2/hm2/genetic_map_GRCh37_chr', CHR,'.txt.gz'))[,CHR:=gsub("chr","",Chromosome)][, Chromosome:=NULL]-> rec
		colnames(rec)<-c('POS', 'RATE_cM_Mb', 'MAP_cM', 'CHR')
		rec[, MAP_bp:=POS-POS[1]]
		rec$MAP_bp->x
		rec$MAP_cM->y
		f <- approxfun(x, y, rule = 2:1) #interpolate physical and genetic distance.
		tp %>% dplyr::arrange(p) %>% as.data.table-> tp
		tp[, Dist:=abs(POS-tp[1,POS])]
		tp[,Dist_cM:=f(Dist)]
	} else if (dis == 'LD'){
		system('cd plink2')
		setwd(paste0(parentdir, opt$dataset1, "/", opt$dataset2, "/plink2/output/"))
#run PLINK with these parameters #TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO
		if(dataset1=='sib_betas'){
		a.str<- paste(paste0("plink \\", "--noweb \\", "--bfile chr", CHR, " \\",  "--clump betas_for_plink.txt \\", paste0("--clump-p1 ", p_thresh, " \\"), paste0("--clump-r2 ", r2, " \\"),  paste0("--clump-kb ", W/1000," \\"), paste0("--out giant_out_clump_p1_", p_thresh, "_r2_", r2, "_kb_", W/1000, "_chr_", CHR), sep="\n"))
		} else if (dataset1=='gwas'){
		a.str<- paste(paste0("plink \\", "--noweb \\", "--bfile chr", CHR, " \\",  "--clump betas_for_plink.txt \\", paste0("--clump-p1 ", p_thresh, " \\"), paste0("--clump-r2 ", r2, " \\"),  paste0("--clump-kb ", W/1000," \\"), paste0("--out giant_out_clump_p1_", p_thresh, "_r2_", r2, "_kb_", W/1000, "_chr_", CHR), sep="\n"))
		}
		system(a.str)
		Sys.sleep(60)
		fread(paste0("giant_out_clump_p1_", p_thresh, "_r2_", r2, "_kb_", W/1000, "_chr_", CHR, ".clumped"))->a
		vec<-a$SP2
		vec2<-a$SNP
	} else if (dis == "LD_block"){ #need to test
		fread(paste0('/project/mathilab/bbita/gwas_admix/height/LD_blocks/nygcresearch-ldetect-data-ac125e47bf7f/AFR/fourier_ls-chr', CHR, '.bed'))-> blck #using AFR blocks
		fread(paste0('/project/mathilab/bbita/gwas_admix/height/LD_blocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-chr', CHR, '.bed'))-> blck_eu #using EUR blocks
		colnames(blck)<-c("CHR", "POS1","POS2")
		colnames(blck_eu)<-c("CHR", "POS1","POS2")
		blck[, CHR:=as.integer(gsub("chr", "", CHR))]
		blck_eu[, CHR:=as.integer(gsub("chr", "", CHR))]
		blck[, ID:=paste0(CHR, "|", POS1, "|", POS2)]
		blck_eu[, ID:=paste0(CHR, "|", POS1, "|", POS2)]
		blck[, POS1:=as.integer(POS1)] #this isbecause chr11 has a problem one of the pos1 is absent (NA)
		blck[, POS2:=as.integer(POS2)]
		blck_eu[, POS1:=as.integer(POS1)] #this isbecause chr11 has a problem one of the pos1 is absent (NA)
		blck_eu[, POS2:=as.integer(POS2)]
		na.omit(blck)-> blck #because of misisng coords
		na.omit(blck_eu)-> blck_eu
		tp[, POS1:=POS]
		tp[, POS2:=POS]
		setkey(tp, CHR, POS1, POS2)
		setkey(blck, CHR, POS1, POS2)
		setkey(blck_eu, CHR, POS1, POS2)
		foverlaps(tp, blck, type='within')-> tmp
		tmp %>% dplyr::group_by(ID) %>% dplyr::summarise(top_SNP=min(p))-> tmp2
		foverlaps(tp, blck_eu, type='within')-> tmp_eu
		tmp_eu %>% dplyr::group_by(ID) %>% dplyr::summarise(top_SNP=min(p))-> tmp2_eu
		vec<-c();vec2<-c();
		vec_eu<-c(); vec2_eu<-c();
		for (i in 1: nrow(tmp2)){
			append(vec2, tmp[ID== tmp2$ID[i] & p==tmp2$top_SNP[i]][,MarkerName])-> vec2
			cat('block ')
			cat(i)
			cat(' done\n')
			}
		
		for (i in 1: nrow(tmp2_eu)){
			append(vec2_eu, tmp_eu[ID== tmp2_eu$ID[i] & p==tmp2_eu$top_SNP[i]][,MarkerName])-> vec2_eu
			cat('block ')
			cat(i)
			cat('done\n')
			}
		
		}
		cat('checkpoint 3\n')
		if(dis == 'phys'){
			vec<-c();vec2<-c();
			while(nrow(tp[p<=p_thresh])>1){
       				cat('Starting another round\n')
				tp %>% dplyr::arrange(p) %>% as.data.table-> tp
       		 		ind<-tp[1,]
        			append(vec2, ind$MarkerName)-> vec2
        			cat(paste0('Index SNP is ', ind$MarkerName, "\n"))
        			tp %>% dplyr::arrange(p) %>% as.data.table-> tp
        			tp[, Dist:=abs(POS-ind[,POS])]
				#physcial window size
        			append(vec, tp[Dist<=abs(W)][, MarkerName])-> vec #assign these SNPs to a list of LD with the target SNP
        			cat(paste0('Physical window of ', W, 'kb around ', ind$MarkerName, ' removed\n'))
        			tp[Dist>abs(W)]-> tp
        			cat(paste0(nrow(tp), ' SNPs left to check\n'))
			}
		} else if (dis == 'genet'){
			vec<-c();vec2<-c();
			while(nrow(tp[p<=p_thresh])>1){
       				print('Starting another round')
				tp %>% dplyr::arrange(p) %>% as.data.table-> tp
        			ind<-tp[1,]
        			append(vec2, ind$MarkerName)-> vec2
        			print(paste0('Index SNP is ', ind$MarkerName))
        			tp %>% dplyr::arrange(p) %>% as.data.table-> tp
				tp[, Dist:=abs(POS-ind$POS)]
				tp[, Dist_cM:=f(Dist)]
				#genetic window size
        			append(vec, tp[Dist_cM<=W][, MarkerName])-> vec #assign these SNPs to a list of LD with the target SNP
        			print(paste0('Genetic window of ', W, 'cM around ', ind$MarkerName, ' removed'))
        			tp[Dist_cM>abs(W)]-> tp
        			print(paste0(nrow(tp), ' SNPs left to check'))
        		}
		}
		if (dis == "LD_block"){
		return(list(AFR=list(keep=vec2,rem=vec), EUR=list(keep=vec2_eu, rem=vec_eu)))
		} else {
		return(list(keep=vec2, rem=vec))
		}
}
#******
#*END *
#******
