#!/usr/bin/env Rscript
############################
############################
#Load packages
suppressMessages(library(getopt)) # Package designed to be used with Rscript to write #! shebang scripts that accept short and long ﬂags/options.
suppressMessages(library("optparse"))
options(scipen=10)
options(digits=10)
#

currentdir="~/height_prediction/scripts/"
parentdir=gsub("scripts/","", currentdir)

TIME.start <- Sys.time() #record time spent running this script.

spec<-matrix(c('chromosome', 'chr', 1, "numeric",
               'windowsize', 'w', 2, "numeric",
               'pvalue', 'p', 2,  "numeric",
               'method', 'm',1, "character",
               'r2', 'r', 2, "numeric",
		'dataset1', 'dataset1',1, "character",
		'dataset2', 'dataset2',1, "character",
               'help', 'h', 0, "logical"
               ),ncol=4,byrow=T)

getopt(spec)-> opt

HELP.MESSAGE <- paste(c(paste0("The script <<",self = commandArgs()[4],">> requires the arguments:"),
		"[-chromosome | -i]\t\t\t\t\t<chromosome to be analyzed>", 
		"[-windowsize | -w][OPTIONAL]\t\t\t\t<number of bp to be used for physical distance or LD clumping (if method is genet, it is the number of centimorgans (cM)>", 
		"[-pvalue | -p][OPTIONAL]\t\t\t\t<p value threshold for top SNPs>", "[method | -m]\t\t\t\t\t\t<Method used for clumping: physical distance (phys), genetic distance (genet), r2-based (LD)>", 
		"[r2  | -r][OPTIONAL]\t\t\t\t\t<If method is LD, this option will tell what r2 threshold to use for clumping>"),collapse="\n")

## If help was asked for.
if ( !is.null(opt$help) ) {
      ##get the script name (only works when invoked with Rscript).
	cat(HELP.MESSAGE)
	cat('\n')
	q(status=1); # quit
}

a <- 0

if ( is.null(opt$chromosome)) { cat("flag [-chromosome| -chr] is required\n")} else {a=a+1}
if ( is.null(opt$method)) { cat("flag [-method | -b]> is required for performing some kind of clumping\n")} else {a=a+1}

if(a < 2) {
	cat(paste0("ERROR: not all necessary arguments were specified.\nFor help type:\n $",sub("--file=","",commandArgs()[4])," -h\n"))
	q(status=1); # quit if not all arguments are specified.
}

cat(paste0('Chromosome is ', opt$chromosome, " and clumping method is ", opt$method, '\n'))

if (is.null(opt$windowsize)){
	cat("  <windowsize> was not specified [-w] so no distance based clumping will be implemented\n")
} else {
	cat('You provided a window size. Distance based clumping will occur if method is phys or genet .\n')
	W <- opt$w
	if(opt$method == "phys"){
		cat(paste0('Window size is ', W, " base pairs\n"))
	}
	if(opt$method == 'LD'){
		cat(paste0('Window size is ', W, " kbp\n"))
	}
	else if (opt$method ==  "genet"){
		cat(paste0("Window size is ", W, " cM\n"))
	}
}
if (is.null(opt$pvalue)){
	cat("  <pvalue> was not specified [-p] so no pvalue based clumping will be implemented\n")
} else {
        cat('You provided a p value. This will be used for clumping.\n')
}

if (is.null(opt$r2)){
	cat("  <r2> was not specified [-r] so no r2 based clumping will be implemented\n")
} else {
	cat('You provided an r2 value. This will be used for clumping IF method is LD.\n')
}

#****************
#Load packages **
#*****************
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(biomaRt))
suppressMessages(library(rlist))

source('~/height_prediction/scripts/LD_prun.R') #clumping function.

cat ('Beginning clumping process\n')

if(opt$method == "LD"){
	try(LD_prun(dis=opt$method, W=opt$windowsize, p_thresh=opt$pvalue, CHR=opt$chromosome, r2=opt$r2, dataset1=opt$dataset1, dataset2=opt$dataset2))->res
	
} else {
	try(LD_prun(dis=opt$method, W=opt$windowsize, p_thresh=opt$pvalue, CHR=opt$chromosome, dataset1=opt$dataset1, dataset2=opt$dataset2))->res
}

cat('Another checkpoint\n')

if(opt$method=="LD"){
#	system('cd ..')
#	setwd('/project/mathilab/bbita/gwas_admix/sib_gwas/JHS/')
	saveRDS(res,file=paste0(parentdir, opt$dataset1, "/", opt$dataset2, "/prunned_1kg/LD_prunned_hei_chr", opt$chromosome, '_', opt$method, '_', opt$windowsize, '_', opt$pvalue, '_', opt$r2, '.Rds'))

} else if (opt$method=="LD_block"){
	saveRDS(res[['AFR']], file=paste0(parentdir,opt$dataset1, "/", opt$dataset2, "/prunned_1kg/LD_prunned_hei_chr", opt$chromosome, '_', opt$method, '_', opt$windowsize, '_', opt$pvalue, '_AFR.Rds'))
	saveRDS(res[['EUR']], file=paste0(parentdir, opt$dataset1, "/", opt$dataset2, "/prunned_1kg/LD_prunned_hei_chr", opt$chromosome, '_', opt$method, '_', opt$windowsize, '_', opt$pvalue, '_EUR.Rds'))
} else {
	saveRDS(res,file=paste0(parentdir, opt$dataset1, "/", opt$dataset2, "/prunned_1kg/LD_prunned_hei_chr", opt$chromosome, '_', opt$method, '_', opt$windowsize, '_', opt$pvalue, '.Rds'))
}

print(paste("Elapsed time is", round(as.numeric(difftime(Sys.time(), TIME.start,units="mins")), 2), "minutes"))

#******
# END *
#******
