#!/usr/bin/env Rscript
#######################################################################################################################################################
#Script Name    :gnomad_freq.R
#Description    :
#Args           :a JSON output from a query from gnomac broser
#Author         :Bárbara D Bitarello
#Email          :barbarabitarello@gmail.com
#Date created   :22-07-2020
#Requirements   :
#USage          :
########################################################################################################################################################
#load packages
suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
	library(rjson)
        library(optparse)
	library(httr)
	library(jsonlite)
	
})
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
#save start time
ptm <- proc.time()
######
js<-args[1]
js<-'results-all.json'
#
cat('Read in JSON file...\n')
json_data <-fromJSON(paste0('input/', js))
setDT(json_data)
js_tidy<-json_data %>% group_by(id) %>% summarise(REF=refAllele, ALT=altAllele, Freq_Alt=altAlleleFreq, Study=study) %>% rename(ID=id) %>% as.data.table
cat('Done\n')












######
tim<-proc.time()[[3]] - ptm[[3]]
cat('Total elapsed time is ', tim, ' seconds\n')
#End

