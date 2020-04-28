library(data.table)
library(dplyr)
for(chr in 22:1){
fread(paste0('output/1000g_chr', chr, '.bim'))-> bim
bim2<-unique(bim, by=c("V1", "V2", "V4")) %>% as.data.table
fwrite(bim2, file=paste0('output/1000g_chr', chr, '.bim'),  col.names=F)
cat(chr, '\n')
}

