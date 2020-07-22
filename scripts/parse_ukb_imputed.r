 library(data.table)
 library(dplyr)
 for (chr in 21:1) {
ukb <- fread(paste0("/project/mathilab/data/UKB/imputed/ukb_imp_chr", chr, "_eur.bim"))
+ hrs_eur <- fread("output/HRS_EUR.ldpred.bim")[V1 == chr]
+ setkey(hrs_eur, V1, V4)
+ setkey(ukb, V1, V4)
+ sincr <- hrs_eur[ukb, nomatch = 0]
+ assign(paste0("keep_chr_", chr), sincr)
+ cat(chr, "done\n")
+ }
KEEP <- rbind(keep_chr_1, keep_chr_2, keep_chr_3, keep_chr_4, 
keep_chr_5, keep_chr_6, keep_chr_7, keep_chr_8, keep_chr_9, 
keep_chr_10, keep_chr_11, keep_chr_12, keep_chr_13, keep_chr_14, 
keep_chr_15, keep_chr_16, keep_chr_17, keep_chr_18, keep_chr_19, 
keep_chr_20, keep_chr_21, keep_chr_22)
fwrite(as.data.frame(KEEP[, i.V2]), file = "~/height_prediction/ldpred/output/keep_hrs_eur_RS.txt", quote = F, col.names = F)
remove(list = ls())

for (chr in 22:1) {
ukb <- fread(paste0("/project/mathilab/data/UKB/imputed/ukb_imp_chr", chr, "_eur.bim"))
hrs_afr <- fread("output/HRS_AFR.ldpred.bim")[V1 == chr]
setkey(hrs_afr, V1, V4)
setkey(ukb, V1, V4)
sincr <- hrs_afr[ukb, nomatch = 0]
assign(paste0("keep_chr_", chr), sincr)
cat(chr, "done\n")
}
KEEP <- rbind(keep_chr_1, keep_chr_2, keep_chr_3, keep_chr_4, 
keep_chr_5, keep_chr_6, keep_chr_7, keep_chr_8, keep_chr_9, 
keep_chr_10, keep_chr_11, keep_chr_12, keep_chr_13, keep_chr_14, 
keep_chr_15, keep_chr_16, keep_chr_17, keep_chr_18, keep_chr_19, 
keep_chr_20, keep_chr_21, keep_chr_22)
fwrite(as.data.frame(KEEP[, i.V2]), file = "~/height_prediction/ldpred/output/keep_hrs_afr_RS.txt", quote = F, col.names = F)
remove(list = ls())

for (chr in 22:1) {
ukb <- fread(paste0("/project/mathilab/data/UKB/imputed/ukb_imp_chr", chr, "_eur.bim"))
whi <- fread("output/WHI.ldpred.bim")[V1 == chr]
whi$V1 <- as.integer(whi$V1)
setkey(ukb, V1, V4)
setkey(whi, V1, V4)
sincr <- whi[ukb, nomatch = 0]
assign(paste0("keep_chr_", chr), sincr)
cat(chr, "done\n")
}
KEEP <- rbind(keep_chr_1, keep_chr_2, keep_chr_3, keep_chr_4, 
keep_chr_5, keep_chr_6, keep_chr_7, keep_chr_8, keep_chr_9, 
keep_chr_10, keep_chr_11, keep_chr_12, keep_chr_13, keep_chr_14, 
keep_chr_15, keep_chr_16, keep_chr_17, keep_chr_18, keep_chr_19, 
keep_chr_20, keep_chr_21, keep_chr_22)

fwrite(as.data.frame(KEEP[, i.V2]), file = "~/height_prediction/ldpred/output/keep_whi_RS.txt", quote = F, col.names = F)
remove(list = ls())
for (chr in 22:1) {
ukb <- fread(paste0("/project/mathilab/data/UKB/imputed/ukb_imp_chr", chr, "_eur.bim"))
jhs <- fread("output/JHS.ldpred.bim")[V1 == chr]
jhs$V1 <- as.integer(jhs$V1)
setkey(jhs, V1, V4)
setkey(ukb, V1, V4)
sincr <- jhs[ukb, nomatch = 0]
assign(paste0("keep_chr_", chr), sincr)
cat(chr, "done\n")
}

KEEP <- rbind(keep_chr_1, keep_chr_2, keep_chr_3, keep_chr_4,
keep_chr_5, keep_chr_6, keep_chr_7, keep_chr_8, keep_chr_9,
keep_chr_10, keep_chr_11, keep_chr_12, keep_chr_13, keep_chr_14,
keep_chr_15, keep_chr_16, keep_chr_17, keep_chr_18, keep_chr_19,
keep_chr_20, keep_chr_21, keep_chr_22)

fwrite(as.data.frame(KEEP[, i.V2]), file = "~/height_prediction/ldpred/output/keep_jhs_RS.txt", quote = F, col.names = F)
