## Preparing genotype data

It is assumed and required that you have access to the genotype data used in this work. We are not allowed to redistribute it, but they are all available through dbgap.


## Health and Retirement Study (HRS)

Follow instructions for [HRS_eur](HRS_eur/README.md)  and [HRS_afr](HRS_afr/README.md). 
## UK Biobank (UKB)

Follow instructions for [UKB_eur](ukb_eur/README.md)  and [UKB_afr](ukb_afr/README.md).
## Women's Health Initivave (WHI)

Follow instructions for [WHI_afr](WHI/README.md).

## Jackson Heart Study (JHS)
Follow instructions for [JHS_afr](JHS/README.md).


## Summary statistics

Download
```
wget https://www.dropbox.com/s/od6dr8kdrrornuz/50_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 #  50_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0
```

Filter out low confidence variants

```
zcat 50_raw.gwas.imputed_v3.both_sexes.tsv.bgz\?dl\=0 |head -1 > header.txt
touch 50_raw_filtered.txt
cat header.txt > 50_raw_filtered.txt
zgrep FALSE 50_raw.gwas.imputed_v3.both_sexes.tsv.bgz\?dl\=0 >> 50_raw_filtered.gz
Rscript --vanilla ../scripts/format_ukbb_height.R #generates input files for plink.

```
