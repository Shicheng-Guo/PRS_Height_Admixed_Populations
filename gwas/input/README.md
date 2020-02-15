wget https://www.dropbox.com/s/od6dr8kdrrornuz/50_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 #  50_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0

#filter out low confidence variants
zcat 50_raw.gwas.imputed_v3.both_sexes.tsv.bgz\?dl\=0 |head -1 > header.txt
touch 50_raw_filtered.txt
cat header.txt > 50_raw_filtered.txt
zgrep FALSE 50_raw.gwas.imputed_v3.both_sexes.tsv.bgz\?dl\=0 >> 50_raw_filtered.gz
