#!/bin/bash


path="$1"
dtset="$2"
echo $path
#

for i in {1..22};
do
rm $path/output/hei_SNPs_chr${i}.vcf
touch $path/output/hei_SNPs_chr${i}.vcf
cat ~/height_prediction/input/header_${dtset}.txt |tail -n 1|sed 's/#//' >  $path/output/hei_SNPs_chr${i}.vcf
cat <(sort -n -k2 $path/temp/temp2_chr$i.txt) >> $path/output/hei_SNPs_chr${i}.vcf; #vcf files need to be sorted in order to be indexed.
echo 'sorting'
bgzip -c $path/output/hei_SNPs_chr${i}.vcf > $path/output/hei_SNPs_chr${i}.vcf.gz
echo 'compressing'
done

#
