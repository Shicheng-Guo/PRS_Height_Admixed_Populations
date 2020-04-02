*Prepare input files for UKB_AFR PCA*

Note: we use the non-imputed dataset to get the PCAs but we use the imputed data for the GWAS (see PCA_and_GWAS/UKB_AFR_imputed)
```
awk '{b=$1"_"$1;print b,"UKA"}' /project/mathilab/data/UKB/UKB_AFR.fam|grep -v '6007195'  > panel_file2.txt
touch "pop2.txt"
echo "UKA" > pop2.txt
#edit arguments file
#remove individual '6007195':
awk 'NR==28{print $0}' ~/height_prediction/input/ukb_afr/UKB_AFR.bas.vcf|tr "\t" "\n"|grep -n '6007195' #line 8801
awk '{$8801="";print $0; OFS="\t"}' ~/height_prediction/input/ukb_afr/UKB_AFR.bas.vcf > test.vcf
bgzip test.vcf
mv test.vcf.gz UKB_AFR.bas.vcf.gz

bsub -n 10 -M 40000 -R "span [hosts=1] rusage [mem=40480]" -o logPCA -e logPCA bash ~/height_prediction/PCA_and_GWAS/generatePCA.sh ~/height_prediction/PCA_and_GWAS/UKB_AFR/generatePCA_arg.txt
```

