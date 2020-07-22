##Try to replicate R2xEUR_ANC plots using LDpred.
##FOr now, using this tutorial as reference: https://choishingwan.github.io/PRS-Tutorial/ldpred/
##And also the ldpred tutorial here: 
source ~/ENV/bin/activate

1. Preprocessing the base data file and coordinating the data
The first step is a data synchronization step, where two or three data sets, genotypes and summary statistics are synchronized. This generates a HDF5 file which contains the synchronized genotypes. This step can be done by running

```
ldpred coord
```

use --help for detailed options. This step requires at least one genotype file (the LD reference genotypes), where we recommend at least 1000 unrelated individuals with the same ancestry make-up as the individuals for which summary statistics datasets are obtained from. Another genotype file can also be given if the user intends to validate the predictions using a separate set of genotypes.

```
1000G EUR (except FIN) as LD panel

MY_PATH=$(pwd)

grep EUR /project/mathilab/data/1kg/20130502_phase3_final/integrated_call_samples_v3.20130502.ALL.panel |grep -v FIN|awk 'OFS="\t"{print $1}' > EUR_samples.txt
cat EUR_samples.txt | tr ” ” “\n” > EUR_sample_byline.txt

```

```
PATH_TO=/project/mathilab/data/1kg/20130502_phase3_final
for chr in {22..1};
do
bcftools view -Oz -S EUR_sample_byline.txt -m2 -M2 -v snps ${PATH_TO}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz |bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz > output/1000g_chr${chr}.vcf.gz
done

#combine all chromosomes
bcftools concat output/1000g_chr1.vcf.gz output/1000g_chr2.vcf.gz output/1000g_chr3.vcf.gz output/1000g_chr4.vcf.gz output/1000g_chr5.vcf.gz output/1000g_chr6.vcf.gz output/1000g_chr7.vcf.gz output/1000g_chr8.vcf.gz output/1000g_chr9.vcf.gz output/1000g_chr10.vcf.gz output/1000g_chr11.vcf.gz output/1000g_chr12.vcf.gz output/1000g_chr13.vcf.gz output/1000g_chr14.vcf.gz output/1000g_chr15.vcf.gz output/1000g_chr16.vcf.gz output/1000g_chr17.vcf.gz output/1000g_chr18.vcf.gz output/1000g_chr19.vcf.gz output/1000g_chr20.vcf.gz output/1000g_chr21.vcf.gz output/1000g_chr22.vcf.gz -Oz > output/1000g_all.vcf.gz

#convert to plink
plink --vcf output/1000g_all.vcf.gz --keep-allele-order --make-bed --out output/1000g_all
#change ID column to be compatible with sumstats
awk '{$2=$1"_"$4;print $0}' output/1000g_all.bim > output/tmp && mv output/tmp output/1000g_all.bim


#UKB_EUR imputed (1000g does not have enough people for LD panel)

Rscript --vanilla ../scripts/parse_ukb_imputed.r #generates files with positions to keep

#print only IDs that are not duplicated
sort output/keep_hrs_eur_RS.txt |uniq -u > output/tmp && mv output/tmp output/keep_hrs_eur_RS.txt
sort output/keep_hrs_afr_RS.txt |uniq -u > output/tmp && mv output/tmp output/keep_hrs_afr_RS.txt
sort output/keep_whi_RS.txt |uniq -u > output/tmp && mv output/tmp output/keep_whi_RS.txt
sort output/keep_jhs_RS.txt |uniq -u > output/tmp && mv output/tmp output/keep_jhs_RS.txt

PATH2=project/mathilab/data/UKB/imputed
for chr in {22..1}; 
do
plink --bfile ${PATH2}/ukb_imp_chr${chr}_eur --extract output/keep_hrs_eur_RS.txt --make-bed --out output/UKB_EUR_imp_chr${chr}
done

for chr in {22..1};
do
plink --bfile ${PATH2}/ukb_imp_chr${chr}_eur --extract output/keep_whi_RS.txt --make-bed --out output/UKB_EUR_whi_imp_chr${chr}
done


for chr in {22..1};
do
plink --bfile ${PATH2}/imputed/ukb_imp_chr${chr}_eur --extract output/keep_jhs_RS.txt --make-bed --out output/UKB_EUR_jhs_imp_chr${chr}
done

rm output/mergelist.txt
for chr in {1..22}
do
echo ${MY_PATH}/output/UKB_EUR_jhs_imp_chr${chr} >> output/mergelist.txt
done

bsub -M 95000 -o ${MY_PATH}/logs/logplink -e ${MY_PATH}/logs/logplink "plink --merge-list ${MY_PATH}/output/mergelist.txt --make-bed --out ${MY_PATH}/output/UKB_EUR_jhs_imp_all"

rm output/mergelist.txt
for chr in {1..22}
do
echo ${MY_PATH}/output/UKB_EUR_whi_imp_chr${chr} >> ${MY_PATH}/output/mergelist.txt
done

bsub -M 95000 -o ${MY_PATH}/logs/logplink -e ${MY_PATH}/logs/logplink "plink --merge-list ${MY_PATH}/output/mergelist.txt --make-bed --out ${MY_PATH}/ldpred/output/UKB_EUR_whi_imp_all"


awk 'OFS="\t"{$2=$1"_"$4;print $0}' output/UKB_EUR_jhs_imp_all.bim > output/tmp && mv output/tmp output/UKB_EUR_jhs_imp_all.bim
sed -i 's/\s/\t/g'  output/UKB_EUR_jhs_imp_all.bim

awk 'OFS="\t"{$2=$1"_"$4;print $0}' output/UKB_EUR_whi_imp_all.bim > output/tmp && mv output/tmp output/UKB_EUR_whi_imp_all.bim
sed -i 's/\s/\t/g'  output/UKB_EUR_whi_imp_all.bim

```
```
#Summary statistics file
Rscript --vanilla ../scripts/format_sumstat.R #format summary statistics file

#LDpred does not support filtering of samples and SNPs, so therefore we must generate a new QCed genotype file using plink:

plink2 \
    --bfile /project/mathilab/data/UKB/UKB_EUR \ 
    --keep /project/mathilab/data/UKB/UKB_EUR_IDS \
    --make-bed \
    --out  ${MY_PATH}/output/UKB_EUR.ldpred 
#
#UKB_AFR

plink2 \
    --bfile /project/mathilab/data/UKB/UKB_AFR \
    --keep /project/mathilab/data/UKB/UKB_AFR_IDS \
    --make-bed \
    --out ${MY_PATH}/output/UKB_AFR.ldpred
#
#HRS_AFR

plink2 \
    --bfile /project/mathilab/data/HRS/data/HRS_AFR_b37_strand_include \
    --keep /project/mathilab/data/HRS/data/HRS_AFR_IDS.fam \
    --make-bed \
    --out ${MY_PATH}/output/HRS_AFR.ldpred
#
#HRS_EUR

plink2 \
    --bfile /project/mathilab/data/HRS/data/HRS_EUR_b37_strand_include \
    --keep /project/mathilab/data/HRS/data/HRS_EUR_IDS.fam \
    --make-bed \
    --out ${MY_PATH}/output/HRS_EUR.ldpred
#WHI

plink2 \
    --bfile /project/mathilab/data/WHI/data/WHI_b37_strand_include \
    --make-bed \
     --autosome \
    --out ${MY_PATH}/output/WHI.ldpred

#JHS

plink2 \
    --bfile JHS_b37_strand \
    --make-bed \
    --autosome \
    --out ${MY_PATH}/output/JHS.ldpred
```

```
**COORD**
#Preprocessing the base data file:

# There are 360,388 samples in the Height GWAS
Rscript --vanilla merge_hrs_ukb_v2.R #change SNP ID column
#UKB_EUR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 360388  \
    --vgf ${MY_PATH}/UKB_EUR.ldpred \
    --ssf ${MY_PATH}/output/Height.QC.gz \ 
    --out ${MY_PATH}/output/UKB_EUR.coord \
    --gf ${MY_PATH}/output/UKB_EUR.ldpred > ${MY_PATH}/logs/log_coord_ukb_eur.log
#HRS_AFR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 360388 \
    --vgf ${MY_PATH}/output/HRS_AFR.ldpred \
    --ssf ${MY_PATH}/output/Height.QC.gz \
    --out ${MY_PATH}/output/HRS_AFR.coord \
    --gf ${MY_PATH}/output/UKB_EUR_imp_all > ${MY_PATH}/logs/log_coord_hrs_afr.log
#HRS_EUR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 360388 \
    --vgf ${MY_PATH}/output/HRS_EUR.ldpred \
    --ssf ${MY_PATH}/output/Height.QC.gz \
    --out ${MY_PATH}/output/HRS_EUR.coord \
    --gf  ${MY_PATH}/output/UKB_EUR_imp_all  > ~/height_prediction/ldpred/logs/log_coord_hrs_eur.log
#UKB_AFR
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 360388 \
    --vgf ${MY_PATH}/output/UKB_AFR.ldpred \
    --ssf ${MY_PATH}/output/Height.QC.gz \
    --out ${MY_PATH}/output/UKB_AFR.coord \
    --gf ${MY_PATH}/output/UKB_EUR.ldpred > ${MY_PATH}/logs/log_coord_ukb_afr.log
#WHI
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 360388 \
    --vgf ${MY_PATH}/output/WHI.ldpred \
    --ssf ${MY_PATH}/output/Height.QC.gz \
    --out ${MY_PATH}/output/WHI.coord \
    --gf  ${MY_PATH}/output/UKB_EUR_whi_imp_all > ${MY_PATH}/logs/log_coord_whi.log
#JHS
ldpred coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos POS \
    --chr CHR \
    --pval PVAL \
    --eff BETA \
    --ssf-format CUSTOM \
    --eff_type LINREG\
    --N 360388 \
    --vgf ${MY_PATH}/output/JHS.ldpred \
    --ssf ${MY_PATH}/output/Height.QC.gz \
    --out ${MY_PATH}/output/JHS.coord \
    --gf ${MY_PATH}/output/UKB_EUR_jhs_imp_all > ${MY_PATH}/logs/log_coord_jhs.log
```

2. Adjust the effect size estimates

```
# LDpred recommend radius to be Total number of SNPs in target / 3000. E.g. 784256/3000=261
#Regarding choice of the LD panel, its LD structure should ideally be similar to the training data for which the summary statistics are calculated.(UKB_EUR)

#UKB_EUR 506723/3000~169
ldpred gibbs  --cf ${MY_PATH}/output/UKB_EUR.coord  --ldr 169 --ldf ${MY_PATH}/output/UKB_EUR.ldpred --out ${MY_PATH}/output/UKB_EUR.weight
#UKB_AFR 479048/3000~159
ldpred gibbs  --cf ${MY_PATH}/output/UKB_AFR.coord  --ldr 159 --ldf ${MY_PATH}/output/UKB_AFR.ldpred --out ${MY_PATH}/output/UKB_AFR.weight
#HRS_eur 1058065/3000~357
ldpred gibbs  --cf ${MY_PATH}/output/HRS_EUR.coord  --ldr 357 --ldf ${MY_PATH}/output/HRS_EUR.ldpred --out ${MY_PATH}/output/HRS_EUR.weight
#HRS_afr 960309/3000~320
ldpred gibbs  --cf ${MY_PATH}/output/HRS_AFR.coord  --ldr 320 --ldf ${MY_PATH}/output/HRS_AFR.ldpred --out ${MY_PATH}/output/HRS_AFR.weight
#WHI 735278/3000 ~ 246
ldpred gibbs  --cf ${MY_PATH}/output/WHI.coord  --ldr 245 --ldf ${MY_PATH}/output/WHI.ldpred --out ${MY_PATH}/output/WHI.weight
#JHS 699041/3000 233
ldpred gibbs  --cf ${MY_PATH}/output/JHS.coord --ldr 233 --ldf ${MY_PATH}/output/JHS.ld --out ${MY_PATH}/output/JHS.weight

#P+T
ldpred p+t --cf ${MY_PATH}/output/UKB_EUR.coord  --ldr 169 --out ${MY_PATH}/output/UKB_EUR_pt.weight > ${MY_PATH}/logs/log_pt_ukb_eur.log
ldpred p+t --cf ${MY_PATH}/output/HRS_EUR.coord  --ldr 357 --out ${MY_PATH}/output/HRS_EUR_pt.weight  > ${MY_PATH}/logs/log_pt_hrs_eur.log
ldpred p+t --cf ${MY_PATH}/output/HRS_AFR.coord  --ldr 320 --out ${MY_PATH}/output/HRS_AFR_pt.weight > ${MY_PATH}/logs/log_pt_hrs_afr.log
ldpred p+t --cf ${MY_PATH}/output/WHI.coord  --ldr 245 --out ${MY_PATH}/output/WHI_pt.weight > ${MY_PATH}/logs/log_pt_whi.log
ldpred p+t --cf ${MY_PATH}/output/JHS.coord --ldr 233 --out ${MY_PATH}/output/JHS_pt.weight > ${MY_PATH}/logs/log_pt_JHS.log
ldpred p+t --cf ${MY_PATH}/output/UKB_AFR.coord --ldr 159 q--out ${MY_PATH}/output/UKB_AFR_pt.weight  > ${MY_PATH}/logs/log_pt_ukb_afr.log
```

3. Calculate the PRS
```
#For just PRS (no validation), run below code (UKB_EUR)
ldpred score \
    --gf output/UKB_EUR.ldpred \
    --rf output/UKB_EUR.weight \
    --out output/UKB_EUR.score \
    --only-score 
ldpred score \
    --gf output/UKB_EUR.ldpred \
    --rf output/UKB_EUR_pt.weight \
    --out output/UKB_EUR.score \
    --only-score
#ukb afr
ldpred score \
    --gf output/UKB_AFR.ldpred \
    --rf output/UKB_AFR.weight \
    --out output/UKB_AFR.score \
    --only-score
ldpred score \
    --gf output/UKB_AFR.ldpred \
    --rf output/UKB_AFR_pt.weight \
    --out output/UKB_AFR.score \
    --only-score
#HRS_EUR
ldpred score \
    --gf output/HRS_EUR.ldpred \
    --rf output/HRS_EUR.weight \
    --out output/HRS_EUR.score \
    --only-score 
ldpred score \
    --gf output/HRS_EUR.ldpred \
    --rf output/HRS_EUR_pt.weight \
    --out output/HRS_EUR.score \
    --only-score
#HRS_afr
ldpred score \
    --gf output/HRS_AFR.ldpred \
    --rf output/HRS_AFR.weight \
    --out output/HRS_AFR.score \
    --only-score
ldpred score \
    --gf output/HRS_AFR.ldpred \
    --rf output/HRS_AFR_pt.weight \
    --out output/HRS_AFR.score \
    --only-score
#WHI
ldpred score \
    --gf output/WHI.ldpred \
    --rf output/WHI.weight \
    --out output/WHI.score \
    --only-score 
ldpred score \
    --gf output/WHI.ldpred \
    --rf output/WHI_pt.weight \
    --out output/WHI.score \
    --only-score
#JHS
ldpred score \
    --gf output/JHS.ldpred \
    --rf output/JHS.weight \
    --out output/JHS.score \
    --only-score 
ldpred score \
    --gf output/JHS.ldpred \
    --rf output/JHS_pt.weight \
    --out output/JHS.score \
    --only-score
```
