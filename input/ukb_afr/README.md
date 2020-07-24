## Getting HRS_afr data ready for all downstream analyses

Goal: go from plink to vcf format.

*copy files to my local directory*
```
PATH_to_data=/project/mathilab/data/UKB #change this accordingly.


cp /project/mathilab/data/UKB/UKB_AFR* .
rm *prune*
#File exclude_biobank is list of samples that as of Nov 2018 should not be included anymore.

```
plink --allow-extra-chr --allow-no-sex --autosome \
--bfile UKB_AFR \
--recode vcf \
--keep-allele-order \
--out UKB_AFR.bas
```

head -n 30 UKB_AFR.bas.vcf | awk '/^ *#/ { print; }'  > header_ukb_afr.txt


for chr in {1..22};
do
touch chr${chr}_bas_afr.vcf
cat header_ukb_afr.txt >  chr${chr}_bas_afr.vcf
awk -v i=$chr '$1==i' UKB_AFR.bas.vcf >> chr${chr}_bas_afr.vcf
echo 'chr'
echo $chr
echo 'done'
done
```

```
for chr in {1..22};
do
bgzip chr${chr}_bas_afr.vcf  &&
tabix -p vcf chr${chr}_bas_afr.vcf.gz
done
```



