## Getting HRS_afr data ready for all downstream analyses

Goal: go from plink to vcf format.

*copy files to my local directory*
```
PATH_to_data=/project/mathilab/data/UKB #change this accordingly.

*copy files to my local directory*
```
cp ${PATH_to_data}/UKB_EUR* .
rm UKB_EUR.afreq
```
*convert to vcf format*


```
plink --allow-extra-chr --allow-no-sex --autosome \
--bfile UKB_EUR \
--recode vcf \
--keep-allele-order \
--out UKB_EUR.bas


head -n 30 UKB_EUR.bas.vcf | awk '/^ *#/ { print; }'  > header_ukb_eur.txt

for chr in {1..22};
do
touch chr${chr}_bas_eur.vcf
cat header_ukb_eur.txt >  chr${chr}_bas_eur.vcf
awk -v i=$chr '$1==i' UKB_EUR.bas.vcf >> chr${chr}_bas_eur.vcf
echo 'chr'
echo $chr
echo 'done'
done

for chr in {1..22};
do
bgzip chr${chr}_bas_eur.vcf  &&
tabix -p vcf chr${chr}_bas_eur.vcf.gz
done
```


