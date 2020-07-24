## Getting HRS_eur data ready for all downstream analyses

Goal: go from plink to vcf format. 

*copy files to my local directory*
```
PATH_to_data=/project/mathilab/data/HRS/ #change this accordingly.

cp ${PATH_to_data}/data/HRS_EUR_b37_strand_include* . #plink files
cp ${PATH_to_data}/data/HRS_EUR_phenotypes.txt .
#cp ${PATH_to_data}/admixture/HRS_AFR_b37_strand_prune_include.2.Q . ##order comes from file below
cp ${PATH_to_data}/data/HRS_EUR_b37_strand_include.fam .  #the order of samples
```
*convert to vcf format*

```
plink --allow-extra-chr --allow-no-sex --autosome \
--bfile HRS_EUR_b37_strand_include \
--recode vcf \
--keep-allele-order \
--out HRS_EUR_b37_strand_include.bas
```

```
awk '/^ *#/ { print; }' HRS_EUR_b37_strand_include.bas.vcf > header_HRS_eur.txt

for chr in {1..22};
do
touch chr${chr}_bas.vcf
cat header_HRS_eur.txt >  chr${chr}_bas.vcf
awk -v i=$chr '$1==i' HRS_EUR_b37_strand_include.bas.vcf >> chr${chr}_bas.vcf
echo 'chr'
echo $chr
echo 'done'
done
```

*compress/index*
```
for chr in {1..22};
do
bgzip -c chr${chr}_bas.vcf > chr${chr}_bas.vcf.gz &&
tabix -p vcf chr${chr}_bas.vcf.gz
done

bgzip -c HRS_EUR_b37_strand_include.bas.vcf > HRS_EUR_b37_strand_include.bas.vcf.gz
