##Getting WHI_afr data ready for all downstream analyses

Goal: go from plink to vcf format.

*copy files to my local directory*
```
PATH_to_data=/project/mathilab/data/WHI #change this accordingly.

*copy files to my local directory*
```
cp ${PATH_to_data}/data/WHI_b37_strand_include* .
cp ${PATH_to_data}/data/WHI_phenotypes.txt .
cp ${PATH_to_data}/admxiture/WHI_b37_strand_prune_include.2.Q . ##order comes from file below
cp ${PATH_to_data}/data/WHI_b37_strand_prune_include.fam .  #the order of samples
```
*convert to vcf format*

```
plink --allow-extra-chr --allow-no-sex --autosome \
--bfile WHI_b37_strand_include \
--recode vcf \
--keep-allele-order \
--out WHI_b37_strand_include.bas
```

```
awk '/^ *#/ { print; }' WHI_b37_strand_include.bas.vcf > header.txt

for chr in {1..22};
do
touch chr${chr}_bas.vcf
cat header.txt >  chr${chr}_bas.vcf
awk -v i=$chr '$1==i' WHI_b37_strand_include.bas.vcf >> chr${chr}_bas.vcf
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

bgzip -c WHI_b37_strand_include.bas.vcf > WHI_b37_strand_include.bas.vcf.gz
```
