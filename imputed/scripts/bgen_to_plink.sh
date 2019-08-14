for chr in {1..22};
do
plink2 --bgen /project/mathilab/data/UKB/imputed/ukb_imp_chr${chr}_eur.bgen --make-bed --out /project/mathilab/data/UKB/imputed/ukb_imp_chr${chr}_eur
plink2 --bgen /project/mathilab/data/UKB/imputed/ukb_imp_chr${chr}_afr.bgen --make-bed --out /project/mathilab/data/UKB/imputed/ukb_imp_chr${chr}_afr
done
