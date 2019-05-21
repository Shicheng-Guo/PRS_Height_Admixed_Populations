for chr in {1..22};
do
cp ../../height/ukbiobank/ukb_data_run/chr${chr}_bas_afr.vcf.gz .
done

#File exclude_biobank is list of samples that as of Nov 2018 should not be included anymore.

cp /project/mathilab/bbita/gwas_admix/height/ukbiobank/ukb_data_run/*AFR* .
cp /project/mathilab/bbita/gwas_admix/height/ukbiobank/ukb_data_run/*afr* .


