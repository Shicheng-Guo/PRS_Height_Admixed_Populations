#local ancestry

for i in {1..22};
do
bsub -M 20000 Rscript --vanilla ~/height_prediction/gwas/ukb_afr/scripts/run_AS_GWAS.R ${i}
done

