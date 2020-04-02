for i in {1..22}; do bsub -M 80000 Rscript --vanilla ~/height_prediction/gwas/HRS_afr/scripts/make_anc.R ${i}; done
