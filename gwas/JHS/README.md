for i in {1..22}; do bsub -M 80000 Rscript --vanilla ~/height_prediction/gwas/JHS/scripts/make_anc.R ${i}; done
