#Perform a local ancestry specific PRS

```
##Format plink output for UKB_afr (imputed) GWAS
Rscript --vanilla combine_plink_files.R
```

```
#Run local ancestry specific PRS for each chromosome and using different alpha values: Rscript --vanilla run_PRS_v2.R <pruned_set> <chr> <alpha>
for chr in {1..22};
do
for alpha in 0 0.05 0.1 0.15 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1;
do
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} alpha
done
done
```

```
#Evaluate prediction power
run_partial_r2.r
```
