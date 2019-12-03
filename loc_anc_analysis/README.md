#Perform a local ancestry specific PRS

```
##Format plink output for UKB_afr (imputed) GWAS
Rscript --vanilla combine_plink_files.R
```

```
#Run local ancestry specific PRS for each chromosome and using different alpha values: Rscript --vanilla run_PRS_v2.R <pruned_set> <chr> <alpha>
for chr in {1..22};
do
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.05
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.1
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.15
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.2
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.3
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.4
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.5
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.6
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.7
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.8
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.9
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 1
done
```

```
#Evaluate prediction power
run_partial_r2.r
```
