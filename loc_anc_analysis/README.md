#Perform a local ancestry specific PRS

```
##Format plink output for UKB_afr (imputed) GWAS
Rscript --vanilla combine_plink_files.R
```

```
#Run local ancestry specific PRS for each chromosome and using different alpha values: Rscript --vanilla run_PRS_v2.R <pruned_set> <chr> <alpha>
for chr in {1..22};
do
for alpha in $(seq 0 0.01 1);
do
bsub -M 60000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} ${alpha}
done
done
```

```
#Evaluate prediction power
run_partial_r2.r
```


```
#Run local ancestry specific PRS for each chromosome and using different alpha values: Rscript --vanilla run_PRS_v3.R <pruned_set> <chr> <alpha>
#This version uses only SNPs in EUR ancestry blocks.
for chr in {1..22};
do
bsub -M 60000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v3.R phys_100000_0.0005 ${chr} 0
done
```


```
#JHS

for chr in {1..22};
do
for alpha in $(seq 0 0.01 1);
do
bsub -M 60000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2_JHS.R phys_100000_0.0005 ${chr} ${alpha}
done
done

for chr in {1..22};
do
bsub -M 60000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v3_JHS.R phys_100000_0.0005 ${chr} 0
done
```


```
#HRS
for chr in {1..22};
do
for alpha in $(seq 0 0.01 1);
do
bsub -M 50000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2_HRS_afr.R phys_100000_0.0005 ${chr} ${alpha}
done
done

for chr in {1..22};
do
bsub -M 60000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v3_HRS_afr.R phys_100000_0.0005 ${chr} 0
done
```


```
#UKB_afr

for chr in {1..22};
do
bsub -M 60000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v3_UKB_afr.R phys_100000_0.0005 ${chr} 0
done



```
#Plots:
#Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_partial_r2_v3.r #makes partial-r2vsEUR_ANC plots using LA PRS using only EUR segments
#Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_partial_r2_v2.R

```
