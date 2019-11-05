#Perform a local ancestry weughted PRS.

```
Rscript --vanilla combine_plink_files.R
```

```
for chr in {1..22};
do
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.05
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.1
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.15
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.2
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.3
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.4
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.5
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.6
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.7
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.8
bsub -M 30000 -o ~/height_prediction/loc_anc_analysis/split_chr${chr} -e ~/height_prediction/loc_anc_analysis/split_chr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_PRS_v2.R phys_100000_0.0005 ${chr} 0.9
done
```

```
run_partial_r2.r
```
