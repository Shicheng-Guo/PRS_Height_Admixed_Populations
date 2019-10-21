#Perform a local ancestry weughted PRS.

```
Rscript --vanilla combine_plink_files.R
```

```
for chr in {1..22};
bsub -M 20000 Rscript --vanilla run_LA_PRS.r phys_100000_0.0005 ${chr}
done
```

```
for chr in {1..22};
do
bsub -M 80000 -o ~/height_prediction/loc_anc_analysis/logchr${chr} -e ~/height_prediction/loc_anc_analysis/logchr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_LA_PRS.r phys_100000_0.0005 ${chr} 0.10
bsub -M 80000 -o ~/height_prediction/loc_anc_analysis/logchr${chr} -e ~/height_prediction/loc_anc_analysis/logchr${chr} Rscript --vanilla ~/height_prediction/loc_anc_analysis/run_LA_PRS.r phys_100000_0.0005 ${chr} 0.20
done
```
