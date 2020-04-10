##

1. Take pruned vcf files for each dataset for some but not all clumping strategies.

2. Run PRS without beta weight (write a slightly modified script for this)

```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/unweighted_prs/calc_PGS.sh gwas $D
done
#combine
for D in JHS WHI ukb_afr ukb_eur HRS_afr HRS_eur; 
do 
rm ~/height_prediction/unweighted_prs/${D}/test2.txt; 
rm ~/height_prediction/unweighted_prs/${D}/run_this_PGS.sh;
~/height_prediction/unweighted_prs/combine_Rds_PGS.sh unweighted_prs $D
done
```

3. Run plots for each dataset

for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
bsub Rscript --vanilla /home/bbita/height_prediction/unweighted_prs/${D}/Plots_${D}.R
done
4. Combine all datasets and make a figure like Figure 1 (again, just modifiy previous script)

```
Rscript --vanilla combine_datasets.R
```

