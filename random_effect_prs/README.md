##

1. Take pruned vcf files for each dataset for some but not all clumping strategies.

2. Run PRS without beta weight (write a slightly modified script for this)

```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/random_effect_prs/calc_PGS.sh gwas $D
done
#combine
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/random_effect_prs/combine_Rds_PGS.sh random_effect_prs $D
done
```
3. Combine all datasets and make a figure like Figure 1 (again, just modifiy previous script)
