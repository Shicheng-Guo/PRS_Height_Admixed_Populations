####
*Select SNPs and prepare the data*
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp sib_betas $D
done
``
*Prune using different methods* #[ongoing]
```
*Prune

```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
~/height_prediction/scripts/LD_prun.bash sib_betas $D
done
```

*Combine
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
~/height_prediction/scripts/combine_Rds_v2.sh sib_betas $D
done
```

*Run polygenic scores*

```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
~/height_prediction/scripts/calc_PGS.sh sib_betas $D
done
```

*Combine PGS results

for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
~/height_prediction/scripts/combine_Rds_PGS.sh sib_betas $D
done
```
