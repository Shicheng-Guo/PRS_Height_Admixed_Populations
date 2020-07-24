### Prune/clump, calculate PRS
*Select SNPs and prepare the data*
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
Rscript --vanilla ../scripts/make_vcf.R temp gwas $D
done
```

*Prune*

```
for D in JHS WHI ukb_afr ukb_eur  HRS_eur HRS_afr;;
do
bash ../scripts/LD_prun.bash gwas $D
done
```

*Combine*

```
for D in JHS WHI ukb_afr ukb_eur  HRS_eur HRS_afr;;
do
bash ../scripts/combine_Rds_v2.sh gwas $D
done
```


*Run polygenic scores*

```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
bash ../scripts/calc_PGS.sh gwas $D
done
```

*Combine PGS results*

```
for D in JHS WHI ukb_afr ukb_eur  HRS_eur HRS_afr;
do
bash ../scripts/combine_Rds_PGS.sh gwas $D
done
```
