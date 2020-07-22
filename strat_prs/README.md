Here we evaluate the effect of recombination rates on prediction power of SNPs used in the PRS.

*Stratify genome into 4 quantiles of rec rate (for AA and CEU maps, separately).

*Using one chosen pruning strategy, take only pruned SNPs within a given quantile and calculate PRS

```
MY_PATH=".."
for i in strat_prs_ukb_afr_v2.R strat_prs_JHS_v2.R strat_prs_WHI_v2.R strat_prs_HRS_afr_v2.R strat_prs_HRS_eur_v2.R; #for each script
do
for k in 3000 6000 10000 20000 40000 100000; #Window sizes for recombination rate
do
bsub -M 90000 Rscript --vanilla ${MY_PATH}/scripts/${i} ${k} AA gwas
bsub -M 90000 Rscript --vanilla ${MY_PATH}/scripts/${i} ${k} CEU gwas
done
done
```

Plot

```
cd scripts/
Rscript --vanilla ${MY_PATH}/scripts/Plot_strat_prs.R
```

