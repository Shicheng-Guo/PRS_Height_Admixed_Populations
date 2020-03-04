Here we evaluate the effect of recombination rates on prediction power of SNPs used in the PRS.

*Stratify genome into 4 quantiles of rec rate (for AA and CEU maps, separately).

*Using one chosen pruning strategy, take only pruned SNPs within a given quantile and calculate PRS

```
for i in strat_prs_ukb_afr.R strat_prs_ukb_eur.R strat_prs_penn_afr.R strat_prs_penn_eur.R strat_prs_JHS.R strat_prs_WHI.R strat_prs_HRS_afr.R strat_prs_HRS_eur.R;
do
for k in 20000 40000 100000;
do
bsub -M 80000 Rscript --vanilla ~/height_prediction/strat_prs/scripts/${i} ${k} AA sib_betas
bsub -M 80000 Rscript --vanilla ~/height_prediction/strat_prs/scripts/${i} ${k} AA gwas
bsub -M 80000 Rscript --vanilla ~/height_prediction/strat_prs/scripts/${i} ${k} CEU sib_betas
bsub -M 80000 Rscript --vanilla ~/height_prediction/strat_prs/scripts/${i} ${k} CEU gwas
done
done
```

See if there are differences.

```
cd scripts
Rscript --vanilla Plot.R gwas
```

```
for i in strat_prs_ukb_afr_v2.R strat_prs_JHS_v2.R strat_prs_WHI_v2.R strat_prs_HRS_afr_v2.R strat_prs_HRS_eur_v2.R;
do
for k in 3000 6000 10000 20000 40000 100000;
do
bsub -M 90000 Rscript --vanilla ~/height_prediction/strat_prs/scripts/${i} ${k} AA gwas
bsub -M 90000 Rscript --vanilla ~/height_prediction/strat_prs/scripts/${i} ${k} CEU gwas
done
done
```
Plot:

```
cd scripts/
Rscript --vanilla Plots_v2.R
```

