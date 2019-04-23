Here we evaluate the effect of recombination rates on prediction power of SNPs used in the PRS.

*Stratify genome into 4 quantiles of rec rate (for AA and CEU maps, separately).

*Using one chosen pruning strategy, take only pruned SNPs within a given quantile and calculate PRS

```
for i in strat_prs_ukb_afr.R strat_prs_ukb_eur.R strat_prs_penn_afr.R strat_prs_penn_eur.R strat_prs_JHS.R strat_prs_WHI.R;
do
for k in 5000 10000 20000;
do
bsub -M 60000 -o /project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/logs/ Rscript --vanilla /project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/scripts/${i} ${k} AA
bsub -M 60000 -o /project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/logs/ Rscript --vanilla /project/mathilab/bbita/gwas_admix/height_prediction/strat_prs/scripts/${i} ${k} CEU
done
done
```

See if there are differences.

```
cd scripts
Rscript --vanilla Plot.R'
```



