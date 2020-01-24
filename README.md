Date created: April 18th 2019

In this directory you will find all scripts to recreate analyses described in our paper.

Intro

Height is a very polygenic trait and well-studied in humans. GWAS summary statistics for height used hundreds of thousands of individuals of European ancestry. It is unclear how well polygenic risk scores (PRS) predict height in non-Europeans in comparison to Europeans. Here we explore and quantify this, and then provide improvements.

*Subdirectories in this repo

WHI, JHS, etc: names of datasets

-gwas: analysis of predictive power of height PRS using UK Biobank GWAS summary statistics.

-sib_betas: analysis of predictive power of height PRS using effect sizes estimated from pairs of white British sibblings from the UK Biobank.

-strat_prs: analysis of predictive power of PRS as a function of recombination rates of SNPs.

*outfiles: where Rds files and such are stored, not to be pushed to repo.

*figs: where figures are stores, not to be pushed to repo.

*input files: where modified input files are stored, not to be pushed to repo.

*scripts and READMEs: should be pushed to repo.

Note: This will be updated as needed.
####
*Select SNPs and prepare the data*
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp sib_betas $D
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp gwas $D
done
``
*Prune using different methods* 
```
*Prune

```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/LD_prun.bash sib_betas $D
~/height_prediction/scripts/LD_prun.bash gwas $D
done
```

*Combine
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/combine_Rds_v2.sh sib_betas $D
~/height_prediction/scripts/combine_Rds_v2.sh gwas $D
done
```

*Run polygenic scores*

```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/calc_PGS.sh sib_betas $D
~/height_prediction/scripts/calc_PGS.sh gwas $D
done
```

*Combine PGS results

for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/combine_Rds_PGS.sh sib_betas $D
~/height_prediction/scripts/combine_Rds_PGS.sh gwas $D
done
```


*Plots

bsub -M 100000 Rscript --vanilla ~/height_prediction/sib_betas/WHI/scripts/Plots_WHI.R
bsub -M 100000 -e ~/height_prediction/gwas/WHI/logs/logplot -o ~/height_prediction/gwas/WHI/logs/logplot Rscript --vanilla ~/height_prediction/gwas/WHI/scripts/Plots_WHI.R
bsub -M 100000 Rscript --vanilla ~/height_prediction/sib_betas/JHS/scripts/Plots_JHS.R
