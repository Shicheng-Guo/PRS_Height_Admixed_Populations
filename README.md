**In this directory you will find all scripts to recreate analyses described in our paper**

*Intro*

Height is a very polygenic trait and well-studied in humans. GWAS summary statistics for height used hundreds of thousands of individuals of European ancestry. It is unclear how well polygenic risk scores (PRS) predict height in non-Europeans in comparison to Europeans. Here we explore and quantify this, and then provide improvements.

*Subdirectories in this repo*

-WHI, JHS, etc: names of datasets

-gwas: analysis of predictive power of height PRS using UK Biobank GWAS summary statistics.

-sib_betas: analysis of predictive power of height PRS using effect sizes estimated from pairs of white British sibblings from the UK Biobank.

-strat_prs: analysis of predictive power of PRS as a function of recombination rates of SNPs.

-unweighted_prs: analysis using an unweighted version of the PRS calculation, where effect sizes are ignored and replaced by +1 or -1 (for positive and negative effects, respectively)

-PCA_and_GWAS:PCA analysis of UKB_afr and GWAS for height in UKB_afr_imputed (subdirectory)

*output: where Rds files and such are stored, not to be pushed to repo.

*figs: where figures are stored, not to be pushed to repo.

*input files: where modified input files are stored, not to be pushed to repo.

*scripts and READMEs: should be pushed to repo.

Note: This will be updated as needed.
####
**RECREATING ANALYSES IN THE PAPER:**

*Prepare input data*

```
cd input/
```

Within each dataset's directory you will find a README.md with instructions on how to prepare input data. 

**Note: there is an assumption that you have access to the datasets used in this paper. We are not allowed to share the raw data.

*Get data ready for clumping/pruning*

Once the input data is formatted, we can do some pruning/clumping using both the GWAS effect sizes ('gwas') and the sibling-estimated effect sizes ('sib_betas'). In both cases, p-values used for clumping come from thefull UKB GWAS.

```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;  #for each dataset
do
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp sib_betas $D
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp gwas $D
done
```

*Prune/clump using different methods* 
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/LD_prun.bash sib_betas $D
~/height_prediction/scripts/LD_prun.bash gwas $D
done
```

*Combine these 80 sets*

```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/combine_Rds_v2.sh sib_betas $D
~/height_prediction/scripts/combine_Rds_v2.sh gwas $D
done
```

*Run polygenic scores*
Now we are ready to calculate polygenic risk scores for each set of SNPs (gwas, sib_betas and unweighted_prs which, as the name suggest, is the unweighted version of the PRS):

```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/calc_PGS.sh sib_betas $D
~/height_prediction/scripts/calc_PGS.sh gwas $D
~/height_prediction/unweighted_prs/calc_PGS.sh gwas $D
done
```

*Combine PRS results*

Combine all PRS results per dataset:
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur HRS_eur HRS_afr;
do
~/height_prediction/scripts/combine_Rds_PGS.sh sib_betas $D
~/height_prediction/scripts/combine_Rds_PGS.sh gwas $D
~/height_prediction/unweighted_prs/combine_Rds_PGS.sh unweighted_prs $D
done
```


*Plots*

```
Rscript --vanilla ~/height_prediction/sib_betas/WHI/scripts/Plots_WHI.R
Rscript --vanilla ~/height_prediction/sib_betas/JHS/scripts/Plots_JHS.R
Rscript --vanilla ~/height_prediction/sib_betas/ukb_afr/scripts/Plots_ukb_afr.R
Rscript --vanilla ~/height_prediction/sib_betas/ukb_eur/scripts/Plots_ukb_eur.R
Rscript --vanilla ~/height_prediction/sib_betas/HRS_afr/scripts/Plots_HRS_afr.R
Rscript --vanilla ~/height_prediction/sib_betas/HRS_eur/scripts/Plots_HRS_eur.R
Rscript --vanilla ~/height_prediction/gwas/WHI/scripts/Plots_WHI.R
Rscript --vanilla ~/height_prediction/gwas/JHS/scripts/Plots_JHS.R
Rscript --vanilla ~/height_prediction/gwas/ukb_afr/scripts/Plots_ukb_afr.R
Rscript --vanilla ~/height_prediction/gwas/ukb_eur/scripts/Plots_ukb_eur.R
Rscript --vanilla ~/height_prediction/gwas/HRS_afr/scripts/Plots_HRS_afr.R
Rscript --vanilla ~/height_prediction/gwas/HRS_eur/scripts/Plots_HRS_eur.R
Rscript --vanilla ~/height_prediction/unweighted_prs/WHI/scripts/Plots_WHI.R
Rscript --vanilla ~/height_prediction/unweighted_prs/JHS/scripts/Plots_JHS.R
Rscript --vanilla ~/height_prediction/unweighted_prs/ukb_afr/scripts/Plots_ukb_afr.R
Rscript --vanilla ~/height_prediction/unweighted_prs/ukb_eur/scripts/Plots_ukb_eur.R
Rscript --vanilla ~/height_prediction/unweighted_prs/HRS_afr/scripts/Plots_HRS_afr.R
Rscript --vanilla ~/height_prediction/unweighted_prs/HRS_eur/scripts/Plots_HRS_eur.R
```


*Combine datasets into Fig 1, Fig S6, Fig S7*

These scripts will produce plots for each pruning/clumping strategy. Throughout the paper we show the one called "phys_100000_0.0005":

```
Rscript --vanilla ~/height_prediction/scripts/combine_datasets.R gwas
Rscript --vanilla ~/height_prediction/combine_datasets.R  sib_betas
Rscript --vanilla~/height_prediction/unweighted_prs/combine_datasets.R
```
