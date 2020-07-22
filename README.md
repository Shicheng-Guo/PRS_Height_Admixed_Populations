**In this directory you will find all scripts to recreate analyses described in our manuscript "Polygenic scores for height in admixed populations**

Link to manuscript: https://www.biorxiv.org/content/10.1101/2020.04.08.030361v2


*Intro*

Height is a very polygenic trait and well-studied in humans. GWAS summary statistics for height used hundreds of thousands of individuals of European ancestry. It is unclear how well polygenic risk scores (PRS) predict height in non-Europeans in comparison to Europeans. Here we explore and quantify this, and then provide improvements.


*Subdirectories in this repo*

-WHI, JHS, etc: names of datasets

-gwas: analysis of predictive power of height PRS using UK Biobank GWAS summary statistics.

-sib_betas: analysis of predictive power of height PRS using effect sizes estimated from pairs of white British sibblings from the UK Biobank.

-strat_prs: analysis of predictive power of PRS as a function of recombination rates of SNPs.

-unweighted_prs: analysis using an unweighted version of the PRS calculation, where effect sizes are ignored and replaced by +1 or -1 (for positive and negative effects, respectively)

-PCA_and_GWAS:PCA analysis of UKB_afr and GWAS for height in UKB_afr_imputed (subdirectory)

-ldpred: analyses with LDpred1 for comparison with C+T

-imputed: analyses with imputed genotype files when available (HRS, UKB)

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
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;  #for each dataset
do
Rscript --vanilla scripts/make_vcf.R temp sib_betas $D
Rscript --vanilla scripts/make_vcf.R temp gwas $D
done
```

*Prune/clump using different methods* 
```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
scripts/LD_prun.bash sib_betas $D
scripts/LD_prun.bash gwas $D
done
```

*Combine these 80 sets*

```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
scripts/combine_Rds_v2.sh sib_betas $D
scripts/combine_Rds_v2.sh gwas $D
done
```

*Run polygenic scores*
Now we are ready to calculate polygenic risk scores for each set of SNPs (gwas, sib_betas and unweighted_prs which, as the name suggest, is the unweighted version of the PRS):

```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
scripts/calc_PGS.sh sib_betas $D
scripts/calc_PGS.sh gwas $D;
unweighted_prs/calc_PGS.sh gwas $D
done
```
##

*LDpred analyses*

See [README.md in the ldpred directory](ldpred/README.md)

##

*Combine PRS results*

Combine all PRS results per dataset:
```
for D in JHS WHI ukb_afr ukb_eur HRS_eur HRS_afr;
do
rm gwas/${D}/scripts/test2.txt;
rm sib_betas/${D}/scripts/test2.txt;
rm gwas/${D}/scripts/run_this_PGS.sh;
rm sib_betas/${D}/scripts/run_this_PGS.sh;
rm unweighted_prs/${D}/run_this_PGS.sh;
rm unweighted_prs/${D}/test2.txt;
scripts/combine_Rds_PGS.sh sib_betas $D
scripts/combine_Rds_PGS.sh gwas $D
unweighted_prs/combine_Rds_PGS.sh unweighted_prs $D
done
```


*Plot & Combine datasets*

These scripts will produce plots for each pruning/clumping strategy. Throughout the paper we show the one called "phys_100000_0.0005":

```
for J in gwas sib_betas unweighted_prs;
do
Rscript --vanilla ${J}/WHI/scripts/Plots_WHI.R
Rscript --vanilla ${J}JHS/scripts/Plots_JHS.R
Rscript --vanilla ${J}ukb_afr/scripts/Plots_ukb_afr.R
Rscript --vanilla ${J}ukb_eur/scripts/Plots_ukb_eur.R
Rscript --vanilla ${J}HRS_afr/scripts/Plots_HRS_afr.R
Rscript --vanilla ${J}HRS_eur/scripts/Plots_HRS_eur.R
done
for J in gwas sib_betas;
do
Rscript --vanilla scripts/combine_datasets.R ${J
done
Rscript --vanilla unweighted_prs/combine_datasets.R
```


*Make all figures and supplementary tables in the paper*

```
Rscript --vanilla scripts/Fig1.R
Rscript --vanilla scripts/Fig2.R
Rscript --vanilla scripts/Fig3.R
Rscript --vanilla scripts/Fig4.R
Rscript --vanilla scripts/Fig5.R
Rscript --vanilla scripts/make_all_tables.R
```
