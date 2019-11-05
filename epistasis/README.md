#########################################################################################################################################################
## For this analysis, we want to test whether population-specific epistasis patterns drive the PRS predictive differences between eur and afr ancestries.
#########################################################################################################################################################


1. Take set of SNPs used in PRS
2. Define 2 or 20Kb window around each SNP
3. Compute causal allele frequency in each dataset
4. Compute the differences in these frequencies across datasets (particularly eur-afr pairs)
5. Check if (b_1-b_2)^2 is proportional to [sum((f_1-f_2)^2)]/N_snps


for i in WHI HRS_eur HRS_afr JHS ukb_afr ukb_eur;
do
Rscript --vanilla ~/height_prediction/epistasis/freq_window.R ${i}
done

