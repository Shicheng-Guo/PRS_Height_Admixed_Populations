*select intersection of UKBiobank and WHI SNPs*
```
Rscript --vanilla scripts/giant_sib.R temp
```
*Prune using different methods*

```
bash /project/mathilab/bbita/gwas_admix/height_prediction/sib_gwas/JHS/scripts/LD_prun_sib.bash
```


```
bash /project/mathilab/bbita/gwas_admix/height_prediction/sib_gwas/JHS/scripts/combine_Rds_v2.sh
bash /project/mathilab/bbita/gwas_admix/height_prediction/sib_gwas/JHS/scripts/run_this.sh
```
*Run polygenic scores*
```
bash /project/mathilab/bbita/gwas_admix/height_prediction/sib_gwas/JHS/scripts/calc_PGS.sh
bash /project/mathilab/bbita/gwas_admix/height_prediction/sib_gwas/JHS/scripts/combine_Rds_PGS.sh	
bash /project/mathilab/bbita/gwas_admix/height_prediction/sib_gwas/JHS/scripts/run_this_PGS.sh
Rscript --vanilla sib_gwas/JHS/scripts/Plots_JHS.R
```


