*select intersection of UKBiobank and WHI SNPs*
```
Rscript --vanilla scripts/make_vcf.R temp
```
*Prune using different methods*
```
bash scripts/LD_prun_sib.bash
bash scripts/combine_Rds_v2.sh
bash scripts/run_this.sh
```
*Run polygenic scores*
```
bash scripts/calc_PGS.sh
bash scripts/combine_Rds_PGS.sh	
bash scripts/run_this_PGS.sh
Rscript --vanilla sib_gwas/JHS/scripts/Plots_JHS.R
```


