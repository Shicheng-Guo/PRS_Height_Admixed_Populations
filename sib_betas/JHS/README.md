*Select intersection of UKBiobank and WHI SNPs*
```
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp sib_betas JHS
~/height_prediction/scripts/combine_Rds_v2.sh sib_betas JHS
```
*Prune using different methods*
```
cd scripts/
bash LD_prun.bash
bash combine_Rds_v2.sh
bash run_this.sh
cd ..
```
*Run polygenic scores*
```
cd scripts/
bash calc_PGS.sh
bash combine_Rds_PGS.sh	
bash run_this_PGS.sh
Rscript --vanilla $cur/Plots_JHS.R
cd ..
```


