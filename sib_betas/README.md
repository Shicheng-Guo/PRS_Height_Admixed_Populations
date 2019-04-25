####
*Select SNPs and prepare the data*
```
for D in WHI JHS pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp sib_betas $D
done
``
*Prune using different methods* #[ongoing]
```
#combine stuff from above first
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
Rscript --vanilla $cur/Plots_WHI.R
cd ..
```
