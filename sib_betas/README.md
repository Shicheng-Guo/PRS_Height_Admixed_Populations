####
*Select SNPs and prepare the data*
```
for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
Rscript --vanilla ~/height_prediction/scripts/make_vcf.R temp sib_betas $D
done
``
*Prune using different methods* #[ongoing]
```
#combine stuff from above first

for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
~/height_prediction/scripts/LD_prun.bash sib_betas $D
done

for D in JHS WHI pennBB_afr pennBB_eur ukb_afr ukb_eur;
do
~/height_prediction/scripts/combine_Rds_v2.sh sib_betas $D
done
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
