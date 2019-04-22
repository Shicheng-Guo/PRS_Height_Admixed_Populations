*Select intersection of UKBiobank and WHI SNPs*
```
cd scripts/
Rscript --vanilla make_vcf.R temp
cd ..
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
Rscript --vanilla Plots_JHS.R
cd ..
```

*Run ancestry specific beta estimation[ongoing]*
```
for chr in {1..22};
do
#bsub -M 20000 Rscript --vanilla /project/mathilab/bbita/gwas_admix/new_height/WHI/run_AS_GWAS.R ${chr}
Rscript --vanilla run_AS_GWAS.R ${chr}
echo ${chr}
done
```
*Plot results for ancestry specific beta estimation[ongoing]*
```
cd scripts/
for i in 100 200 1000 2000 10000 20000 100000 200000;
do
#bsub -M 40000 -o logLocAnc${i} -e logLocAnc${i} Rscript --vanilla /project/mathilab/bbita/gwas_admix/new_height/WHI/plot_diff_map.R phys_100000_0.0005 ${i}
Rscript --vanilla plot_diff_map.R phys_100000_0.0005 ${i}
echo ${i}
done
cd ..
```
