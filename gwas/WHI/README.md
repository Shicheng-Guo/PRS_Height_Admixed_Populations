
for i in {1..22};
do
bsub -M 20000 Rscript --vanilla ~/height_prediction/gwas/WHI/scripts/run_AS_GWAS.R ${i}
done


for i in 100 200 1000 2000 10000 20000 100000 200000;
do
bsub -M 40000 -o logLocAnc${i} -e logLocAnc${i} Rscript --vanilla ~/height_prediction/gwas/WHI/scripts/plot_diff_map.R phys_100000_0.0005 ${i}
echo ${i}
done

for i in {1..22}; do bsub -M 80000 Rscript --vanilla ~/height_prediction/gwas/WHI/scripts/make_anc.R ${i}; done
