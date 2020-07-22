for PHYS in 5000 10000 25000 50000 75000 100000 500000 1000000;
do
for P in 0.0005 0.00005 0.000005 0.0000005 0.00000005;
do
bsub -M 40000 Rscript --vanilla ~/height_prediction/scripts/density_1000G_ukb.R phys_${PHYS}_${P}
done
echo ${PHYS}_${P}
done

for GENET in 0.1 0.15 0.2 0.25 0.3 0.5 1
do
for P in 0.0005 0.00005 0.000005 0.0000005 0.00000005;
do
bsub -M 40000 Rscript --vanilla ~/height_prediction/scripts/density_1000G_ukb.R genet_${GENET}_${P}
done
echo ${GENET}_${P}
done
