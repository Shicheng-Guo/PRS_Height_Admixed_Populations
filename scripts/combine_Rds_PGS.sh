#!/bin/sh

home='~/height_prediction'
eval home=$home
MY=${home}/$1/$2

echo $MY



touch $MY/scripts/run_this_PGS.sh
touch $MY/scripts/test2.txt
echo '#!/bin/sh' > $MY/scripts/run_this_PGS.sh
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005;
do
for k in 1000000 500000 100000 75000 50000 25000 10000 5000;
do
echo  "phys_${k}_${j}" >> ${MY}/scripts/test2.txt
done
for l in 1 0.5 0.3 0.25 0.2 0.15 0.1;
do
echo "genet_${l}_${j}" >> ${MY}/scripts/test2.txt
done
done

echo "LD_250000_0.01_0.5" >> ${MY}/scripts/test2.txt
echo "LD_100000_0.01_0.5" >> ${MY}/scripts/test2.txt
echo "LD_50000_0.01_0.5" >> ${MY}/scripts/test2.txt
echo "LD_block_0_0_AFR" >> ${MY}/scripts/test2.txt
echo "LD_block_0_0_EUR" >> ${MY}/scripts/test2.txt


echo 'Rscript --vanilla ~/height_prediction/scripts/combine_Rds_PGS.R \' >> ${MY}/scripts/run_this_PGS.sh
echo ${1} ${2} '\'  >> ${MY}/scripts/run_this_PGS.sh
awk 'ORS=" "{print $1}'  ${MY}/scripts/test2.txt >> ${MY}/scripts/run_this_PGS.sh
chmod 777 ${MY}/scripts/run_this_PGS.sh

bash $MY/scripts/run_this_PGS.sh
