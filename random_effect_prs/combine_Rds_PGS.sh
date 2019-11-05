#!/bin/sh

home='~/height_prediction'
eval home=$home
MY=${home}/random_effect_prs/$2

echo $MY



touch $MY/run_this_PGS.sh
touch $MY/test2.txt
echo '#!/bin/sh' > $MY/run_this_PGS.sh
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005;
do
for k in 1000000 500000 100000 75000 50000 25000 10000 5000;
do
echo  "phys_${k}_${j}" >> ${MY}/test2.txt
done
done

echo 'Rscript --vanilla ~/height_prediction/random_effect_prs/combine_Rds_PGS.R \' >> ${MY}/run_this_PGS.sh
echo ${1} ${2} '\'  >> ${MY}/run_this_PGS.sh
awk 'ORS=" "{print $1}'  ${MY}/test2.txt >> ${MY}/run_this_PGS.sh
chmod 777 ${MY}/run_this_PGS.sh

bash $MY/run_this_PGS.sh
