#!/bin/bash

home='~/height_prediction'
eval home=$home

dtset1=$1
dtset2=$2

for i in 1000000 500000 100000 75000 50000 25000 10000 5000;
do
echo 'Window size is '
echo $i
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005;
do
echo $j
bsub -M 90240 -o ~/height_prediction/unweighted_prs/logs/log_${dtset1}_${dtset2}_${i}_${j}_prs  -e ~/height_prediction/unweighted_prs/logs/log_${dtset1}_${dtset2}_${i}_${j}_prs Rscript --vanilla ~/height_prediction/unweighted_prs/run_PRS.R ${dtset1} ${dtset2} phys_${i}_${j}
echo $j
echo 'done'
done
echo $i
echo 'done'
done

