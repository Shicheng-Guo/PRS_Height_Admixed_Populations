#!/bin/bash

home='~/height_prediction'
eval home=$home

dtset1=$1
dtset2=$2

for i in 250000 100000 50000;
do
echo 'Window size is '
echo $i
bsub -M 90240 -o ~/height_prediction/logs/log_${dtset1}_${dtset2}_${i}_${j}_prs -e ~/height_prediction/logs/log_${dtset1}_${dtset2}_${i}_${j}_prs  Rscript --vanilla ~/height_prediction/scripts/run_PRS.R ${dtset1} ${dtset2} LD_${i}_0.01_0.5
echo $i
echo 'done'
done

bsub -M 90240 -o ~/height_prediction/logs/log_${dtset1}_${dtset2}_LD_block_AFR_prs -e ~/height_prediction/logs/log_${dtset1}_${dtset2}_LD_block_AFR_prs Rscript --vanilla ~/height_prediction/scripts/run_PRS.R ${dtset1} ${dtset2} LD_block_0_0_AFR


bsub -M 90240 -o ~/height_prediction/logs/log_${dtset1}_${dtset2}_LD_block_EUR_prs -e ~/height_prediction/logs/log_${dtset1}_${dtset2}_LD_block_EUR_prs  Rscript --vanilla ~/height_prediction/scripts/run_PRS.R ${dtset1} ${dtset2} LD_block_0_0_EUR


for i in 1000000 500000 100000 75000 50000 25000 10000 5000;
do
echo 'Window size is '
echo $i
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005;
do
echo $j
bsub -M 90240 -o ~/height_prediction/logs/log_${dtset1}_${dtset2}_${i}_${j}_prs  -e ~/height_prediction/logs/log_${dtset1}_${dtset2}_${i}_${j}_prs Rscript --vanilla ~/height_prediction/scripts/run_PRS.R ${dtset1} ${dtset2} phys_${i}_${j}
echo $j
echo 'done'
done
echo $i
echo 'done'
done

for i in 1 0.5 0.3 0.25 0.2 0.15 0.1;
do
echo 'Window size is '
echo $i
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005;
do
echo $j
bsub -M 90240 -o ~/height_prediction/logs/log_${i}_${j}_prs -e ~/height_prediction/logs/log_${i}_${j}_prs  Rscript --vanilla ~/height_prediction/scripts/run_PRS.R ${dtset1} ${dtset2} genet_${i}_${j}
echo $j
echo 'done'
done
echo $i
echo 'done'
done

