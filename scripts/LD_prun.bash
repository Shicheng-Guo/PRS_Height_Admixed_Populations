#!/bin/sh

dtset1=$1
dtset2=$2
echo $dtset1
echo $dtset2

for i in {1..22};
do
echo 'Chr is '
echo $i
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005; 
do
echo $j
for w in 5000 10000 25000 50000 75000 100000 500000 1000000;
do 
#bsub  -M 20240 -o $parentdir/logs/logphys${j}_${i}_${w} Rscript --vanilla $currentdir/LD_prunning.R --chromosome $i --method phys --windowsize $w --pvalue $j
Rscript --vanilla LD_prunning.R --chromosome $i --method phys --windowsize $w --pvalue $j --dataset1 $dtset1 --dataset2 $dtset2
echo $w
echo 'done'
done
for w in 0.1 0.15 0.2 0.25 0.3 0.5 1;
do
#bsub  -M 20240 -o $parentdir/logs/loggenet${j}_${i}_${w} Rscript --vanilla $currentdir/LD_prunning.R --chromosome $i --method genet --windowsize w --pvalue $j
Rscript --vanilla LD_prunning.R --chromosome $i --method genet --windowsize $w --pvalue $j --dataset1 $dtset1 --dataset2 $dtset2
echo $w
echo 'done'
done
echo $j
echo 'done'
done
echo $i
done

#alicia martin  (peding r2 filter, to do)

for i in {1..22};
do
echo 'Chr is '
echo $i
for w in 250000 100000 50000;
do
#bsub -M 20240 -o $parentdir/logs/logmartin_${w} -e $parentdir/logs/logmartin_${w} Rscript --vanilla $currentir/LD_prunning.R --chromosome $i --method LD --windowsize $w --pvalue 0.01 --r2 0.5
Rscript --vanilla LD_prunning.R --chromosome $i --method LD --windowsize $w --pvalue 0.01 --r2 0.5 --dataset1 $dtset1 --dataset2 $dtset2
echo $w
echo 'done'
done
echo $i
echo '$done'
done

for i in {1..22};
do
echo 'Chr is '
echo $i
#bsub -M 20240 -o $parentdir/logs/logblock -e $parentdir/logs/logblock Rscript --vanilla $currentdir/LD_prunning.R --chromosome $i --method LD_block --windowsize 0 --pvalue 0
Rscript --vanilla LD_prunning.R --chromosome $i --method LD_block --windowsize 0 --pvalue 0 --dataset1 $dtset1 --dataset2 $dtset2
echo $j
echo 'done'
done

