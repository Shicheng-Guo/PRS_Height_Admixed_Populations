#!/bin/bash

rm /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
sleep 3
touch /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005;
do
for k in 1000000 500000 100000 75000 50000 25000 10000 5000;
do
echo  "phys_${k}_${j}" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
done
for l in 1 0.5 0.3 0.25 0.2 0.15 0.1;
do
echo "genet_${l}_${j}" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
done
done

echo "LD_250000_0.01_0.5" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
echo "LD_100000_0.01_0.5" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
echo "LD_50000_0.01_0.5" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
echo "LD_block_0_0_AFR" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt
echo "LD_block_0_0_EUR" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt

rm /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this_PGS.sh
sleep 3

touch /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this_PGS.sh
echo '#!/bin/bash' > /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this_PGS.sh
echo 'Rscript --vanilla /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/combine_Rds_PGS.R \' >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this_PGS.sh
awk 'ORS=" "{print $1}' /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/test2.txt >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this_PGS.sh
chmod 777 /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this_PGS.sh

