#!/bin/sh

MYPWD=`dirname $0`
touch $MYPWD/run_this.sh
echo '#!/bin/bash' > $MYPWD/run_this.sh
for j in 0.00000005 0.0000005 0.000005 0.00005 0.0005;
do
for k in 1000000 500000 100000 75000 50000 25000 10000 5000;
do
#echo  "bsub -M 60000 Rscript --vanilla /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/combine_Rds.R phys_${k}_${j}" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this.sh
echo "Rscript --vanilla $MYPWD/combine_Rds.R phys_${k}_${j}" >> $MYPWD/run_this.sh
done
for l in 1 0.5  0.3 0.25 0.2 0.15 0.1;
do
#echo "bsub -M 60000 Rscript --vanilla /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/combine_Rds.R genet_${l}_${j}" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this.sh
echo "Rscript --vanilla $MYPWD/combine_Rds.R genet_${l}_${j}" >> $MYPWD/run_this.sh
done
done

#these three work differently. See run_PRS.R
for k in 250000 100000 50000;
do
#echo "bsub -M 60000 Rscript --vanilla /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/combine_Rds.R LD_${k}_0.01_0.5" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this.sh
echo "Rscript --vanilla $MYPWD/combine_Rds.R LD_${k}_0.01_0.5" >> $MYPWD/run_this.sh
done

#echo "bsub -M 60000 Rscript --vanilla /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/combine_Rds.R LD_block_0_0_AFR" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this.sh
#echo "bsub -M 60000 Rscript --vanilla /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/combine_Rds.R LD_block_0_0_EUR" >> /project/mathilab/bbita/gwas_admix/sib_gwas/JHS/run_this.sh


echo "Rscript --vanilla $MYPWD/combine_Rds.R LD_block_0_0_AFR" >> $MYPWD/run_this.sh
echo "Rscript --vanilla $MYPWD/combine_Rds.R LD_block_0_0_EUR" >> $MYPWD/run_this.sh

chmod 777 $MYPWD/run_this.sh

#echo $MYPWD
