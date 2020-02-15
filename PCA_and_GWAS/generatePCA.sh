#!/bin/sh
## Generate first 10 principal components from VCFs of modern population and
#` samples to be projected.
#' Directory should have panel.txt, add_pop_labels.R, list of populations and
#` VCF files required.

source $1 # arguments in template of generatePCA_arg.txt

cd $DIR
VCF=${VCF_FILENAME%.vcf.gz*}

module load python/2.7.9
echo "python ~/packages/gdc/vcf2eigenstrat.py -v $VCF_FILENAME -o $VCF"
python ~/packages/gdc/vcf2eigenstrat.py -v $VCF_FILENAME -o $VCF

module load R/3.2.1
Rscript --vanilla add_pop_labels.R $DIR/$VCF $DIR/$PANEL_FILENAME

#write smartpca.par file
printf D1:\\t%s\\nindivname:\\tD1/%s_renamed_pops.ind\\nsnpname:\\tD1/%s\
.snp\\ngenotypename:\\tD1/%s.geno\\nnumoutliter:\\t0\\nnumoutevec:\
\\t10\\nevecoutname:\\tD1/%s.evec\\nphylipname:\\tD1/pca_fst.phyl\
\\npoplistname:\\tD1/%s\\nlsqproject:\\tYES \
$DIR $VCF $VCF $VCF $VCF $POPLIST_FILENAME > $DIR/smartpca.par

smartpca -p smartpca.par
