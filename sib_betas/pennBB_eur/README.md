*Prune using different methods*
```
cd scripts/
bash LD_prun.bash
bash combine_Rds_v2.sh
bash run_this.sh
cd ..
```
*Run polygenic scores*
```
cd scripts/
bash calc_PGS.sh
bash combine_Rds_PGS.sh
bash run_this_PGS.sh
Rscript --vanilla $cur/Plots_JHS.R
cd ..
```

