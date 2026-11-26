#!/bin/bash
#SBATCH -J partvar
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o partvar_%j.log
# Variance-route partition simulation on PanUKB_EUR_ld blocks 1-100 with the real annotations (run from
# ~/partition_sim2; plink4r in ./lib; scripts of sim/partition here). Two designs (h2 shares; residual = unannotated
# SNPs, 20%, t(3) effects), 3 seeds each; PI = 0.01 causal, 6:3:1 slab mixture, N = 300000, h2 = 0.3.
set -euo pipefail
LD=/mnt/vstor/SOM_EPBI_XXZ10/yxy1234/LDpanel/PanUKB_EUR_ld
AN=/mnt/vstor/SOM_EPBI_XXZ10/yxy1234/LDpanel/Annotation
cd ~/partition_sim2
[ -f ref_1-100.tsv ] || python3 make_ref.py $LD 1-100 ref_1-100.tsv
CODING=High+Moderate+Low+Splice_Region+UTR_Core+UTR5_Other+UTR3_Other
declare -A DES
DES[lipid]="coding=$CODING:0.5,hepatocyte=Hepatocyte_Cholangiocyte:0.3"
DES[bp]="coding=$CODING:0.2,endothelial=Endothelial:0.3,smc=Stromal_Pericyte_SMC:0.3"
for d in lipid bp; do
  for seed in 1 2 3; do
    D=var_${d}_s$seed
    ANNOT=$AN/annot_9M_all_20261122.parquet DESIGN="${DES[$d]}" PI=0.01 NGWAS=300000 BLOCKS=1-100 \
      python3 gen_var.py $LD ref_1-100.tsv $D $seed | tee $D.gen.txt
    LDSCORE=$AN/ldscore_9M_all_20261122.parquet BINS=ldbin0.98_1-100.bin Rscript run_var.R $LD ref_1-100.tsv $D | tee $D/run.txt
    python3 eval_var.py $LD $D V A1 AMm AM AML | tee $D/eval.txt
  done
done
