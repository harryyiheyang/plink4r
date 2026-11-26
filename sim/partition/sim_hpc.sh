#!/bin/bash
#SBATCH -J partsim
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH -o partsim_%j.log
# Partition simulation on PanUKB_EUR_ld blocks 1-100 (run from ~/partition_sim; plink4r in ./lib).
# Two settings: ENRICH=1 (tissue SNPs as likely causal as the others) and ENRICH=5; 3 seeds each.
set -euo pipefail
LD=/mnt/vstor/SOM_EPBI_XXZ10/yxy1234/LDpanel/PanUKB_EUR_ld
cd ~/partition_sim
[ -f ref_1-100.tsv ] || python3 make_ref.py $LD 1-100 ref_1-100.tsv
for en in 1 5; do
  for seed in 1 2 3; do
    D=sim2_e${en}_s$seed
    ENRICH=$en NGWAS=300000 BLOCKS=1-100 python3 gen_ideal.py $LD ref_1-100.tsv $D $seed
    Rscript run_ideal.R $LD ref_1-100.tsv $D | tee $D/run.txt
    python3 eval_ideal.py $LD $D $D/V.tsv $D/A1.tsv $D/AM.tsv $D/AML.tsv | tee $D/eval.txt
  done
done
