#!/bin/bash
#SBATCH -J partvar
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --array=0-5
#SBATCH -o partvar_%a.log
# Variance-route partition simulation on PanUKB_EUR_ld blocks 1-100 with the real annotations (run from
# /home/yxy1234/plink4r_test/partition_sim; plink4r installed per task (below); scripts of sim/partition here). One array task per
# design and seed (2 designs x 3 seeds; task = 3 * design + seed - 1), all in parallel. Two designs (h2 shares;
# intercept = 20%, t(3) effects on the base causal SNPs outside every design category); PI = 0.01 causal (at least MIN_CAUSAL = 100 per category),
# 6:3:1 slab mixture, N = 300000, h2 = 0.3.
# Methods ORACLE, SSL, MASH, SSL_PI, MASH_PI (run_var.R; eigen cutoff 1), directory var_<design>_s<seed>. GIBBS_DRAWS (env, default 300), PI_BY (env,
# default elbo; pseudo_r: directory <dir>_pr) and GRID_EXTEND (env, default 4) are passed on to run_var.R.
# Submit: sbatch sim_var.sh
# Shared across tasks and cases (made once: the ref table and LD bins under a lock, annot.tsv and snps.tsv renamed
# into place): ref_1-100.tsv, ldbin0.98_1-100.bin, annot.tsv, snps.tsv.
# Each directory keeps only gen.txt, run.txt and eval.txt; a rerun overwrites.
set -euo pipefail
LD=/mnt/vstor/SOM_EPBI_XXZ10/yxy1234/LDpanel/PanUKB_EUR_ld
AN=/mnt/vstor/SOM_EPBI_XXZ10/yxy1234/LDpanel/Annotation
cd /home/yxy1234/plink4r_test/partition_sim
# plink4r installed per task on its own node (R compiles with -march=native: a build runs only on that CPU type),
# into lib_<job id>, removed at the end
LIB=/home/yxy1234/plink4r_test/lib_${SLURM_JOB_ID}; mkdir -p $LIB; trap 'rm -rf $LIB' EXIT
# a clean copy of the source (git archive: no object files of another node, no clash between parallel tasks)
mkdir -p $LIB/pkg; git -C /home/yxy1234/plink4r_test/src archive HEAD | tar -x -C $LIB/pkg
R CMD INSTALL -l $LIB $LIB/pkg > /dev/null; export R_LIBS=$LIB
export SHARED=$PWD
DS=(lipid bp); d=${DS[$((SLURM_ARRAY_TASK_ID / 3))]}; seed=$((SLURM_ARRAY_TASK_ID % 3 + 1))
CODING=High+Moderate+Low+Splice_Region+UTR_Core+UTR5_Other+UTR3_Other
declare -A DES
DES[lipid]="coding=$CODING:0.5,hepatocyte=Hepatocyte_Cholangiocyte:0.3"
DES[bp]="coding=$CODING:0.2,endothelial=Endothelial:0.3,smc=Stromal_Pericyte_SMC:0.3"
D=var_${d}_s$seed; M="ORACLE SSL MASH SSL_PI MASH_PI"
[ "${PI_BY:-elbo}" = pseudo_r ] && D=${D}_pr   # pi chosen by pseudo r: its own directory
BINS=ldbin0.98_1-100.bin
rm -rf $D; mkdir -p $D
flock shared.lock bash -c "[ -f ref_1-100.tsv ] || python3 make_ref.py $LD 1-100 ref_1-100.tsv"
ANNOT=$AN/annot_9M_all_20261122.parquet DESIGN="${DES[$d]}" PI=0.01 NGWAS=300000 BLOCKS=1-100 \
  python3 gen_var.py $LD ref_1-100.tsv $D $seed > $D/gen.txt
run() { METHODS="${M// /,}" LDSCORE=$AN/ldscore_9M_all_20261122.parquet BINS=$BINS \
  Rscript run_var.R $LD ref_1-100.tsv $D > $D/run.txt; }
if [ -f $BINS ]; then run; else flock shared.lock bash -c "$(declare -p M D BINS LD AN); $(declare -f run); run"; fi
python3 eval_var.py $LD $D $M > $D/eval.txt
find $D -type f ! -name gen.txt ! -name run.txt ! -name eval.txt -delete
