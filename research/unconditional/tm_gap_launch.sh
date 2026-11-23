#!/bin/bash
cd "D:/Dev/math/research/unconditional"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="C:/Users/hooya/AppData/Local/Programs/Python/Python312/python.exe"
# L=4 jobs: nconf=20000 ; L=6 jobs: nconf=6000
run() { # L beta nconf seed
  local L=$1 beta=$2 nc=$3 seed=$4
  echo "START L=$L beta=$beta nconf=$nc seed=$seed $(date)"
  $PY tm_gap_run.py --L $L --beta $beta --ntherm 1000 --gapevery 4 \
      --nconf $nc --levels 0,4,8 --seed $seed \
      --out tm_gap_data/L${L}_b${beta}.npz > tm_gap_logs/L${L}_b${beta}.log 2>&1
  echo "DONE  L=$L beta=$beta $(date) rc=$?"
}
run 4 2.2 20000 11 &
run 4 2.5 20000 12 &
run 4 2.8 20000 13 &
run 6 2.2 6000 21 &
run 6 2.5 6000 22 &
run 6 2.8 6000 23 &
wait
echo "ALL DONE $(date)"
