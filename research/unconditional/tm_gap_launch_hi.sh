#!/bin/bash
cd "D:/Dev/math/research/unconditional"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="C:/Users/hooya/AppData/Local/Programs/Python/Python312/python.exe"
run() {
  local L=$1 beta=$2 nc=$3 seed=$4
  echo "START-hi L=$L beta=$beta nconf=$nc seed=$seed $(date)"
  $PY tm_gap_run.py --L $L --beta $beta --ntherm 800 --gapevery 4 \
      --nconf $nc --levels 0,8,16 --seed $seed \
      --out tm_gap_data/L${L}_b${beta}_hi.npz > tm_gap_logs/L${L}_b${beta}_hi.log 2>&1
  echo "DONE-hi L=$L beta=$beta $(date) rc=$?"
}
run 6 2.2 3500 31 &
run 6 2.5 3500 32 &
run 6 2.8 3500 33 &
wait
echo "ALL-HI DONE $(date)"
