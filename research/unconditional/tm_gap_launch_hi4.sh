#!/bin/bash
cd "D:/Dev/math/research/unconditional"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="C:/Users/hooya/AppData/Local/Programs/Python/Python312/python.exe"
run() {
  local L=$1 beta=$2 nc=$3 seed=$4
  echo "START-hi4 L=$L beta=$beta nconf=$nc seed=$seed $(date)"
  $PY tm_gap_run.py --L $L --beta $beta --ntherm 800 --gapevery 4 \
      --nconf $nc --levels 0,8,16 --seed $seed \
      --out tm_gap_data/L${L}_b${beta}_hi.npz > tm_gap_logs/L${L}_b${beta}_hi.log 2>&1
  echo "DONE-hi4 L=$L beta=$beta $(date) rc=$?"
}
run 4 2.2 8000 41 &
run 4 2.5 8000 42 &
run 4 2.8 8000 43 &
wait
echo "ALL-HI4 DONE $(date)"
