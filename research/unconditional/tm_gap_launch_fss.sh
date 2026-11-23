#!/bin/bash
cd "D:/Dev/math/research/unconditional"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY="C:/Users/hooya/AppData/Local/Programs/Python/Python312/python.exe"
run() { # L beta nconf seed levels outtag
  local L=$1 beta=$2 nc=$3 seed=$4 lv=$5 tag=$6
  echo "START $tag L=$L beta=$beta nconf=$nc $(date)"
  $PY tm_gap_run.py --L $L --beta $beta --ntherm 1000 --gapevery 4 \
      --nconf $nc --levels $lv --seed $seed \
      --out tm_gap_data/${tag}.npz > tm_gap_logs/${tag}.log 2>&1
  echo "DONE  $tag $(date) rc=$?"
}
run 4 2.3 15000 51 "0,4,8"        L4_b2.3 &
run 6 2.3  8000 52 "0,4,8"        L6_b2.3 &
run 8 2.2  5000 61 "0,4,8,12"     L8_b2.2 &
run 8 2.3  5000 62 "0,4,8,12"     L8_b2.3 &
run 8 2.5  5000 63 "0,4,8,12"     L8_b2.5 &
run 6 2.2  8000 64 "0,4,8,12,16"  L6_b2.2_deep &
wait
echo "ALL-FSS DONE $(date)"
