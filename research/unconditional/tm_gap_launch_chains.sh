#!/bin/bash
cd "D:/Dev/math/research/unconditional"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY="C:/Users/hooya/AppData/Local/Programs/Python/Python312/python.exe"
mkdir -p tm_gap_data/chains tm_gap_logs/chains
launch() { # L beta nconf seed levels
  local L=$1 beta=$2 nc=$3 seed=$4 lv=$5
  local out="tm_gap_data/chains/L${L}_b${beta}_s${seed}.npz"
  local log="tm_gap_logs/chains/L${L}_b${beta}_s${seed}.log"
  $PY tm_gap_run.py --L $L --beta $beta --ntherm 1000 --gapevery 4 \
      --nconf $nc --levels $lv --seed $seed --out $out > $log 2>&1 &
}
# confined points, many independent chains
for s in 201 202 203 204 205 206 207 208; do launch 8 2.2 4000 $s "0,4,8,12"; done
for s in 221 222 223 224 225 226;         do launch 8 2.3 4000 $s "0,4,8,12"; done
for s in 241 242 243 244;                 do launch 6 2.2 8000 $s "0,4,8";    done
for s in 261 262 263 264;                 do launch 6 2.3 8000 $s "0,4,8";    done
wait
echo "ALL-CHAINS DONE $(date)"
