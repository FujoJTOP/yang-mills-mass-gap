#!/bin/bash
PY=/c/Users/hooya/AppData/Local/Programs/Python/Python312/python
cd "D:/Dev/math/research/probe-center"
run() { "$PY" run_probe.py --npz "$1" --label "$2" --rmax "$3"; }
R6="D:/Dev/math/research/lattice-cpb/runs"
run "$R6/L6_b1.700_s12.npz" b1.7_L6 3
run "$R6/L6_b2.200_s22.npz" b2.2_L6 3
run "$R6/L6_b2.600_s51.npz" b2.6_L6 3
run "$R6/L6_b2.800_s52.npz" b2.8_L6 3
run "$R6/L6_b3.000_s53.npz" b3.0_L6 3
run "$R6/L6_b3.300_s54.npz" b3.3_L6 3
run "$R6/L6_b5.700_s31.npz" b5.7_L6 3
run "$R6/L6_b6.000_s32.npz" b6.0_L6 3
run "$R6/L8_b1.700_s13.npz" b1.7_L8 3
run "$R6/L8_b2.200_s23.npz" b2.2_L8 3
run "runs_L8_b2.600.npz" b2.6_L8 4
run "runs_L8_b3.000.npz" b3.0_L8 4
run "$R6/L8_b6.000_s33.npz" b6.0_L8 4
run "$R6/L10_b6.000_s41.npz" b6.0_L10 4
echo "ALL PROBES DONE"
