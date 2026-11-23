#!/bin/bash
PY=/c/Users/hooya/AppData/Local/Programs/Python/Python312/python
cd "D:/Dev/math/research/probe-center"
run() { "$PY" run_probe.py --npz "$1" --label "$2" --rmax "$3"; }
run "D:/Dev/math/research/lattice-cpb/runs/L8_b1.700_s13.npz" b1.7_L8 3
run "D:/Dev/math/research/lattice-cpb/runs/L8_b2.200_s23.npz" b2.2_L8 3
run "D:/Dev/math/research/lattice-cpb/runs/L8_b6.000_s33.npz" b6.0_L8 4
run "D:/Dev/math/research/lattice-cpb/runs/L10_b6.000_s41.npz" b6.0_L10 4
