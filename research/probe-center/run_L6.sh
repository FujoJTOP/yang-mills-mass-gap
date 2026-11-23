#!/bin/bash
PY=/c/Users/hooya/AppData/Local/Programs/Python/Python312/python
B=/c/Users/hooya/AppData/Local/Programs/Python/Python312/python
cd "D:/Dev/math/research/probe-center"
run() { "$PY" run_probe.py --npz "$1" --label "$2" --rmax "$3"; }
# primary grid L=6
run "D:/Dev/math/research/lattice-cpb/runs/L6_b2.600_s51.npz" b2.6_L6 3
run "D:/Dev/math/research/lattice-cpb/runs/L6_b2.800_s52.npz" b2.8_L6 3
run "D:/Dev/math/research/lattice-cpb/runs/L6_b3.000_s53.npz" b3.0_L6 3
run "D:/Dev/math/research/lattice-cpb/runs/L6_b3.300_s54.npz" b3.3_L6 3
# strong-coupling conservative
run "D:/Dev/math/research/lattice-cpb/runs/L6_b1.700_s12.npz" b1.7_L6 3
run "D:/Dev/math/research/lattice-cpb/runs/L6_b2.200_s22.npz" b2.2_L6 3
# weak control (scaling region)
run "D:/Dev/math/research/lattice-cpb/runs/L6_b5.700_s31.npz" b5.7_L6 3
run "D:/Dev/math/research/lattice-cpb/runs/L6_b6.000_s32.npz" b6.0_L6 3
