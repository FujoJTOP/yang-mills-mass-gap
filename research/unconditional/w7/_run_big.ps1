$ErrorActionPreference='Continue'
cd D:\Dev\math\research\unconditional\w7
"=== L=12 (therm1500 + 9000 sweeps, every 2) ==="
py -3.12 -I block_mc.py --L 12 --beta 2.5 --therm 1500 --sweeps 9000 --every 2 --seed 2031 --out bmc_L12_b2.500.npz 2>&1 | Select-Object -Last 8
"=== L=16 (therm1200 + 5400 sweeps, every 2) ==="
py -3.12 -I block_mc.py --L 16 --beta 2.5 --therm 1200 --sweeps 5400 --every 2 --seed 2032 --out bmc_L16_b2.500.npz 2>&1 | Select-Object -Last 8
"ALL DONE"
