$ErrorActionPreference='Continue'
cd D:\Dev\math\research\unconditional\w7
"=== L=6  生产跑（只存相关函数）==="
py -3.12 -I block_mc.py --L 6 --beta 2.5 --therm 3000 --sweeps 20000 --every 1 --seed 2026 --out bmc_L6_b2.500.npz 2>&1 | Select-Object -Last 9
"=== L=8  生产跑 ==="
py -3.12 -I block_mc.py --L 8 --beta 2.5 --therm 3000 --sweeps 14000 --every 2 --seed 2027 --out bmc_L8_b2.500.npz 2>&1 | Select-Object -Last 9
"=== L=10 生产跑 ==="
py -3.12 -I block_mc.py --L 10 --beta 2.5 --therm 2500 --sweeps 11000 --every 2 --seed 2028 --out bmc_L10_b2.500.npz 2>&1 | Select-Object -Last 9
"=== L=6  RTC 跑（存 q 场）==="
py -3.12 -I block_mc.py --L 6 --beta 2.5 --therm 3000 --sweeps 16000 --every 4 --seed 2029 --store-fields --out bmc_L6_b2.500_fields.npz 2>&1 | Select-Object -Last 9
"ALL DONE"
