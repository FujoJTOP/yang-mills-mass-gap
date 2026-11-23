import subprocess, sys, os, time
env = dict(os.environ)
env['NUMBA_NUM_THREADS'] = '1'
env['OPENBLAS_NUM_THREADS'] = '1'
env['OMP_NUM_THREADS'] = '1'
jobs = [
    (6, 2.2, 3000, 1500, 401),
    (6, 2.6, 3000, 1500, 402),
    (6, 3.0, 3000, 1500, 404),
    (6, 3.3, 3000, 1500, 405),
    (6, 6.0, 3000, 1200, 406),
    (8, 3.0, 3000, 1200, 412),
]
for L, b, sw, th, seed in jobs:
    out = f'runs/L{L}_b{b:.3f}_s{seed}.npz'
    if os.path.exists(out):
        print('skip', out, flush=True); continue
    t = time.time()
    print(f'=== L={L} beta={b} sweeps={sw} therm={th} seed={seed}', flush=True)
    r = subprocess.run([sys.executable, 'probe.py', '--L', str(L), '--beta', str(b),
                        '--sweeps', str(sw), '--therm', str(th), '--meas', '2',
                        '--seed', str(seed), '--out', out], env=env)
    print(f'=== done rc={r.returncode} {time.time()-t:.0f}s', flush=True)
print('ALL DONE')
