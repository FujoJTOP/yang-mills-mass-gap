import subprocess, sys, os, time
env = dict(os.environ)
env['NUMBA_NUM_THREADS']='1'; env['OPENBLAS_NUM_THREADS']='1'; env['OMP_NUM_THREADS']='1'
for b, seed in [(3.6,421),(4.0,422),(4.4,423),(5.0,424)]:
    out=f'runs/L6_b{b:.3f}_s{seed}.npz'
    if os.path.exists(out): print('skip',out,flush=True); continue
    t=time.time(); print(f'=== scan beta={b}',flush=True)
    r=subprocess.run([sys.executable,'probe.py','--L','6','--beta',str(b),
                      '--sweeps','1000','--therm','800','--meas','2','--seed',str(seed),
                      '--out',out],env=env)
    print(f'=== done rc={r.returncode} {time.time()-t:.0f}s',flush=True)
print('SCAN DONE')
