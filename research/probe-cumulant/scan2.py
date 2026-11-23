import subprocess, sys, os, time
env=dict(os.environ); env['NUMBA_NUM_THREADS']='1'; env['OPENBLAS_NUM_THREADS']='1'; env['OMP_NUM_THREADS']='1'
for b,seed in [(4.7,425),(5.3,426),(5.6,427)]:
    out=f'runs/L6_b{b:.3f}_s{seed}.npz'
    if os.path.exists(out): print('skip',out,flush=True); continue
    t=time.time(); print('=== beta',b,flush=True)
    r=subprocess.run([sys.executable,'probe.py','--L','6','--beta',str(b),'--sweeps','1000',
                      '--therm','800','--meas','2','--seed',str(seed),'--out',out],env=env)
    print('=== done',r.returncode,f'{time.time()-t:.0f}s',flush=True)
print('SCAN2 DONE')
