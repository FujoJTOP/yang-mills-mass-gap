"""GP9 parallel accumulator: sums per-block correlators over the 8 procs of one wave."""
import numpy as np, glob, sys, os
w = sys.argv[1]; OBS=[0,3,4]; NB=20; Lt=20
tot = np.zeros((NB,Lt,3,3)); vac = 0.0; nsum = 0; plaqs=[]
for k,f in enumerate(sorted(glob.glob(f'gh3_x{w}_*.npz'))):
    d=np.load(f); S=d['S']; N=int(d['n'])
    A=np.ascontiguousarray(S[:,0][:,OBS,:]).astype(np.float64); del S
    Ab=A.mean(axis=(0,2))
    X=A-Ab[None,:,None]; del A
    blocks=np.array_split(np.arange(N),NB)
    for b,idx in enumerate(blocks):
        Xb=X[idx]
        for t in range(Lt):
            tot[b,t]+=np.einsum('nct,ndt->ncd',Xb,np.roll(Xb,-t,axis=2),optimize=True)/Lt
    vac=Ab*np.sqrt(N); nsum+=N; plaqs.append(float(d['plaq'].mean()))
    print(f'  {os.path.basename(f)} N={N} plaq={plaqs[-1]:.6f}', flush=True)
np.save(f'gp9_par_w{w}.npy', tot); np.save(f'gp9_par_vac_w{w}.npy', vac)
print(f'wave {w}: N_total={nsum}  sum C00={tot[:,0,0,0].sum():.6f}  plaq={np.mean(plaqs):.6f}')
