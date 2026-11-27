"""GP9 -- per-block correlator accumulator: tiny output, so the 10x-statistics run needs no big npz."""
import numpy as np, json, os, sys
OBS=[0,3,4]; NB=20
d=np.load('gh3_probe.npz'); S,Lt,N=d['S'],int(d['Lt']),int(d['n'])
A=np.ascontiguousarray(S[:,0][:,OBS,:]).astype(np.float64); del S
Abar=A.mean(axis=(0,2))
X=A-Abar[None,:,None]                      # connected
blocks=np.array_split(np.arange(N),NB)
acc=np.zeros((NB,Lt,3,3))
for b,idx in enumerate(blocks):
    Xb=X[idx]
    for t in range(Lt):
        acc[b,t]=np.einsum('nct,ndt->ncd',Xb,np.roll(Xb,-t,axis=2),optimize=True)/Lt
np.save('gp9_block_acc.npy',acc)
np.save('gp9_vac.npy',Abar*np.sqrt(N))     # store scaled by sqrt(N) so sums combine
print(f'chunk N={N} Lt={Lt} -> gp9_block_acc.npy  sum C00={acc[:,0,0,0].sum():.6f}  plaq={float(d["plaq"].mean()):.6f}')
