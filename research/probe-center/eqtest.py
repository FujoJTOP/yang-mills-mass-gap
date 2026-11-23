import numpy as np, sys, time
sys.path.insert(0,'D:/Dev/math/research/lattice-cpb')
import su3_numba as SN
# 1) P of stored snapshots
for f in ['L6_b2.200_s22','L8_b2.200_s23','L6_b2.600_s51','L6_b3.000_s53','L8_b6.000_s33']:
    d=np.load(f'D:/Dev/math/research/lattice-cpb/runs/{f}.npz',allow_pickle=True)
    L=int(d['L']); U=d['snapshots'][-1]
    print(f, 'L',L,'beta',float(d['beta']),'P(snapshot last)=',SN.plaquette_mean(U,L),'P(recorded)=',d['plaq'].mean())
# 2) hot start MC at L8 beta=2.6
L=8; beta=2.6; kappa=beta/6
rng=np.random.default_rng(5)
U=np.zeros((4,L,L,L,L,3,3),np.complex128)
a=rng.normal(size=(4,L,L,L,L,3,3))+1j*rng.normal(size=(4,L,L,L,L,3,3))
q,r=np.linalg.qr(a); dg=np.einsum('...ii->...i',r); dg/=np.abs(dg); q*=dg[...,None,:]
det=np.linalg.det(q); U=q/(det**(1/3))[...,None,None]
print("\nHOT start L8 b=2.6:")
t0=time.time()
for k in range(2000):
    SN.sweep_par(U,L,kappa)
    if (k+1)%200==0: print(f"  {k+1}: P={SN.plaquette_mean(U,L):.5f} t={time.time()-t0:.0f}s",flush=True)
