import numpy as np, sys, time
sys.path.insert(0,'D:/Dev/math/research/lattice-cpb')
import su3_numba as SN
import su3 as S

def hot(L, beta, seed):
    rng=np.random.default_rng(seed)
    a=rng.normal(size=(4,L,L,L,L,3,3))+1j*rng.normal(size=(4,L,L,L,L,3,3))
    q,r=np.linalg.qr(a); dg=np.einsum('...ii->...i',r); dg/=np.abs(dg); q*=dg[...,None,:]
    det=np.linalg.det(q); return q/(det**(1/3))[...,None,None]

beta=2.6; kappa=beta/6
for L in [6,8]:
    U=hot(L,beta,11)
    print(f"\nNUMBA par L={L} b={beta}")
    for k in range(1500):
        SN.sweep_par(U,L,kappa)
        if (k+1)%300==0: print(f"  {k+1}: P={SN.plaquette_mean(U,L):.5f}",flush=True)
# sequential numba at L=6
U=hot(6,beta,12)
print(f"\nNUMBA seq L=6 b={beta}")
for k in range(800):
    SN.sweep(U,6,kappa)
    if (k+1)%200==0: print(f"  {k+1}: P={SN.plaquette_mean(U,6):.5f}",flush=True)
