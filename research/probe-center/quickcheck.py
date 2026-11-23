import numpy as np, sys
sys.path.insert(0,'D:/Dev/math/research/lattice-cpb')
import su3_numba as SN
def hot(L): 
    rng=np.random.default_rng(3); a=rng.normal(size=(4,L,L,L,L,3,3))+1j*rng.normal(size=(4,L,L,L,L,3,3))
    q,r=np.linalg.qr(a); d=np.einsum('...ii->...i',r); d/=np.abs(d); q*=d[...,None,:]
    det=np.linalg.det(q); return q/(det**(1/3))[...,None,None]
for L in [6]:
    U=hot(L); kap=2.6/3.0
    for k in range(600): SN.sweep_par(U,L,kap)
    print(f"L={L} kappa=beta/3=2.6/3: P={SN.plaquette_mean(U,L):.5f}  (existing L6 config=0.172)")
