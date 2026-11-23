"""Addendum: evaluate T1 quantities at the SU(2) q~1 point (beta~6, xi~1; see STATUS.md)."""
import t1_kappa as T
import numpy as np

print("A'  character coefficients at beta = 6.0, 6.3 (q~1, xi~1):")
for beta in [6.0, 6.3]:
    vals=[T.c_r_beta(r,beta) for r in range(7)]
    print(f"  beta={beta}: c_r(r=0..6) = "+"  ".join(f"{v:.4f}" for v in vals)+f"   min={min(vals):.3e}")

print("C'  exact signed-vs-abs surface-weight enumeration (12-plaquette 4-valent 2-chain):")
for beta in [6.0, 6.3]:
    res=T.enumerate_chain(beta,L=2)
    print(f"  beta={beta}: terms={res['ntot']} neg={res['nneg']} min_w={res['minw']:.3e} "
          f"z_signed/z_abs={res['ratio']:.12f}  kappa_c={res['ratio']**(1.0/res['nplaquette']):.12f}")

print("D'  SU(2) checkerboard heat-bath <P> at beta=6.0 (L=6):")
m,eb,es=T.run_plaq(6.0,L=6,nwarm=400,nmeas=300,nbin=15)
print(f"  beta=6.0: <P> = {m:.4f} +/- {eb:.4f} (binned)")
