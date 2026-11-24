"""STAR-NO-BZ witness: pure YM has no IR fixed point perturbatively; the conformal window shows
AF-UV + conformal-IR is a REAL phenomenon, so no-BZ for pure YM is dynamical (not kinematic)."""
import numpy as np

def beta_coeffs(N, Nf):
    CF = (N*N - 1)/(2*N)
    b0 = 11*N/3 - 2*Nf/3
    b1 = 34*N*N/3 - (10*N/3 + 2*CF)*Nf
    return b0, b1

def beta(g, b0, b1):
    # beta = -b0 g^3 - b1 g^5  (2-loop, sign convention AF = negative)
    return -b0*g**3 - b1*g**5

print("=== pure YM: beta(g) < 0 for ALL g>0 (perturbative, no zero) ===")
for N in (2,3,4):
    b0,b1 = beta_coeffs(N,0)
    g = np.linspace(1e-3, 2.0, 5)
    vals = beta(g,b0,b1)
    print(f" SU({N}) pure: b0={b0:.3f} b1={b1:.3f}  beta(g>0) all<0: {np.all(vals<0)}")
    print(f"   beta at g=0.3,1.0,2.0: {beta(0.3,b0,b1):.4g}, {beta(1.0,b0,b1):.4g}, {beta(2.0,b0,b1):.4g}")

print("\n=== conformal window: SU(N) with N_f near 11N/2 has a Banks-Zaks zero g*>0 ===")
for (N,Nf) in [(3,0),(3,8),(3,12),(3,16),(3,16.5)]:
    b0,b1 = beta_coeffs(N,Nf)
    has_bz = (b0>0 and b1<0)
    gstar = np.sqrt(-b0/b1) if has_bz else None
    print(f" SU({N}) N_f={Nf}: b0={b0:.3f} b1={b1:.3f} BZ_zero={has_bz}" + (f" g*={gstar:.3f}" if has_bz else " (beta<0 all g)"))

print("\n=== trace anomaly at a putative fixed point: beta(g*)=0 => <theta>=0 (CONSISTENT with CFT) ===")
N,Nf=3,12; b0,b1=beta_coeffs(N,Nf); gs=np.sqrt(-b0/b1)
print(f" at g*={gs:.3f}: beta(g*)={beta(gs,b0,b1):.2e} => <theta>=(beta/2g)<G^2> = 0  (no contradiction)")
