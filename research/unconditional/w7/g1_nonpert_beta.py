"""
g1_nonpert_beta.py -- witness for G1-NONPERT-BETA.md  (py -3.12 -I, numpy only)

Shows:
 (1) pure-YM 2-loop beta < 0 for all g>0 (no perturbative zero) [T-small-g]
 (2) a-theorem allows BOTH a_IR=0 (gapped) and 0<a_IR<a_UV (CFT) -> no sign of beta
 (3) an artificial FULL beta = -b0 g^3 (1 - g^2/g*^2) has an IR-stable zero while its
     2-loop part -b0 g^3 has none  => 2-loop beta<0 cannot exclude a strong-coupling
     IR fixed point (2-loop beta<0 is NOT a nonperturbative sign theorem)
 (4) conformal window is REAL: SU(3) N_f=12 has b0>0, b1<0 => a Banks-Zaks zero g*>0
     (an AF theory that is gapless+conformal in the IR)  => 'AF => gapped' is false in general
"""
import numpy as np

def beta2(g, b0, b1):
    return -b0*g**3 - b1*g**5

def coeffs(N, Nf):
    b0 = (11*N - 2*Nf)/(48*np.pi**2)          # >0 for AF
    b1 = (34/3)*N**2 - (13/3)*N*Nf + Nf/N      # 2-loop (with 1/(16pi^2)^2 folded in convention)
    return b0, b1

print("="*70); print("(1) pure YM (N_f=0) 2-loop beta<0 for all g>0 (no positive zero)")
for N in (2,3,4):
    b0,b1 = coeffs(N,0)
    g = np.linspace(1e-3, 2.0, 4000)
    b = beta2(g,b0,b1)
    sign_changes = int(np.sum(np.diff(np.sign(b)) != 0))
    print(f"  SU({N}): b0={b0:.4f}>0 b1={b1:.4f}>0  max(beta)={b.max():.3e} (<0)  sign_changes={sign_changes}")

print("="*70); print("(2) a-theorem: a_UV >= a_IR allows BOTH gapped (a_IR=0) and CFT (a_IR>0)")
a_UV = (3**2-1)*31/180.0   # free-vector a per gauge boson (N=3)
for label,a_IR in [("gapped (a_IR=0)",0.0),("CFT (a_IR=0.5 a_UV)",0.5*a_UV),("CFT (a_IR=0.9 a_UV)",0.9*a_UV)]:
    print(f"  {label:24s} a_UV={a_UV:.4f} a_IR={a_IR:.4f}  satisfies a_UV>=a_IR? {a_UV>=a_IR-1e-12}")

print("="*70); print("(3) artificial full beta with IR zero vs its 2-loop part")
b0 = coeffs(3,0)[0]; gstar = 1.0   # illustrative
g = np.linspace(1e-3,2.0,4000)
beta_full = -b0*g**3*(1 - g**2/gstar**2)
beta_2loop = -b0*g**3
# full beta has a zero at g=gstar; 2-loop does not
idx = np.argmin(np.abs(g-gstar))
print(f"  at g*=1: beta_full={beta_full[idx]:.3e} (~0), beta_2loop={beta_2loop[idx]:.3e} (<0)")
print(f"  => a strong-coupling IR-stable zero coexists with 2-loop beta<0 everywhere")

print("="*70); print("(4) conformal window is REAL: SU(3) N_f=12 has Banks-Zaks zero")
for Nf in (0,12,16):
    b0,b1 = coeffs(3,Nf)
    has_bz = (b0>0 and b1<0)
    gstar2 = (-b0/b1) if has_bz else float('nan')
    print(f"  SU(3) N_f={Nf:2d}: b0={b0:.5f} b1={b1:.4f}  BZ-zero? {has_bz}"
          + (f"  g*={np.sqrt(gstar2):.4f}>0 (gapless conformal IR!)" if has_bz else "  (no zero)"))
print("  => 'AF => gapped' is FALSE in general; pure YM (N_f=0) gapped is DYNAMICAL, not kinematic.")
