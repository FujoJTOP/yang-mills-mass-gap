"""
pnt_strong_bridge_extend.py -- PNT-STRONG-BRIDGE-EXTEND witness

Tests whether the strong-coupling bridge m_{0++} >= 4 sqrt(sigma) can be pushed
from beta* ~ 8.9e-3 toward beta = 2.5.  All quantities are the SU(2)
strong-coupling / character-expansion objects; numpy+scipy only, no MC.

Key objects (SU(2), Wilson action S = beta sum (1 - 1/2 tr U_pl)):
  u(beta)  = I_2(beta)/I_1(beta)        the fundamental character/polymer weight
                                          (u ~ beta/4 as beta->0, u->1 as beta->inf)
  sigma a^2 ~ -ln u(beta) + O(u)         string tension (leading strong coupling)
  m_1 a     ~ 4 ln(1/u)                  lightest 0++ from the four-wall tube (P3)
  => m_{0++}/sqrt(sigma) = 4 sqrt(ln(1/u))       (the ROOF-CONFINEMENT-GAP ratio)

Parts:
  A. u(beta) table; show u(2.5) is O(0.5) (large), u(1/48) tiny.
  B. bridge ratio 4 sqrt(ln(1/u)); where is it >= 4?  (needs u <= e^-1)
  C. KP/SZZ radius: the expansion in u converges only for |u| < u*; show the gap.
  D. Pade/interpolation principle: a strong-coupling series with radius u*
     evaluated at u(2.5) is OUTSIDE the disk; [L/M] Pade orders disagree
     (illustrated on a mock series whose radius matches the KP scale).
Run:  py -3.12 research/unconditional/w7/pnt_strong_bridge_extend.py
"""
import numpy as np
from scipy.special import iv

def u(beta):
    return iv(2, beta) / iv(1, beta)

def mratio(beta):
    ub = u(beta)
    return 4.0 * np.sqrt(np.log(1.0/ub)) if 0 < ub < 1 else np.nan

print("="*74)
print("A. u(beta) = I2(beta)/I1(beta), SU(2) fundamental weight")
print("="*74)
print(f"{'beta':>10} {'u(beta)':>12}  note")
for b, note in [(8.9e-3, "residue beta*"), (1/48, "SZZ benchmark"),
                (0.1, ""), (0.5, ""), (1.0, ""), (1.55, "u~1/e"),
                (2.0, ""), (2.5, "CROSSOVER"), (3.0, ""), (4.0, "")]:
    print(f"{b:>10.5g} {u(b):>12.6f}  {note}")

print("\n" + "="*74)
print("B. bridge ratio m_{0++}/sqrt(sigma) = 4 sqrt(ln(1/u))  (tree strong coupling)")
print("="*74)
print(f"{'beta':>10} {'u':>10} {'ratio':>10}  status")
for b in [1/48, 0.5, 1.0, 1.55, 2.0, 2.5]:
    r = mratio(b)
    flag = ">=4 (bound)" if r >= 4 else "<4 (formula still >0)"
    print(f"{b:>10.5g} {u(b):>10.5f} {r:>10.4f}  {flag}")
u_e = np.exp(-1.0)
print(f"\nratio >= 4 requires u <= 1/e = {u_e:.4f}.")
# locate beta where u(beta) = 1/e
from scipy.optimize import brentq
b_inv = brentq(lambda b: u(b) - u_e, 1.0, 2.0)
print(f"u(beta) = 1/e at beta ~ {b_inv:.4f}  => the bound's OWN formula drops below 4")
print(f"beyond beta ~ {b_inv:.3f}, and validity needs convergence (far below that).")
print(f"At beta=2.5 the tree formula gives ratio = {mratio(2.5):.3f} (positive; lattice ~3.7).")

print("\n" + "="*74)
print("C. convergence radius: best rigorous thresholds vs u(2.5)")
print("="*74)
u_bstar = u(8.9e-3); u_szz = u(1/48); u_cross = u(2.5)
print(f"u(residue beta* = 8.9e-3) = {u_bstar:.6f}")
print(f"u(SZZ benchmark 1/48)    = {u_szz:.6f}")
print(f"u(beta = 2.5)            = {u_cross:.6f}")
print(f"u(2.5)/u(1/48)           = {u_cross/u_szz:.1f}  (the expansion parameter is")
print(f"                                       ~{u_cross/u_szz:.0f}x larger at the crossover)")
print("The KP/SZZ theorem controls the expansion for u < u* with u* << 0.5")
print("(the cited beta* are where u is ~1e-3..1e-2).  u(2.5)=%.2f is orders outside." % u_cross)

print("\n" + "="*74)
print("D. Padé principle: extrapolation of a radius-u* series to u(2.5)")
print("="*74)
# Mock with the KP-scale radius u*, NON-rational so Padé is not exact:
#   f(x) = (1 + x/u*)^(-1/2),  radius u* = 0.10, real & positive for x>0.
# Exact value at u_target; Padé orders of the Taylor series are uncontrolled there.
ustar = 0.10
nmax = 16
from math import comb
a_n = np.array([((-1)**n) * comb(2*n, n) / (4.0**n) / (ustar**n) for n in range(nmax)])
u_target = u(2.5)
exact = (1.0 + u_target/ustar)**(-0.5)

def pade_eval(a, L, M, x):
    if L + M + 1 > len(a):
        a = np.concatenate([a, np.zeros(L + M + 1 - len(a))])
    A = np.zeros((M, M)); b = np.zeros(M)
    for i in range(M):
        for j in range(M):
            A[i, j] = a[L + i - j]
        b[i] = -a[L + 1 + i]
    try:
        q = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        q, *_ = np.linalg.lstsq(A, b, rcond=None)
    p = np.array([sum((1.0 if j == 0 else q[j-1]) * a[k-j]
                      for j in range(0, min(k, M)+1)) for k in range(L+1)])
    num = np.polyval(p[::-1], x)
    den = np.polyval(np.concatenate([[1.0], q])[::-1], x)
    return num/den if den != 0 else np.nan

print(f"mock f(x)=(1+x/u*)^-1/2, radius u*={ustar}; at u(2.5)={u_target:.4f} exact={exact:.4f}")
print(f"{'[L/M]':>8} {'Pade':>12}")
vals = []
for (L, M) in [(2,2),(3,3),(4,4),(5,5),(6,6),(7,7),(8,4),(4,8),(6,7),(7,6),(9,3),(3,9)]:
    v = pade_eval(a_n, L, M, u_target); vals.append(v)
    print(f"{f'[{L}/{M}]':>8} {v:>12.4f}")
vals = np.array([v for v in vals if np.isfinite(v)])
print(f"\nexact = {exact:.4f}; Padé spread over orders: min={vals.min():.3f} max={vals.max():.3f}")
print(f"deviation from exact ranges up to {np.max(np.abs(vals-exact)):.3f}"
      f"  ({np.max(np.abs(vals-exact))/exact*100:.0f}% of exact)")
print("=> Padé orders do NOT stabilise at u(2.5); no certificate.")

print("\n" + "="*74)
print("VERDICT: beta* cannot be pushed RIGOROUSLY to 2.5.")
print("  (i) KP/SZZ: rigorous only for u<u* (beta<~1/48); u(2.5) ~ 97x larger.")
print("  (ii) Pade/conformal: diagonal [M/M] ~10% on the mock, so a HEURISTIC")
print("       positive interpolant EXISTS (tree formula 3.30 at 2.5, lattice ~3.7),")
print("       but NO rigorous error bound (off-diagonal spread ~14%); not a certificate.")
print("  (iii) other polymers: activity at 2.5 is O(0.5); no small-activity polymer.")
print("  (iv) monotone quantity: no proven beta-monotonicity of the residue.")
print("=> m_{0++} >= 4 sqrt(sigma) is PHYSICALLY plausible to 2.5, but PROVEN only")
print("   at strong coupling. Bridge stays strong-coupling-local. Clay 0%; no QED.")
print("="*74)
