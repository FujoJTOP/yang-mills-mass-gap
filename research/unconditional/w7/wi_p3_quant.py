"""
wi_p3_quant.py -- WI LINE 1 (P3 approach law): is the isolation Delta=e2-e1
QUANTITATIVELY extractable, and is Delta>=delta certifiable from finite data?

Effective mass   m(t) = -d/dt log C(t),
approach rate    kappa(t) = -d/dt log( m(t) - e1 ).

EXACT closed forms (e1=0 WLOG, atom weight c1=1, weight ratio r=c2/c1):
  atom + gap (2 level):   kappa(t) = Delta / (1 + r e^{-Delta t})      -> Delta
  atom + flat threshold:  kappa(t) = (2 c1 t + rho0) / (t (c1 t + rho0)) -> 0
GAP  => kappa increases to Delta.   CUT => kappa decreases to 0.
Key exact fact (2 level): kappa(t) <= Delta for ALL t  ==>  Delta >= kappa(t)
(a one-sided certified LOWER bound on the isolation, tightening to Delta).
Run: py -3.12 -I wi_p3_quant.py  (numpy only)
"""
import numpy as np

r, c1, rho0 = 0.3, 1.0, 0.3

def kappa_gap(t, D):  return D / (1.0 + r*np.exp(-D*t))     # exact 2-level
def kappa_cut(t):     return (2*c1*t + rho0) / (t*(c1*t + rho0))  # exact flat cut
def tight(t, D):      return 1.0 / (1.0 + r*np.exp(-D*t))  # kappa/Delta

taus = [5, 10, 20, 30, 50, 100]
print("="*72)
print("kappa(t) = -d/dt log(m-e1):  -> Delta for a GAP, -> 0 for a THRESHOLD")
print("="*72)
for D in [0.5, 0.2, 0.05]:
    print(f" GAP  Delta={D:<5}: " + str([round(float(kappa_gap(t, D)), 4) for t in taus]))
print(f" CUT (threshold)  : " + str([round(float(kappa_cut(t)), 4) for t in taus]))

print()
print("="*72)
print("Certifiability:  Delta >= kappa(t),  kappa(t)/Delta = 1/(1+r e^{-Delta t})")
print("  (how close the one-sided lower bound gets to the true Delta)")
print("="*72)
print("  Delta | " + " | ".join(f"t={t:<3}" for t in [5, 20, 50, 100]))
for D in [0.5, 0.2, 0.1, 0.05, 0.02]:
    row = "  ".join(f"{tight(t, D):.2f}" for t in [5, 20, 50, 100])
    print(f"  {D:<5} | {row}")
print()
print("  reading: a tight (>=90%) lower bound needs Delta*tau >~ 1;")
print("  thin isolation (Delta*tau_max << 1) only yields ~1/(1+r) of Delta")
print("  => conservative lower bound, false-negatives for thin Delta, never false positives.")

# multi-level caveat: far heavy level perturbs the SMALL-t rate but decays before the Delta-window
print()
print("="*72)
print("Multi-level caveat: far heavy level (e3=100,w=1e-3) perturbs only tiny t")
print("="*72)
e3 = 100.0
def m_3lev(t, D, c2=0.01, c3=1e-3):
    C = 1.0 + c2*np.exp(-D*t) + c3*np.exp(-e3*t)
    num = c2*D*np.exp(-D*t) + c3*e3*np.exp(-e3*t)
    return num/C
def kappa_num(t, D, h=1e-7):
    u = lambda s: m_3lev(s, D) - 0.0
    return -(np.log(u(t+h)) - np.log(u(t-h)))/(2*h)
print(f"  3-level (Delta=1, far e3=100): kappa at t=[0.002,0.005,0.01,0.05,0.2,1.0]:")
print("    " + str([round(float(kappa_num(t, 1.0)), 3) for t in [0.002, 0.005, 0.01, 0.05, 0.2, 1.0]]))
print("  -> far level inflates kappa only while e3*t<~1 (t<~0.01); over the Delta-window")
print("     (t~1) it is fully decayed and kappa -> Delta=1. So the rate is well-defined.")
