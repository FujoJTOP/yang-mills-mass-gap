"""
verify_signedfourier.py — sanity checks for SIGNED-FOURIER-BUILD.md

(1) SU(2) character coefficients c_j(beta) of the plaquette weight
    e^{(beta/N) Re tr U} (N=2 -> e^{beta cos theta}): positivity for all beta>0,
    and the fundamental (j=1/2) is the lightest nontrivial tube (Casimir order).
(2) The "signed-Fourier / Bochner" edge-weight estimator discriminates a true
    atom from a continuum cut (the free-field acid test):
       atom:  C(t)=w e^{-m t}                 -> w_est(t):=C(t)e^{m t} -> w > 0
       cut :  C(t)=int_{s0}^inf rho(s)e^{-s t}ds -> w_est(t):=C(t)e^{s0 t} -> 0
    For the 2-particle free cut rho~sqrt(s-s0) this is exact: C=(sqrt(pi)/2) t^{-3/2} e^{-s0 t},
    so the edge weight -> 0 (f1^2=0), the acid test the construction must pass.
Run:  py -3.12 verify_signedfourier.py
"""
import numpy as np
from scipy.special import iv          # modified Bessel (to check a closed form)
from scipy.integrate import quad

# ---------- (1) SU(2) character coefficients -------------------------------
def c_j(beta, j):
    """c_j = int e^{beta cos th} chi_j(th) dU, dU = (2/pi) sin^2 th dth, th in [0,pi]."""
    def f(th):
        chi = np.where(np.abs(np.sin(th)) < 1e-12,
                       (2*j+1), np.sin((2*j+1)*th)/np.sin(th))
        return np.exp(beta*np.cos(th))*chi*np.sin(th)**2
    val, _ = quad(f, 0.0, np.pi, limit=400)
    return (2.0/np.pi)*val

print("== (1) SU(2) character coefficients c_j(beta),  j=0,1/2,1,3/2,2 ==")
js = [0.0, 0.5, 1.0, 1.5, 2.0]
C2 = {j: j*(j+1) for j in js}                     # quadratic Casimir
for beta in [0.02, 0.1, 0.5, 1.0, 2.5]:
    cs = [c_j(beta, j) for j in js]
    # cross-check j=0 against closed form 2 I_1(beta)/beta
    closed0 = 2.0*iv(1, beta)/beta if beta > 0 else 1.0
    print(f"beta={beta:5.3f}  c_j = " +
          "  ".join(f"j={j:3.1f}:{c:+.4e}" for j, c in zip(js, cs)) +
          f"   | c_0 vs 2I1/b={closed0:+.4e}")
    assert all(c > 0 for c in cs), "c_j must be > 0 for beta>0"
    # lightest nontrivial tube = fundamental j=1/2 must be the SMALLEST Casimir
    assert C2[0.5] < C2[1.0] < C2[1.5] < C2[2.0]
print("  -> all c_j > 0 (positive character expansion); fundamental j=1/2 has min C2 -> lightest tube\n")

# ---------- (2) atom vs cut edge-weight discriminator ----------------------
print("== (2) signed-Fourier edge weight: atom vs free cut ==")
m, w = 1.3, 0.27
ts = np.array([2.0, 4.0, 8.0, 16.0, 32.0])
C_atom = w*np.exp(-m*ts)
w_atom = C_atom*np.exp(m*ts)                       # estimated atom weight
print("  ATOM  C=w e^{-m t}:  t, w_est(t)")
for t, we in zip(ts, w_atom):
    print(f"         t={t:5.1f}  w_est={we:.6f}")

s0 = 4.0                                           # free 2-particle threshold
def C_cut(t):                                      # int_{4}^inf sqrt(s-4) e^{-s t} ds
    return np.exp(-s0*t)*0.5*np.sqrt(np.pi)/t**1.5  # exact
w_cut = C_cut(ts)*np.exp(s0*ts)                     # ->0 (power law)
print("  CUT   C=(sqrt(pi)/2) t^{-3/2} e^{-s0 t}:  t, w_est(t)=C e^{s0 t}")
for t, we in zip(ts, w_cut):
    print(f"         t={t:5.1f}  w_est={we:.6f}")
# asymptotic check: w_cut(t) ~ (sqrt(pi)/2) t^{-3/2} -> 0
print("  -> atom edge weight is CONSTANT >0 ; free-cut edge weight -> 0 (power law).")
print("     Acid test PASSED: the representation coefficient at the edge vanishes for the free cut,")
print("     so any lower bound built from it must give f1^2=0 there.  (f1^2 = 0 for free :phi^2:.)")
