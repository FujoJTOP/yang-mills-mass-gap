"""
wi_continuum_subtract.py -- witness for WI-CONTINUUM-SUBTRACT (numpy only).

Spectral model for the 0++ channel:
    Z_O(tau) = w0 + a_O^2 e^{-E1 tau} + sum_k g_O[k]^2 e^{-E_k tau}
    (atom at E1 weight a_O^2; continuum modes E_k with amplitudes g_O[k])

Checks:
 (A) Densities are QUADRATIC in the operator. For a linear reweighting
     O_c = O - c O', the continuum density is |g_O - c g_O'|^2 (cross terms).
     Minimising it over c leaves a residual > 0 unless g_O || g_O' (Cauchy-Schwarz),
     and if the vectors are not parallel the atom does NOT cancel at the same c.
 (B) A DESIGNED pair with identical continuum vector g but different atom
     amplitude cancels the continuum exactly (c=1) and leaves a pure atom.
 (C) The nearest complex-tau zero recovers E1 exactly for a pure atom, and its
     contamination shrinks as the continuum threshold theta/E1 grows (a far
     continuum = isolation). Quantifies why 1% continuum drifts A.
"""
import numpy as np

E1, w0 = 1.0, 1.0
K, theta, W = 14, 2.5, 2.0
Ek = np.linspace(theta, theta + W, K)
dk = W / (K - 1)
rng = np.random.default_rng(0)

# ---- (A) generic pair ------------------------------------------------------
aO, aOp = 1.0, 0.5
gO = rng.normal(0, 1, K) + 2.0
gOp = rng.normal(0, 1, K) + 2.0        # non-parallel to gO

def cont_resid(g, gp, c):
    return np.sum((g - c * gp) ** 2) * dk

cs = np.linspace(-6, 6, 24001)
res = np.array([cont_resid(gO, gOp, c) for c in cs])
i0 = int(np.argmin(res)); c0 = cs[i0]
total = np.sum(gO ** 2) * dk
atom_at_c0 = abs(aO - c0 * aOp)
print("(A) generic pair (non-parallel continuum vectors)")
print(f"    total continuum weight  = {total:.4f}")
print(f"    min residual over c     = {res[i0]:.4e}  at c*={c0:.3f}")
print(f"    residual / total        = {res[i0]/total:.4e}  (0 iff gO || gOp)")
print(f"    atom amplitude at c*    = {atom_at_c0:.4f}  (cannot be both: cancel cut AND keep atom)")

# ---- (B) designed pair: identical continuum vector, different atom ---------
gO2, a2 = gO.copy(), 3.0
print("\n(B) designed pair (identical continuum vector, different atom)")
print(f"    continuum residual at c=1: {cont_resid(gO, gO2, 1.0):.4e}")
print(f"    atom amplitude  a - c a' :  {aO - 1.0*a2:.4f}  -> pure atom recovered")

# ---- complex-Newton zero of Z(tau) ----------------------------------------
def Zfun(t, amp, E):
    return w0 + np.sum(amp * np.exp(-E * t))

def dZfun(t, amp, E):
    return np.sum(-E * amp * np.exp(-E * t))

def nearest_zero(t0, amp, E, iters=80):
    t = t0
    for _ in range(iters):
        dz = dZfun(t, amp, E)
        if abs(dz) < 1e-300:
            break
        t = t - Zfun(t, amp, E) / dz
    return t

def A_of(a, g, Ek, scale_E=True):
    amp = np.concatenate(([a ** 2], g ** 2))
    E = np.concatenate(([E1], Ek))
    t0 = (np.log(a ** 2 / w0) - 1j * np.pi) / E1   # pure-atom zero as seed
    t = nearest_zero(t0, amp, E)
    return np.pi / abs(t.imag), t

print("\n(C) nearest-zero value A = pi/|Im tau| (atom weight a=1, continuum weight eps)")
for eps in [0.0, 1e-3, 1e-2]:
    g = np.sqrt((eps / K)) * np.ones(K)          # uniform continuum, total weight eps
    A0, _ = A_of(1.0, g, Ek)
    print(f"    eps={eps:8.3f}:  A = {A0:.5f}   (pure atom = 1.00000)")

print("\n    contamination vs threshold theta/E1 (eps=1e-2):")
for th in [1.5, 2.5, 4.0, 8.0]:
    Ekk = np.linspace(th, th + 1.0, K)
    g = np.sqrt(1e-2 / K) * np.ones(K)
    A0, _ = A_of(1.0, g, Ekk)
    print(f"    theta/E1={th:4.1f}:  A = {A0:.5f}   |A-E1| = {abs(A0-1.0):.2e}")

print("\nSUMMARY:")
print("  * continuum cancels ONLY for a parallel pair (gO||gOp), which also cancels the atom;")
print("  * an EXACT subtraction needs a continuum-null operator pair -- the hard nonperturbative object;")
print("  * a generic local 0++ 2-point function HAS the Kallen-Lehmann cut (multi-glueball), so no")
print("    local continuum-null operator exists => no constructible local subtraction.  A->E1 only in")
print("    the far-continuum (isolation) limit.")
