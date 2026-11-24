"""
wi_admissible.py -- witness for research/unconditional/WI-ADMISSIBLE.md

Demonstrates the weight-invariance atom signature on the IR-analytic functional
    A(C) = pi / |Im tau_zero|,   C(tau) = w0 + sum_i f2_i * exp(-E_i tau),
by finding the complex-tau zeros (Newton, numpy only).

Sections
  A. PURE ATOM: A is invariant under weight changes (T1)          -> A = E, constant
  B. TWO-LEVEL / CUT: A MOVES when the second weight changes (T2) -> no invariance
  C. QUASI-ATOM: a weak far mode leaves A ~ E (indistinguishable at finite resolution)
  D. VEV CONFOUND: w0 = 0 (no finite zero) / w0 < 0 (degenerate real zero)

NOT a proof. The invariance is exact for the pure atom (T1) and generic drift is T2.
Run: py -3.12 -I research/unconditional/w7/wi_admissible.py
"""
import numpy as np

def mkf(w0, f2, E):
    f2 = np.array(f2, float); E = np.array(E, float)
    def f(z):  return w0 + np.sum(f2 * np.exp(-E * z))
    def fp(z): return -np.sum(E * f2 * np.exp(-E * z))
    return f, fp

def zeros(w0, f2, E, re=(-10, 10), im=(-20, 20), nr=41, ni=81, tol=1e-11, maxit=80):
    f, fp = mkf(w0, f2, E); out = []
    for r0 in np.linspace(re[0], re[1], nr):
        for i0 in np.linspace(im[0], im[1], ni):
            z = r0 + 1j * i0; ok = False
            for _ in range(maxit):
                fz = f(z)
                if abs(fz) < tol: ok = True; break
                d = fp(z)
                if abs(d) == 0: break
                dz = fz / d; z = z - dz
                if not np.isfinite(z) or abs(z) > 1e6: break
                if abs(dz) < 1e-14: ok = True; break
            if ok and np.isfinite(z) and abs(f(z)) < 1e-7:
                out.append(z)
    uniq = []
    for z in out:
        if not any(abs(z - u) < 1e-6 for u in uniq): uniq.append(z)
    return uniq

def A_of(w0, f2, E):
    zs = [z for z in zeros(w0, f2, E) if abs(z.imag) > 1e-6]
    if not zs: return None, None
    z = min(zs, key=lambda z: abs(z.imag))
    return abs(z.imag), np.pi / abs(z.imag)

def banner(s): print("\n" + "=" * 70 + "\n" + s + "\n" + "=" * 70)

banner("A. PURE ATOM  C = w0 + f2*exp(-E*tau),  E=1, w0=1  ->  A should be E=1 for ALL f2 (T1)")
for f2 in [1.0, 0.3, 5.0]:
    im, A = A_of(1.0, [f2], [1.0])
    print(f"  f2={f2:>4}   |Im tau_zero|={im:.6f}   A=pi/|Im|={A:.6f}   (Re part = ln(f2/w0)/E)")

banner("A'. weights change Re(tau_zero), NOT Im: show Re moves, A fixed (atom)")
for f2 in [0.3, 1.0, 5.0]:
    zs = [z for z in zeros(1.0, [f2], [1.0]) if abs(z.imag) > 1e-6]
    z = min(zs, key=lambda z: abs(z.imag))
    print(f"  f2={f2:>4}   tau_zero = {z.real:+.4f} {z.imag:+.4f}i   A={np.pi/abs(z.imag):.6f}")

banner("B. TWO LEVELS  C = 1 + f2a*exp(-E1 tau) + f2b*exp(-E2 tau),  E=(1,4)  -> A MOVES (T2)")
for f2b in [0.02, 0.2, 1.0, 3.0]:
    im, A = A_of(1.0, [1.0, f2b], [1.0, 4.0])
    print(f"  f2b={f2b:>5}   |Im|={im:.6f}   A={A:.6f}")

banner("C. ANY higher level moves A, even a 1%-weight one -> the test is NOT atom-vs-cut (T2)")
im, A = A_of(1.0, [1.0, 0.01], [1.0, 4.0])
print(f"  |Im|={im:.6f}   A={A:.6f}   (pure atom gave A=E1=1; a 1%-weight far level already pushes A to {A:.3f})")

banner("D. VEV CONFOUND  w0 = 0 -> no finite zero ; w0 < 0 -> real (Im=0) degenerate zero")
im0, A0 = A_of(0.0, [1.0], [1.0])
print(f"  w0=0 : A = {A0}   (e^(-tau) has no finite zero -> undefined)")
zsneg = zeros(-1.0, [1.0], [1.0])
near = [z for z in zsneg if abs(z.imag) < 1e-6]
print(f"  w0=-1: pure-real zeros present? {len(near) > 0}   (nearest zero has Im=0 -> A degenerate)")

print("\nSummary: pure single level -> A invariant (=E1, T1); ANY higher level (even 1%) -> A drifts (T2).")
print("         So the test detects SINGLE-EXPONENTIAL dominance, NOT atom-vs-continuum.")
print("         VEV sign/zero breaks the zero (C1). Downgraded lead; (L,a)-uniform version = H1. No QED.")
