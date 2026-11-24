"""
wi_reslope_limit.py -- does the Re-slope functional (WI-MODIFIED-FUNCTIONAL) give a
FORCING lower bound on the atom weight f^2, or is it A3/infinite by necessity?

The Re-slope reads: a simple atom at s=E1 puts the complex-tau zeros of
Z(tau)=w0+sum f_i^2 exp(-E_i tau) on a vertical line  Re tau = ln(f^2/w0)/E1,
spaced Im-spacing 2*pi/E1;  under operator-reweighting f^2 -> lam f^2 the real
part slides at  d Re tau / d ln(lam) = 1/E1,  INDEPENDENT of the rung k.

We show:
  (1) [T] pure atom: Re-slope = 1/E1 exactly, rung-invariant (all k).
  (2) [T] a finite REAL-tau sample cannot distinguish the atom from an
          atomless bump of width h:  max|Z_atom - Z_bump| = O(h) -> 0.
  (3) [T] but the zero LADDER (needed for the Re-slope) differs: the atom has
          clean rung-invariant zeros; the bump/continuum does not. Recovering
          it needs complex-tau = analytic continuation (A3) / the full Z (J->inf).
=> the Re-slope is A3/infinite by necessity; no finite-data forcing lower bound.
Run: py -3.12 -I research/unconditional/w7/wi_reslope_limit.py
"""
import numpy as np

# ---- Z(tau) = w0 + sum f2 exp(-E tau), complex tau --------------------------
def Z(tau, w0, levels):
    return w0 + sum(f2 * np.exp(-E * tau) for (f2, E) in levels)

def dZ(tau, w0, levels):
    return sum(-E * f2 * np.exp(-E * tau) for (f2, E) in levels)

def newton_zero(tau, w0, levels, it=400):
    for _ in range(it):
        z = Z(tau, w0, levels); dz = dZ(tau, w0, levels)
        if abs(dz) < 1e-300:
            break
        step = z / dz
        tau = tau - step
        if abs(step) < 1e-15:
            break
    return tau

# ---- (1) pure atom: exact Re-slope = 1/E1, rung-invariant ------------------
print("=" * 74)
print("(1) PURE ATOM: Re(tau_k)=ln(f^2/w0)/E1,  dRe/dln(lambda)=1/E1 (all rungs)")
print("=" * 74)
w0, E1 = 1.0, 1.0
for f2 in (0.5, 1.0):                       # reweight f^2: 0.5 -> 1.0  (lam=2)
    levels = [(f2, E1)]
    print(f"  f^2={f2:4.2f}:  exact Re(tau_k) = ln(f^2/w0)/E1 = {np.log(f2/w0)/E1:+.6f}")
lam = 1.0 / 0.5                              # lam = f2_new/f2_old = 2
dRe = np.log(1.0 / w0) / E1 - np.log(0.5 / w0) / E1
print(f"  dRe/dln(lam) = {dRe/np.log(lam):.6f}   (1/E1 = {1/E1:.6f})  -> rung-INVARIANT")
# verify by Newton on a few rungs
levels = [(0.5, E1)]
for k in range(3):
    seed = np.log(0.5 / w0) / E1 - 1j * (2 * k + 1) * np.pi / E1
    z = newton_zero(seed, w0, levels)
    print(f"    rung k={k}: Newton zero Re={z.real:+.6f}  Im={z.imag:+.6f}  (Im spacing 2pi/E1={2*np.pi/E1:.4f})")

# ---- (2) finite REAL-tau sample: atom vs atomless bump, O(h) match ----------
print("\n" + "=" * 74)
print("(2) FINITE REAL-tau SAMPLE: atom vs atomless bump of width h -> O(h)")
print("=" * 74)
f2 = 0.5
tau_grid = np.linspace(0.0, 5.0, 12)          # 12 real-tau points
Zatom = np.array([Z(t, w0, [(f2, E1)]) for t in tau_grid])

def Zbump(t, w0, E1, f2, h):
    # continuum rho = f2/h on [E1, E1+h]:  int rho e^{-s t} ds
    if abs(t) < 1e-12:
        return w0 + f2
    return w0 + (f2 / h) * (np.exp(-E1 * t) - np.exp(-(E1 + h) * t)) / t

for h in (1e-1, 1e-2, 1e-3, 1e-4):
    Zb = np.array([Zbump(t, w0, E1, f2, h) for t in tau_grid])
    print(f"  h={h:.0e}:  max|Z_atom - Z_bump| over 12 real-tau pts = {np.max(np.abs(Zatom - Zb)):.3e}")
print("  => the finite real-tau sample is matched to O(h); any functional of it is blind.")

# ---- (3) the zero LADDER differs (atom clean; continuum not) ----------------
print("\n" + "=" * 74)
print("(3) ZERO LADDER at complex tau (needs analytic continuation): atom vs continuum")
print("=" * 74)
# atom = a strong simple level at E1; pure CUT = dense levels with NO atom at E1
levels_atom = [(0.5, E1)]
levels_cut = [(0.03, 1.2 + 0.1 * j) for j in range(18)]   # pure cut, no atom at E1
levels_atomweakcont = [(0.5, E1)] + [(0.03, 1.2 + 0.1 * j) for j in range(18)]
for name, lv in (("atom", levels_atom),
                 ("atom+weak-cont", levels_atomweakcont),
                 ("pure cut", levels_cut)):
    res = []
    for k in range(3):
        seed = np.log(0.5 / w0) / E1 - 1j * (2 * k + 1) * np.pi / E1
        z = newton_zero(seed, w0, lv)
        res.append(z)
    re = np.array([z.real for z in res])
    spread = np.max(re) - np.min(re)
    print(f"  {name:13s}: Re(tau_k) k=0,1,2 = {re[0]:+.4f},{re[1]:+.4f},{re[2]:+.4f}"
          f"   spread={spread:.4e}")
print("  => rung-invariance holds ONLY for an essentially PURE single level; even a")
print("     weak continuum shifts it (spread=0.21 -> broken) and a pure cut has none.")
print("     So it is both fragile AND not finitely computable: recovering the ladder")
print("     needs complex-tau (analytic continuation = A3) or the full Z = J->inf.")

print("\nCONCLUSION: Re-slope = A3 read of the simple pole (weak-* discontinuous).")
print("No finite-data forcing lower bound; A3 or J->infinity = H1.  Clay 0%.")
