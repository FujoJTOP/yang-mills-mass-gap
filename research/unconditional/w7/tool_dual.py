"""
tool_dual.py — TOOL-DUAL witness: can convex / Legendre-Fenchel / Nehari duality
FORCE a lower bound on the KL atom  f1^2 = nu({m1^2})  from the integral data
(susceptibility chi = ∫ s^{-1} dnu, total weight M = ∫ dnu, support floor m1^2)?

Setup: nu >= 0 supported on [m1^2, inf).  Atom a := nu({m1^2}).
We solve the measure LP  min a  s.t.  ∫dnu = M, ∫s^{-1}dnu = chi, nu>=0,
discretised on a grid containing m1^2, and compare with the Legendre-Fenchel dual.

Honest expectation (this IS the test): the atomless measure is feasible, so the
primal infimum is 0 and — by EXACT LP duality — the dual optimum is also 0.
Duality therefore CANNOT force f1^2 > 0. The witness confirms it numerically.
"""
import numpy as np
from scipy.optimize import linprog

np.set_printoptions(precision=5, suppress=True)


def measure_lp(M, chi, m1sq, smax, n):
    """min w[0]  s.t. sum w = M, sum w/s = chi, w>=0, grid s in [m1sq,smax]."""
    s = np.linspace(m1sq, smax, n)
    s[0] = m1sq                      # atom sits exactly at the floor
    A_eq = np.vstack([np.ones_like(s), 1.0 / s])   # moments p=0 (M) and p=-1 (chi)
    b_eq = np.array([M, chi])
    c = np.zeros(n); c[0] = 1.0      # objective: the atom weight at s=m1sq
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * n, method="highs")
    # dual: max b·y s.t. A_eq^T y <= c  (y = equality marginals)
    y = res.eqlin.marginals if res.eqlin is not None else np.zeros(2)
    dual_obj = float(b_eq @ y)
    w = res.x
    return dict(primal=float(res.fun), dual=dual_obj, y=np.asarray(y),
                atom=float(w[0]), w=w, s=s)


def atomless_solution(M, chi, m1sq):
    """Analytic atomless feasible point: put all mass at s0 = M/chi (>= m1^2 iff chi<=M/m1^2)."""
    s0 = M / chi
    return s0, (s0 >= m1sq - 1e-12)


def extra_moment_lp(M, chi, A, m1sq, smax, n):
    """Add the trace-anomaly moment  ∫ s^{-2} dnu = A.  Still linear -> still blind?"""
    s = np.linspace(m1sq, smax, n); s[0] = m1sq
    A_eq = np.vstack([np.ones_like(s), 1.0 / s, 1.0 / s**2])
    b_eq = np.array([M, chi, A])
    c = np.zeros(n); c[0] = 1.0
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * n, method="highs")
    if not res.success:
        return None
    y = res.eqlin.marginals
    return dict(primal=float(res.fun), dual=float(b_eq @ y), atom=float(res.x[0]))


def main():
    print("=" * 74)
    print("TOOL-DUAL : does convex duality force the KL atom f1^2 > 0 ?")
    print("=" * 74)

    # ---- Data: honestly chosen to be FEASIBLE and gapped (chi <= M/m1^2) ----
    M, chi, m1sq = 1.0, 0.5, 1.0
    print(f"\nData: M={M}, chi={chi}, m1^2={m1sq}  (feasibility chi*m1^2={chi*m1sq} <= M={M}: {chi*m1sq<=M})")

    # ---- (1) The measure LP and its exact dual ----
    for smax, n in [(20.0, 401), (50.0, 801), (200.0, 2001)]:
        r = measure_lp(M, chi, m1sq, smax, n)
        print(f"  grid smax={smax:6.1f} n={n:5d}:  primal(min atom)={r['primal']:.3e}   "
              f"dual={r['dual']:.3e}   y=(M,chi)-multipliers={r['y']}")

    # ---- (2) Explicit atomless feasible point (mass at s0 = M/chi) ----
    s0, ok = atomless_solution(M, chi, m1sq)
    print(f"\nAtomless feasible point: put all mass at s0 = M/chi = {s0:.4f}  (s0>=m1^2: {ok})")
    print(f"  => a = nu({{m1^2}}) = 0 is FEASIBLE.  primal infimum = 0.")

    # ---- (3) Add the trace-anomaly moment: still blind? ----
    # Choose A consistent with an atomless continuous spectrum over [m1^2, smax]: pick a test.
    A = 0.35
    print(f"\nAdd trace-anomaly moment  ∫s^(-2)dnu = A = {A} (3 constraints):")
    for smax, n in [(50.0, 1001), (200.0, 4001)]:
        r = extra_moment_lp(M, chi, A, m1sq, smax, n)
        if r is None:
            print(f"  smax={smax}: infeasible at this A (data infeasible)")
        else:
            print(f"  smax={smax:6.1f} n={n:5d}:   primal(min atom)={r['primal']:.3e}   dual={r['dual']:.3e}")

    # ---- (4) Shadow price of the atom constraint: is there a finite dual witness? ----
    r = measure_lp(M, chi, m1sq, 50.0, 1201)
    print(f"\nDual multipliers at optimum: y = {r['y']}  (both ~0 => no finite dual witness pushes the atom up)")

    # ---- (5) The 'forcing' would need the rate-function J -> infinity ----
    print("\nLegendre-Fenchel reading: forcing the atom needs a source parameter J -> infinity")
    print("  (the max-entropy / exponential-family solution is CONTINUOUS for finite J;")
    print("   an atom is only recovered in the J->inf limit = the full Stieltjes transform).")

    print("\nVERDICT: primal inf = dual optimum = 0.  Duality is EXACT and BLIND.")
    print("  Strong duality HOLDS (gap = 0) --- that is precisely why no dual trick works:")
    print("  the atomless measure is feasible, so both primal and dual optima are 0.")


if __name__ == "__main__":
    main()
