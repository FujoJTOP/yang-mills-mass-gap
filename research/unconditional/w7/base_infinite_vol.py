#!/usr/bin/env python3
"""
base_infinite_vol.py -- witness for research/unconditional/BASE-INFINITE-VOL.md
Tiny, numpy/scipy only. NO lattice MC. One process, <1 s.

Question: can the P3 strong-coupling RESIDUE floor f1^2_L >= r0(beta) > 0 (L-uniform,
I1) be extended in beta toward the crossover beta~2.5 WITHOUT crossing it, using a
beta-UNIFORM structural tool (reflection positivity, Gaussian domination, monotone
coupling, correlation inequality)?

Answer the witness makes arithmetic-explicit:
  (1) The ONLY rigorous lower bound on f1^2 is the polymer/KP expansion. Its convergence
      radius beta* is set by  lambda_a * Phi(1) <= 1/2  (Kotecky-Preiss, a=1).
  (2) beta* is O(1e-3 .. 2e-2), i.e. >=2 orders of magnitude BELOW the crossover.
  (3) The four beta-uniform tools give at best f1^2 >= 0 (RP), an UPPER bound (Gaussian
      domination), or nothing (monotone / correlation for a KL atom). None lowers beta*.

Run:  py -3.12 -I research/unconditional/w7/base_infinite_vol.py
"""
import numpy as np
from scipy.optimize import brentq

N, d = 3, 4
Delta = 20.0                  # plaquette adjacency degree (P3-FULLREP-KP §4.1) = 4*(6-1)
lam_a = np.e * Delta          # crude animal constant ~ 54.4  (§2.2, a=1)
CROSSOVER = 2.5               # beta where alpha_s ~ O(1) (the wall)
SZZ_BENCH = 1.0/48.0          # external benchmark radius (Shen-Zhu-Zhu)

def Phi(beta, a=1.0):
    """per-plaquette representation activity  (§2.2, (star-star))."""
    return np.exp(beta) * (np.exp(beta * np.exp(a)) - 1.0)

def kp_le_half(beta, a=1.0):
    """KP at a=1 is satisfied iff (S_x <= a=1) iff lambda_a*Phi <= 1/2."""
    return lam_a * Phi(beta, a) <= 0.5

def beta_star(a=1.0):
    """largest beta with lambda_a*Phi(beta,a) = 1/2."""
    return brentq(lambda b: lam_a * Phi(b, a) - 0.5, 1e-8, 1.0)

print("="*74)
print("BASE-INFINITE-VOL witness  (residue beta*-radius, crude KP)")
print("="*74)
print(f"  N={N}, d={d}, Delta={Delta:.0f}, lambda_a=e*Delta={lam_a:.3f} (crude, unoptimised)")
print(f"  crossover beta = {CROSSOVER}   benchmark radius 1/48 = {SZZ_BENCH:.5f}")

bs = beta_star()
print("\n(1) KP convergence radius (crude constants):")
print(f"    beta* (crude KP, a=1) = {bs:.3e}")
print(f"    orders below crossover : log10(2.5/beta*) = {np.log10(CROSSOVER/bs):.2f}")
print(f"    orders below 1/48      : log10((1/48)/beta*) = {np.log10(SZZ_BENCH/bs):.2f}")

print("\n(2) activity and KP status across beta:")
print(f"    {'beta':>10} {'Phi(1)':>12} {'lam_a*Phi':>12} {'KP ok?':>8}")
for b in [bs, 1e-3, 1e-2, SZZ_BENCH, 0.1, 1.0, CROSSOVER]:
    ph = Phi(b)
    print(f"    {b:10.4g} {ph:12.4g} {lam_a*ph:12.4g} {str(kp_le_half(b)):>8}")

print("\n(3) even OPTIMISING the grading a>0 cannot close the gap:")
best = min(((lam_a*Phi(1e-3, a), a) for a in np.linspace(0.02, 3.0, 400)))
print(f"    min over a of lam_a*Phi(beta=1e-3,a) = {best[0]:.4g} at a={best[1]:.2f}")
# optimised radius: scale beta until lam*Phi = 1/2 for the best a
bs_opt = brentq(lambda b: lam_a*Phi(b, best[1]) - 0.5, 1e-8, 1.0)
print(f"    optimised-a radius beta* ~ {bs_opt:.3e}  (still {np.log10(CROSSOVER/bs_opt):.1f} orders short)")

print("\n(4) what each beta-UNIFORM tool yields for the residue floor f1^2_L >= ? :")
print("    RP across L          : f1^2 >= 0        (trivial; no beta dependence)")
print("    Gaussian domination  : UPPER bound only (|<...>| <= ...); no lower floor")
print("    monotone coupling    : no proven monotonicity of a KL atom weight in beta")
print("    correlation inequality: sign/expectation bounds, not spectral weights")
print("    => the ONLY lower bound is the polymer/KP expansion, valid iff beta<beta*.")

print("\n(5) acid test (free-field / beta->0):")
u = lambda b: b/(2*N)          # fundamental coefficient c_f ~ beta/(2N) (literal coords)
for b in [1e-4, 1e-2, bs]:
    print(f"    beta={b:.2e}: tube weight ~ u^4 = {u(b)**4:.2e} -> 0  (f1^2 -> 0; free field)")

print("\n"+"="*74)
print("VERDICT: beta* is the strong-coupling/KP radius, UNMOVED by any of the four")
print("beta-uniform tools. Exact beta* ~ %.1e (crude) .. %.5f (benchmark)." % (bs, SZZ_BENCH))
print("It sits >=2 orders of magnitude below the crossover beta=%.1f. The residue" % CROSSOVER)
print("cannot be carried to the crossover by a beta-uniform inequality; the only")
print("exit is a rigorous 4D-gauge RG, which does not exist. No QED.")
print("="*74)
