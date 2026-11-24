"""
cert_chi_sharp.py -- CERT-CHI-SHARP witness.
Computes the marginal inequality `mu(1 - gB*alpha) >= lambda` that is the
chi-half of H1 (BLAST-OPERATOR-INEQ / BLAST-YM-NONLINEAR):
  Var(a_{k->0}) <= 1/(mu(1 - c*alpha)),  c = gB,  mu = lambda_min(H+),
  alpha = phi_hat^T H+^{-1} phi_hat,   valid iff c*alpha < 1.
  H1 (Var <= 1/lambda) <=> mu(1 - c*alpha) >= lambda.

The Savvidy one-loop Landau levels:  lambda_{n,s} = (2n+1) gB + 2 gB s,  s = -1,+1.
The negative mode is (n=0, s=-1) at -gB; c = gB.
Positive modes: (0,+1)=3gB, (1,-1)=gB, (1,+1)=5gB, (2,-1)=3gB, ...

Overlap model: the negative mode is the n=0 Landau wavefunction psi_0.  It overlaps
a positive mode ONLY if that mode has the same n=0 (same psi_0); n>=1 modes are
orthogonal.  So alpha = 1/(eigenvalue of the (0,+1) mode) = 1/(3gB) in this model.
  c*alpha = gB/(3gB) = 1/3,  mu = min positive = gB (the (1,-1) mode)
  => mu(1-c*alpha) = gB*(1-1/3) = (2/3) gB.

Alternative (fully degenerate H+ = mu*I with the neg mode overlapping the mu-mode):
  alpha = 1/mu = 1/gB, c*alpha = 1, mu(1-c*alpha) = 0.   <-- EXACTLY marginal.

So the certificate LHS ranges over [0, (2/3) gB] depending on the overlap/spectral
structure -- always O(gB).  Compare to lambda = 1/chi (lattice).

Run: py -3.12 -I research/unconditional/w7/cert_chi_sharp.py
"""
import numpy as np

gB = 1.0  # work in units of gB
levels = []
for n in range(0, 6):
    for s in (-1, +1):
        levels.append((n, s, (2 * n + 1) * gB + 2 * gB * s))
levels.sort(key=lambda t: t[2])

print("Lowest Savvidy Landau levels (n, s, lambda/gB):")
for n, s, lam in levels[:8]:
    print(f"   n={n}  s={s:+d}   lambda/gB = {lam/gB:+.2f}")

neg = [t for t in levels if t[2] < 0]
pos = [t for t in levels if t[2] > 0]
c = -min(t[2] for t in neg)      # negative-mode magnitude
mu = min(t[2] for t in pos)      # smallest positive eigenvalue of H+
print(f"\nc   (= gB, negative-mode magnitude)      = {c/gB:.4f} gB")
print(f"mu  (= min positive eigenvalue of H+)   = {mu/gB:.4f} gB   [the (1,-1) mode]")

# --- overlap model alpha = sum over positive modes of |<neg|pos>|^2 / lambda_pos,
#     with |<psi_0 | psi_{n}>|^2 = 1 if n=0 else 0.
alpha = sum(1.0 / lam for (n, s, lam) in pos if n == 0)
c_alpha = c * alpha
lhs_overlap = mu * (1 - c_alpha)
print(f"\nalpha (n=0-overlap model)                = {alpha/gB:.4f} / gB")
print(f"c*alpha                                  = {c_alpha:.4f}")
print(f"mu(1 - c*alpha)                          = {lhs_overlap/gB:.4f} gB   <-- LHS (overlap model)")

# --- fully degenerate H+ = mu*I, neg mode overlaps the mu-eigenvector:
alpha_full = 1.0 / mu
c_alpha_full = c * alpha_full
lhs_full = mu * (1 - c_alpha_full)
print(f"\n[fully degenerate H+=mu*I]  c*alpha      = {c_alpha_full:.4f}")
print(f"[fully degenerate]          mu(1-c*alpha) = {lhs_full/gB:.4f} gB  <-- LHS (degenerate model)")

# --- the lambda side: 1/chi = G'' measured on the lattice (corpus, finite volume)
print("\nlambda = 1/chi = G'' on the lattice (corpus finite-volume values):")
for L, Gpp in [(4, 0.28), (8, 0.04)]:
    print(f"   L={L}:  G'' = 1/chi = {Gpp}")
print("   -> in lattice units both LHS = O(gB) and lambda = 1/chi are O(1): SAME order.")

# --- crossing: LHS = mu(1-c*alpha) as a function of the overlap strength
print("\nCertificate LHS vs overlap strength t := <neg|mu-mode> (0..1):")
for t in [0.0, 0.25, 0.5, 0.7071, 0.9, 1.0]:
    # alpha = t^2/mu_eff ; take mu_eff = mu (the lightest positive)
    a_t = (t ** 2) / mu
    print(f"   t={t:+.4f}:  c*alpha={c*a_t:+.4f}   mu(1-c*alpha)/gB = {mu*(1-c*a_t)/gB:+.4f}")

print("\nSUMMARY (honest):")
print(" - The LHS mu(1-gB*alpha) is O(gB), ranging over [0, (2/3)gB] for the Landau")
print("   structure (0 in the fully-degenerate single-mode limit).")
print(" - lambda = 1/chi is O(1) (measured G'' ~ 0.04-0.28 in lattice units).")
print(" - mu, alpha, gB and lambda are ALL nonperturbative O(1) at beta=2.5: no")
print("   hierarchy, no free pass.  The two sides are the SAME order -> marginal.")
print(" - The certificate mu(1-gB*alpha) >= lambda is NOT closeable exactly: it is")
print("   the O(1) competition = H1.  The negative mode eats O(1) of the floor.")
