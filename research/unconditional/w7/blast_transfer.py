"""
blast_transfer.py -- witness for BLAST-TRANSFER.md  (2026-11-24)

Transfer-matrix / variational route to H1: chi = Var(a_{k->0}) = sum_n |c_n|^2/m_n^2 <= 1/lambda.

Parts (numpy only, exact where possible):
 1. The sharp variational UPPER bound on Var(a) = operator-range bound ((a_max-a_min)/2)^2
    (attained by the two-extremal superposition) -- genuine [T] but VACUOUS for the coarse mode.
 2. Free scalar on Z_n: coarse (zero) mode Var = 1/(n m^2); massless -> inf (acid);
    the only finite-L gap is the box gap 4 sin^2(pi/n) -> 0.
 3. Direction obstruction: Dirichlet/subspace bounds on chi = ||(H-E0)^-1 xi||^2 ascend to chi
    from BELOW -- never an upper bound (wrong direction for H1).
 4. The squeeze: chi <= C0/Delta^2 ; chi grows as Delta^2 -> 0, so an UPPER bound on chi needs a
    LOWER bound on the gap Delta^2 (= H1).

Run:  py -3.12 -I research/unconditional/w7/blast_transfer.py
"""
import numpy as np

print("=" * 72)
print("BLAST-TRANSFER witness -- transfer-matrix variational route to H1")
print("=" * 72)

# ---------------------------------------------------------------------------
# Part 1: the sharp variational upper bound = operator-range bound (range/2)^2
# ---------------------------------------------------------------------------
print("\n[1] max_psi Var_psi(a) = ((a_max-a_min)/2)^2  (two-extremal superposition)")
# 2-level operator a = diag(a_min, a_max); for real normalized psi, P_i = psi_i^2
a_min, a_max = -3.0, 7.0
rng = a_max - a_min
a = np.array([a_min, a_max])
def var_state(psi):
    P = psi ** 2                       # probabilities (psi real, normalized)
    m = P @ a                          # <psi|a|psi>
    return P @ (a - m) ** 2            # Var_psi(a)
# maximising state = equal superposition of the two eigenvectors (probs 1/2,1/2)
psi = np.array([1.0, 1.0]) / np.sqrt(2)
print(f"  a_min={a_min}, a_max={a_max}, range={rng}")
print(f"  Var at two-extremal superposition = {var_state(psi):.4f}   (range/2)^2 = {(rng/2)**2:.4f}")
# check it is the max: random states
best = 0.0
r = np.random.default_rng(0)
for _ in range(20000):
    v = r.normal(size=2); v /= np.linalg.norm(v)
    best = max(best, var_state(v))
print(f"  max over 20000 random states = {best:.4f} (== range bound, sharp) "
      f"-> range bound is sharp [T]")
# for the coarse mode a = tr F^2 smeared, range ~ O(field range) => O(Lambda^4 V^2) >> 1/lambda  => VACUOUS.

# ---------------------------------------------------------------------------
# Part 2: free scalar on Z_n -- coarse zero-mode Var = 1/(n m^2); box gap -> 0
# ---------------------------------------------------------------------------
print("\n[2] free scalar on Z_n: zero-mode Var = 1/(n m^2); massless -> inf; box gap 4 sin^2(pi/n)")
print("    n    m^2=0.5 Var       m^2=0 Var        box gap 4sin^2(pi/n)")
for n in [4, 8, 16, 32, 64, 128]:
    Vmass = 1.0 / (n * 0.5)
    boxgap = 4.0 * np.sin(np.pi / n) ** 2
    print(f"  {n:4d}   {Vmass:.3e}        inf              {boxgap:.3e}")
print("    => finite-L 'gap' is the box gap, -> 0 as L->inf; massless chi ~ L^2 -> inf (acid: no gap).")

# ---------------------------------------------------------------------------
# Part 3: Dirichlet/subspace principle gives LOWER bounds on chi (wrong direction)
# ---------------------------------------------------------------------------
print("\n[3] Dirichlet: chi = ||(H-E0)^-1 xi||^2 ; subspace restriction gives LOWER bounds")
evals = np.array([1.0, 2.0, 3.0, 4.0, 8.0])          # m_n^2
xi = np.array([1.0, 0.7, 0.5, 0.3, 0.2])              # |c_n| ~ couplings (schematic)
chi_true = np.sum(xi ** 2 / evals ** 2)
print(f"  true chi = sum xi_n^2 / m_n^2 = {chi_true:.4f}")
ok = True
for k in range(1, 6):
    chi_k = np.sum(xi[:k] ** 2 / evals[:k] ** 2)
    ok = ok and (chi_k <= chi_true + 1e-12)
    print(f"    subspace dim {k}: chi_lower = {chi_k:.4f}")
print(f"  all subspace bounds <= true chi ? {ok}  => variational/Dirichlet yields LOWER bounds,"
      f" H1 needs UPPER => WRONG DIRECTION.")

# ---------------------------------------------------------------------------
# Part 4: the squeeze chi <= C0/Delta^2 ; chi grows as Delta^2 -> 0
# ---------------------------------------------------------------------------
print("\n[4] squeeze chi <= C0/Delta^2 ; chi blows up as the gap Delta^2 -> 0 (needs a gap lower bound)")
f2, w2, m2c = 0.5, 0.5, 4.0          # atom weight, continuum weight, continuum mass^2
C0 = f2 + w2
print("    Delta^2     chi        C0         C0/Delta^2   (chi<=1/lam requires Delta^2>=lam*C0)")
for D2 in [1.0, 0.5, 0.2, 0.05, 0.01]:
    chi = f2 / D2 + w2 / m2c
    print(f"    {D2:6.3f}   {chi:9.3f}   {C0:.3f}      {C0/D2:9.3f}")
print("    => an upper bound on chi (H1) is equivalent to a lower bound on the gap Delta^2 (the gap).")

print("\n" + "=" * 72)
print("Verdict: the only non-vacuous variational upper bound is the squeeze-with-gap; the")
print("unconstrained variational bound (range/2)^2 is vacuous; Dirichlet runs the wrong direction.")
print("BLAST-TRANSFER is CIRCULAR: H1 reduces to the uniform gap = H1.  Clay 0%; no QED.")
print("=" * 72)
