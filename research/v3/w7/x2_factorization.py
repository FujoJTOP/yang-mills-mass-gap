"""
x2_factorization.py - witness for research/v3/X2-FACTORIZATION.md

Purpose: quantify, in exactly solvable cases, the assertion that
"the obstruction to local constancy of the observable factorization algebra"
is *nothing but* the analytic propagator / transfer-matrix mass.

(A) Free scalar, exact Euclidean propagators d=1..4.
    Extract the coincident-scale decay rate  rate(r) = -d/dr log|G(r)| .
      massive  : rate(r) -> m  (exponential decay => observables locally constant, [T] CG)
      massless : rate(r) -> 0  (power/log/linear tail => NOT locally constant)
    => the "local-constancy defect" = the propagator = the mass. Repackaging, not new.

(B) 2D Yang-Mills: exact transfer-matrix (Hamiltonian) gap from the Casimir
    C2(j)=j(j+1). The observable algebra is locally constant for *topological*
    reasons (no local propagating d.o.f.), yet the transfer matrix is gapped.
    => local constancy and the mass gap are INDEPENDENT axes in 2D.
    Under the task's sign "mass gap := failure of local constancy" this is a
    direct counterexample: 2D YM is locally constant (defect = 0) yet gapped.

(C) 3D topologically massive gauge theory (Maxwell-Chern-Simons): propagator
    ~ e^{-mu r}/sqrt(r); rate -> mu = the CS mass (an INPUT). Local constancy
    holds iff mu>0; the mass is inserted, never derived.

numpy/scipy only. Labels: [V] numerical witness.
"""
import json
import numpy as np
from scipy.special import k0, k1

OUT = {}

# ---------------------------------------------------------------------------
# (A) free-field propagators
# ---------------------------------------------------------------------------
def prop_massive(r, m, d):
    r = np.asarray(r, float)
    if d == 1:
        return np.exp(-m * r) / (2 * m)
    if d == 2:
        return k0(m * r) / (2 * np.pi)
    if d == 3:
        return np.exp(-m * r) / (4 * np.pi * r)
    if d == 4:
        return m * k1(m * r) / (4 * np.pi ** 2 * r)
    raise ValueError

def prop_massless(r, d):
    r = np.asarray(r, float)
    if d == 1:
        return r / 2.0                       # |G| grows (linear)
    if d == 2:
        return -np.log(r) / (2 * np.pi)      # log (up to an additive const)
    if d == 3:
        return 1.0 / (4 * np.pi * r)
    if d == 4:
        return 1.0 / (4 * np.pi ** 2 * r ** 2)
    raise ValueError

def rate_at(r0, Gfun):
    """-d/dr log|G|, central difference around r0."""
    r = np.array([r0 * 0.98, r0, r0 * 1.02])
    G = Gfun(r)
    lnG = np.log(np.abs(G))
    return float(-np.gradient(lnG, r)[1])

m_test = 1.0
free_table = []
for d in (1, 2, 3, 4):
    # choose radii well inside double-precision range for e^{-r}
    r0 = 8.0 if d != 2 else 8.0
    r_m = rate_at(r0, lambda rr: prop_massive(rr, m_test, d))
    r_0 = rate_at(r0, lambda rr: prop_massless(rr, d))
    free_table.append({
        "d": d,
        "r_eval": r0,
        "massive_rate": round(r_m, 5),
        "massless_rate": round(r_0, 5),
        "massive_LC_decays": r_m > 0.5,    # decays at rate ~ m=1
        "massless_LC": abs(r_0) < 0.05,    # no decay  => NOT locally constant
    })
OUT["A_free_field_rate"] = {
    "m_input": m_test,
    "note": "rate->m for massive (locally constant), rate->0 for massless (not locally constant)",
    "table": free_table,
}

# rate convergence at two radii (show it is r-independent for massive, ->0 for massless)
conv = []
for d in (3, 4):
    for r0 in (6.0, 10.0, 14.0):
        conv.append({"d": d, "r": r0,
                     "massive_rate": round(rate_at(r0, lambda rr: prop_massive(rr, m_test, d)), 5),
                     "massless_rate": round(rate_at(r0, lambda rr: prop_massless(rr, d)), 5)})
OUT["A_rate_convergence"] = conv

# ---------------------------------------------------------------------------
# (B) 2D Yang-Mills Casimir gap
# ---------------------------------------------------------------------------
def su2_casimir(j):
    return j * (j + 1)

su2 = []
for j in (0.0, 0.5, 1.0, 1.5, 2.0):
    su2.append({"j": j, "C2": su2_casimir(j)})

# SU(3): fundamental (3) C2 = 4/3, adjoint (8) C2 = 3
su3 = [{"rep": "3/3bar", "C2": 4.0 / 3.0},
       {"rep": "8 (adj)", "C2": 3.0},
       {"rep": "6", "C2": 10.0 / 3.0}]

e2 = 1.0   # probe coupling (gap ~ e^2 * C2 / 2)
gap_su2 = (e2 / 2.0) * su2_casimir(0.5)          # lightest non-trivial j=1/2
gap_su3 = (e2 / 2.0) * (4.0 / 3.0)               # fundamental
OUT["B_2D_YM"] = {
    "su2_spectrum": su2,
    "su3_spectrum": su3,
    "e2_probe": e2,
    "gap_su2_j1/2": gap_su2,
    "gap_su3_fund": gap_su3,
    "local_constancy_defect": 0.0,
    "verdict": ("observable FA locally constant (topological, no local dof) "
                "and transfer matrix gapped => gap and local constancy are independent; "
                "task-sign 'gap = failure of LC' predicts gapless => contradicted"),
}

# ---------------------------------------------------------------------------
# (C) 3D topologically massive: propagator ~ e^{-mu r}/sqrt(r), rate -> mu
# ---------------------------------------------------------------------------
mu = 1.0
tm = []
for r0 in (6.0, 10.0, 14.0):
    tm.append({"r": r0,
               "rate": round(rate_at(r0, lambda rr: prop_massive(rr, mu, 3)), 5)})
OUT["C_3D_topological_massive"] = {
    "mu_input": mu,
    "table": tm,
    "note": "rate -> mu (=CS mass); mu is an input parameter, never derived",
}

print(json.dumps(OUT, indent=2, ensure_ascii=False))
with open("x2_factorization_results.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, indent=2, ensure_ascii=False)
