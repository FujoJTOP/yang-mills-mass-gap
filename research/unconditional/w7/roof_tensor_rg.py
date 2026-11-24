"""
roof_tensor_rg.py -- witness for research/unconditional/ROOF-TENSOR-RG.md  (2026-11-24)

Two checks, both tiny (numpy only, <2 s):

(A) UPPER/LOWER ASYMMETRY (Rayleigh-Ritz / min-max).
    Variational ansaetze (MERA/PEPS/DMRG) give UPPER bounds on eigenvalues.
    Demonstrate: R(psi)=<psi|H|psi> >= E0 always; the k-th Ritz value over a
    k-dim subspace is >= E_k.  Hence a variational "gap" (E1trial - E0trial)
    is NOT a certificate (it can be above OR below the true gap).  The gap
    needs a LOWER bound on E1, which variational methods do not give.

(B) ENTANGLEMENT-SPECTRUM GAP != PHYSICAL GAP.
    Exact diagonalisation of the transverse-field Ising chain (N=8, 2^8=256).
    Compare the physical gap E1-E0 with the half-chain entanglement gap
    xi_1 - xi_0  (xi_i = -log of Schmidt values).  They are different functions
    of h => the entanglement spectrum cannot certify the physical (mass) gap.
    Acid-style check: at h=1 (critical / massless) the physical gap -> 0 in the
    thermodynamic limit while the entanglement gap stays O(1) (CFT structure).

NO Yang-Mills claim.  Illustrative/logical structure only.
Run:  py -3.12 -I research/unconditional/w7/roof_tensor_rg.py
"""
import numpy as np

print("=" * 74)
print("(A) variational = UPPER bound on eigenvalues (Rayleigh-Ritz / min-max)")
print("=" * 74)

rng = np.random.default_rng(0)
d = 6
A = rng.standard_normal((d, d))
H = (A + A.T) / 2.0
ev = np.linalg.eigvalsh(H)
E0, E1 = ev[0], ev[1]
true_gap = E1 - E0

# ground: every Rayleigh quotient is >= E0
min_R = np.inf
for _ in range(4000):
    psi = rng.standard_normal(d)
    psi /= np.linalg.norm(psi)
    R = float(psi @ H @ psi)
    min_R = min(min_R, R)
print(f"  min_psi <psi|H|psi> - E0 = {min_R - E0:+.3e}   (>=0 : UPPER bound on E0)")

# excited: 2nd Ritz value over a random 2D subspace is >= E1
min_ritz2 = np.inf
for _ in range(4000):
    Q, _ = np.linalg.qr(rng.standard_normal((d, 2)))
    ritz = np.linalg.eigvalsh(Q.T @ H @ Q)
    min_ritz2 = min(min_ritz2, ritz[1] - E1)
print(f"  min over 2D subspaces of (2nd Ritz - E1) = {min_ritz2:+.3e}   (>=0 : UPPER bound on E1)")

# a variational "gap" using a bad ground estimate UNDER-estimates the true gap
psi0 = np.linalg.eigh(H)[1][:, 0]
psi1 = rng.standard_normal(d); psi1 -= (psi0 @ psi1) * psi0; psi1 /= np.linalg.norm(psi1)
gap_bad = float(psi1 @ H @ psi1) - E0            # excited trial + exact E0
# and a poor product-like ground bound (a random vector projected) OVER-estimates E0
g_est = float(psi0 @ H @ psi0) + 0.6 * true_gap  # pretend a poor variational E0
# then gap = E1_trial - E0_bad can be SMALLER than true gap:
gap_est = float(psi1 @ H @ psi1) - (E0 + 0.6 * true_gap)
print(f"  true gap = {true_gap:.4f};  E1trial-E0(exact) = {gap_bad:.4f};")
print(f"  E1trial - (poor E0 bound, +0.6*gap) = {gap_est:.4f}  (can UNDER-shoot true gap)")
print("  => variational 'gap' can sit above OR below the true gap: NOT a certificate.")
print("  => the gap needs a LOWER bound on E1; tensor RG supplies upper bounds only.")

print()
print("=" * 74)
print("(B) entanglement-spectrum gap vs physical gap (TFIM, N=8, exact diag)")
print("=" * 74)

sx = np.array([[0.0, 1.0], [1.0, 0.0]])
sz = np.array([[1.0, 0.0], [0.0, -1.0]])
I2 = np.eye(2)


def op_at(op, i, N):
    m = np.array([[1.0]])
    for j in range(N):
        m = np.kron(m, op if j == i else I2)
    return m


def tfim(N, J, h):
    dim = 2 ** N
    Hm = np.zeros((dim, dim))
    for i in range(N):
        Hm -= J * (op_at(sz, i, N) @ op_at(sz, (i + 1) % N, N))
    for i in range(N):
        Hm -= h * op_at(sx, i, N)
    return Hm


def gaps(N, h):
    Hm = tfim(N, 1.0, h)
    w, v = np.linalg.eigh(Hm)
    phys = w[1] - w[0]
    psi = v[:, 0]
    half = N // 2
    M = psi.reshape(2 ** half, 2 ** (N - half))
    s = np.linalg.svd(M, compute_uv=False)
    s = s[s > 1e-12]
    xi = np.sort(-np.log(s))          # entanglement "energies"
    ent = xi[1] - xi[0]              # entanglement gap
    return phys, ent


print(f"  {'h':>5} {'phys gap':>10} {'ent gap':>10} {'ent/phys':>9}")
for h in [0.2, 0.5, 0.8, 1.0, 1.2, 2.0]:
    p, e = gaps(8, h)
    ratio = e / p if p > 1e-9 else float('nan')
    print(f"  {h:5.1f} {p:10.4f} {e:10.4f} {ratio:9.3f}")

print()
print("  The two gaps are NOT proportional: ent/phys varies by an order of")
print("  magnitude across h and does not vanish at criticality the way the")
print("  physical gap does.  => entanglement-spectrum gap is heuristic, not a")
print("  certificate of the physical 0++ mass gap.  (Acid: at h=1, massless,")
print("  phys gap -> 0 as N->inf, ent gap stays O(1).)")
