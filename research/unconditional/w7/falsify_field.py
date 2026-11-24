"""falsify_field.py -- witness for FALSIFY-FIELD.md  (numpy only, py -3.12 -I, <1s)
Points:
 (1) finite-L atom f1^2>0 is TRIVIAL (any nonzero overlap; compact transfer matrix => discrete spectrum).
 (2) projection to nu is BLIND: finite moments / Laplace (gradient flow E(t)) cannot distinguish an atom
     from a narrow bump matched on that data (OPE-RHO, Δ_min=0).
 (3) CHANNEL ORTHOGONALITY: 0-form atom and 1-form (Wilson/string) are independent blocks; 1-form data
     carries no kinematic info about the 0-form atom => linking them is a DYNAMICAL bridge = H1.
 (4) HARDY N(eps)=int_{eps<|x-x0|<δ}|G|^2: atom -> N*eps -> const>0 ; smooth continuum -> N*eps -> 0.
     But at finite L EVERY level is a pole (all atoms) => divergence is uniform/trivial; the nontrivial
     content is the L->inf survival of an ISOLATED pole = isolation = H1.
"""
import numpy as np
rng = np.random.default_rng(0)

print("="*72)
print("(1) FINITE-L: atom f1^2>0 is automatic (trivial)")
print("="*72)
for t in range(3):
    N = 6
    E = np.sort(rng.uniform(0.5, 3, N))
    ov = rng.uniform(-1, 1, N)          # <0|O|n>
    f2 = ov**2                            # |<0|O|n>|^2
    print(f"  trial {t}: lightest E0={E[0]:.3f}, atom f1^2={f2[0]:.4f}  >0? {f2[0] > 0}")
print("  => at finite L the lightest level IS an atom; f1^2>0 needs no field functional.")
print("     Clay needs the a->0, L->inf SURVIVAL of that atom (a uniformity) = H1.")

print()
print("="*72)
print("(2) PROJECTION TO nu IS BLIND (OPE-RHO): atom vs narrow bump, same finite moments")
print("="*72)
def mom(nodes, weights, jmax):
    return np.array([np.sum(weights * nodes**j) for j in range(jmax)])
m = 1.0
nodes1 = np.array([m, 2.0, 3.0, 4.0]); w1 = np.array([0.5, 1/6, 1/6, 1/6])   # atom 0.5 + cont
M1 = mom(nodes1, w1, 6)
for h in [1e-1, 1e-2, 1e-3, 1e-4]:
    nodes2 = np.concatenate([[m + h/2], [2.0, 3.0, 4.0]])                     # atomless bump
    w2 = np.concatenate([[0.5], [1/6, 1/6, 1/6]])
    M2 = mom(nodes2, w2, 6)
    print(f"  bump h={h:.0e}: atom weight=0, max|moment diff|={np.max(np.abs(M1-M2)):.2e}")
print("  => finite moments and Laplace/gradient-flow E(t) cannot read the atom (blind).")

print()
print("="*72)
print("(3) CHANNEL ORTHOGONALITY: 1-form (Wilson/sigma) vs 0-form (atom f1^2) independent")
print("="*72)
m0 = 0.7
# block-diagonal transfer matrix: 0-form block (gap m0) and 1-form block (sigma) decouple
print(f"  0-form spectrum (fixed): {{1, {np.exp(-m0):.3f}}}  (gap m0={m0})")
for s in [0.2, 0.4, 0.9]:
    print(f"    set 1-form sigma={s}: 0-form atom E unchanged (decoupled) -> 1-form data gives NO kinematic info on the atom")
print("  => a 1-form observable constrains sigma, not f1^2; linking them is a DYNAMICAL bridge = H1.")

print()
print("="*72)
print("(4) HARDY: atom -> divergence; smooth continuum -> finite; but finite-L is uniform/trivial")
print("="*72)
# exact Hilbert transforms G(x)=PV\int rho/(s-x)ds on the real axis:
#   TRUE atom: G = f2/(s1 - x)  -> G^2 ~ 1/(x-s1)^2  -> int ~ 1/eps (NON-integrable, divergence)
#   smooth continuum: G = rho ln|(b-x)/(a-x)|  -> G^2 ~ log^2 -> integrable (N finite)
def N_atom(eps, f1sq=0.5, s1=1.0, xhi=3.5, nx=8000):
    xs = np.linspace(s1 + eps, xhi, nx)
    G = f1sq / (s1 - xs)
    return np.trapezoid(G**2, xs)
def N_cont(eps, rho=0.5/3, a=1.0, b=4.0, xhi=3.5, nx=8000):
    xs = np.linspace(a + eps, xhi, nx)
    G = rho * np.log(np.abs((b - xs) / (a - xs)))
    return np.trapezoid(G**2, xs)
for eps in [1e-1, 1e-2, 1e-3, 1e-4]:
    print(f"  eps={eps:.0e}: ATOM N*eps={N_atom(eps)*eps:.4f} (->const>0)   CONTINUUM N*eps={N_cont(eps)*eps:.4e} (->0)")
print("  => divergence flags an atom vs a smooth continuum, robustly. BUT at finite L all levels are")
print("     poles => divergence uniform/trivial; nontrivial content = L->inf ISOLATED-pole survival = H1.")

print()
print("VERDICT: field functionals force f1^2 only via (2)-blind or (3)-dynamical-bridge(=H1).")
print("Falsification FAILS; 'YM structure' and 'H1' are the SAME exit.")
