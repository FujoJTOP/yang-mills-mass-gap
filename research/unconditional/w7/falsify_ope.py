"""
falsify_ope.py -- witness for FALSIFY-OPE.
Question: can OPE coefficients of tr F^2 (operator-specific, perturbative-UV) FORCE a lower bound
on the 0++ KL atom f1^2?

Structural fact: Pi(Q^2)=int dnu(s)/(Q^2+s) = sum_k (-1)^k M_{-(k+1)}/Q^{2k} + ...  So the OPE
coefficients of the 1/Q^{2k} terms ARE the moments M_{-k}=int s^{-k} dnu -- LINEAR functionals of nu.
Demonstrate: an ATOMLESS narrow bump at the floor (mass = atom weight, width h) reproduces every
finite OPE moment to O(h) -> 0.  Hence finite OPE cannot force f1^2>0 (the atom is hideable).

Runs: py -3.12 -I  (numpy only)
"""
import numpy as np

s1 = 1.0
f1 = 0.5
cont = np.linspace(2.0, 4.0, 400); cont_w = np.full_like(cont, 0.5 / len(cont))
K = 8
ks = np.arange(1, K + 1, dtype=float)

def moms(nodes, w):
    nodes = np.asarray(nodes, float); w = np.asarray(w, float)
    return np.array([np.sum(w * nodes ** (-k)) for k in ks])

M_atom  = f1 * s1 ** (-ks) + moms(cont, cont_w)          # reference: atom 0.5 at s1 + continuum
M_atomless_base = moms(cont, cont_w)                       # continuum only (no floor mass)

print("reference OPE moments M_{-1}..M_{-%d}:" % K); print("  ", np.round(M_atom, 6))

# atomless replacement of the atom: uniform bump on [s1, s1+2h], total mass f1, NO point mass.
print("\natomless bump (mass f1=0.5, width h) vs atom -- max |Delta M_{-k}| over k=1..%d:" % K)
for h in [1e-1, 1e-2, 1e-3, 1e-4]:
    b = np.linspace(s1, s1 + 2 * h, 2000); bw = np.full_like(b, f1 / len(b))
    M_bump = moms(b, bw) + M_atomless_base
    d = np.max(np.abs(M_bump - M_atom))
    print("  h=%.0e   max|DeltaM| = %.3e   (atom weight forced from OPE: 0)" % (h, d))

print("\n=> as h->0 the atomless bump matches ALL finite OPE moments. So for any finite OPE truncation")
print("   an atomless measure (same floor s1, zero weight at s1) is feasible => OPE cannot force f1^2>0.")
print("   The floor is unchanged (inf supp = s1 in both); only the POINT weight differs.")

# trace anomaly A = M_{-2} is one of these moments (a normalization), not a pole-forcing condition.
print("\ntrace anomaly A = M_{-2} = %.6f  -- already a moment above; an atomless bump matches it too." % M_atom[1])

print("\nThe atom is recovered ONLY by the FULL expansion (all k; equivalently full Pi(Q^2), every Q^2):")
print("  infinite-OPE / Stieltjes inversion = J->infinity = H1; and the perturbative coefficients")
print("  themselves carry the u=2 IR renormalon (ROOF-OPE) = H1.  Operator-specific structure does not")
print("  help: Pi is LINEAR in nu, and the atom is hideable behind a continuum matching any finite OPE.")
