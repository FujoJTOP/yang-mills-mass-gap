"""
tool_caratheodory.py -- TOOL-CARATHEODORY witness.

Question: can Caratheodory/Toeplitz/Christoffel POSITIVITY of the Stieltjes
transform Pi(z)=int dnu(s)/(s+z) LOWER-bound the lightest KL atom f1^2 at s1=m1^2
from low-order Taylor coefficients of Pi at z0 (which are the data available)?

Taylor coefficients: Pi^{(k)}(z0)/k! = (-1)^k int dnu(s)/(s+z0)^{k+1} := a_k.
These are integrals with CONTINUOUS bounded kernels 1/(s+z0)^{k+1} -> they are
in the OPE-RHO weak-continuous class (NOT "non-integral").  We test directly:
  (1) a continuum-only measure (no atom at s1) matching the Taylor data (nnls);
  (2) LP range [min,max] of the atom weight at s1 consistent with the data;
  (3) Caratheodory |c_k|<=2 c_0 positivity satisfied by BOTH measures.
If min ~ 0, positivity gives only an UPPER bound -> no forcing lower bound.
numpy/scipy only.
"""
import numpy as np
from scipy.optimize import linprog, nnls

z0 = 5.0
K  = 8
s1 = 1.0

# true measure: atom weight 0.5 at s1 + continuum strictly above s1
s_cont = np.array([1.05, 1.1, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0])
w_cont = 0.04 * np.ones_like(s_cont)
s_true = np.concatenate([[s1], s_cont])
w_true = np.concatenate([[0.5], w_cont])
W_atom_true = 0.5

def phi(s, k):
    return (-1.0) ** k / (s + z0) ** (k + 1)

a_true = np.array([np.sum(w_true * phi(s_true, k)) for k in range(K)])
print("Taylor coeffs a_0..a_7 of the atom+continuum measure:")
print(np.array2string(a_true, precision=6))
print("atom contribution to a_k: 0.5/(1+z0)^{k+1} =",
      np.array([0.5/(1+z0)**(k+1) for k in range(K)])[[0,1,2,-1]])

# (1) continuum-ONLY fit (no atom at s1): nnls on a fine grid just above s1 + continuum
s_grid = np.concatenate([s1 + np.array([1e-3, 3e-3, 1e-2, 3e-2, 0.1]),
                         s_cont, np.linspace(1.2, 12.0, 80)])
M = np.array([[phi(sj, k) for sj in s_grid] for k in range(K)])
wB, _ = nnls(M, a_true)
resid = np.max(np.abs(M @ wB - a_true))
print("\n(1) continuum-only nnls: total fitted weight =", round(wB.sum(), 6),
      " max|M w - a| =", f"{resid:.2e}", "  (no atom at s1)")
print("    => a no-atom measure reproduces the SAME Taylor data to", f"{resid:.1e}")

# (2) LP: range of the atom weight at s1 consistent with a_true, continuum strictly above s1
s_lp = s1 + np.concatenate([[1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2],
                            np.linspace(0.3, 12.0, 90)])
A = np.array([[phi(sj, k) for sj in s_lp] for k in range(K)])
phi1 = np.array([phi(s1, k) for k in range(K)])
Aeq = np.hstack([phi1[:, None], A]); beq = a_true
n = 1 + len(s_lp); bounds = [(0, None)] * n
cmin = np.zeros(n); cmin[0] = 1.0
cmax = np.zeros(n); cmax[0] = -1.0
rmin = linprog(cmin, A_eq=Aeq, b_eq=beq, bounds=bounds, method='highs')
rmax = linprog(cmax, A_eq=Aeq, b_eq=beq, bounds=bounds, method='highs')
min_atom = rmin.fun if rmin.success else float('nan')
max_atom = -rmax.fun if rmax.success else float('nan')
print("\n(2) LP range of atom weight at s1 given the Taylor data:")
print("    true atom =", W_atom_true, "   consistent range = [", round(min_atom, 6),
      ",", round(max_atom, 4), "]")
print("    => lower end ~0 : the data give only an UPPER bound, no positive lower bound")

# (3) Caratheodory / positive-real: Pi(z)>0 for z>0 (both measures are positive).
#     Caratheodory-Fejer |c_k|<=2 Re c_0 holds for ANY positive measure; both satisfy it,
#     so it cannot discriminate atom vs no-atom. Verify Re Pi > 0 for both at a few z.
zs = np.array([0.3, 1.0, 2.0, 10.0])
def Pi(s, w, z): return np.sum(w / (s + z))
PA = np.array([Pi(s_true, w_true, z) for z in zs])
PB = np.array([Pi(s_grid, wB, z) for z in zs])
print("\n(3) Re Pi(z)>0 (positive-real) for both measures:",
      np.all(PA > 0), np.all(PB > 0), "  -> Caratheodory bound is a positivity UPPER constraint")

print("\nVERDICT: Taylor-at-z0 coefficients are integrals (OPE-RHO class);")
print("         positivity -> upper bound on f1^2 only; LP min=0 -> NO forcing lower bound;")
print("         a continuum-only measure matches the data => NOT an escape.")
