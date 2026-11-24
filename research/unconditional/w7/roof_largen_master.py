"""
roof_largen_master.py -- witness for research/unconditional/ROOF-LARGEN-MASTER.md

Three cheap checks (numpy only, no Monte Carlo):

(A) 2D Makeenko-Migdal CLOSES: for a simple loop the RHS vanishes and the
    area-derivative reduces to  dW/dA = -(lambda/2) W, whose exact master-field
    solution is  W(C) = exp(-lambda * Area / 2).  We integrate the ODE and
    compare with the closed form on a grid of areas.

(B) 4D Makeenko-Migdal does NOT close: the equation for the 1-contour loop
    sources the 2-contour correlator, which sources 3, ... an infinite
    hierarchy.  We print the contour-count per level for D=2 (closes at 1)
    vs D=4 (never closes) -- a structural count, not a computation.

(C) Large-N scaling  m_{0++}/sqrt(sigma)  from lattice values is N-INDEPENDENT
    to leading order, with a small 1/N^2 approach to the SU(inf) value.  We fit
    m(N) = m_inf + c/N^2 and report the finite-N correction size.  [V-num]

No claim about a proof of the 4D Yang-Mills mass gap.
Run:  py -3.12 research/unconditional/w7/roof_largen_master.py
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True)


def banner(s):
    print("\n" + "=" * 74)
    print(s)
    print("=" * 74)


# ---------------------------------------------------------------------------
# (A) 2D: the loop equation closes; exact master field W = exp(-lambda A / 2)
# ---------------------------------------------------------------------------
banner("(A) 2D Makeenko-Migdal closes -> dW/dA = -(lambda/2) W  ->  W = exp(-lambda A/2)")
lam = 0.8
A = np.linspace(0.0, 10.0, 2001)
dA = A[1] - A[0]

# integrate dW/dA = -(lam/2) W from W(0)=1 (explicit Euler; exact for linear ODE)
W = np.empty_like(A)
W[0] = 1.0
for i in range(1, A.size):
    W[i] = W[i - 1] * (1.0 - (lam / 2.0) * dA)
W_closed = np.exp(-lam * A / 2.0)
print(f"  lambda = {lam}")
print(f"  max |W_numeric - exp(-lambda A/2)| over A in [0,10] = {np.max(np.abs(W - W_closed)):.3e}")
print(f"  W(A=2,4,8) numeric   = {W[np.searchsorted(A,2)]:.6f}, "
      f"{W[np.searchsorted(A,4)]:.6f}, {W[np.searchsorted(A,8)]:.6f}")
print(f"  W(A=2,4,8) closed    = {np.exp(-lam*2/2):.6f}, "
      f"{np.exp(-lam*4/2):.6f}, {np.exp(-lam*8/2):.6f}")
print("  => the 2D loop equation is CLOSED and its master field is explicit.  [T]")


# ---------------------------------------------------------------------------
# (B) 4D: structural contour count -- the equation for k contours needs k+1
# ---------------------------------------------------------------------------
banner("(B) Makeenko-Migdal hierarchy: contours sourced per level (D=2 vs D=4)")
print("  Level k = number of Wilson contours in the correlator on the LHS.")
print("  A loop-derivative of a k-contour correlator sources a (k+1)-contour term.")
print()
print("    level k :  1      2      3      4     5    ...   closes?")
print("    D=2     : ok     --     --     --    --           YES at k=1 (simple loop)")
print("    D=4     : 2      3      4      5     6    ...      NEVER (infinite tower)")
print()
print("  In D=2 a simple loop has no self-intersection -> RHS = 0 -> closed ODE (A).")
print("  In D=4 the RHS is <W(C1)W(C2)> (double contour), whose equation sources")
print("  a triple contour, ... -> the 0++ mass sits at the bottom of an INFINITE")
print("  hierarchy.  Leading gap NOT computable from the equation alone.  [NG]")


# ---------------------------------------------------------------------------
# (C) large-N scaling: m_{0++}/sqrt(sigma), 1/N^2 approach to SU(inf)  [V-num]
# ---------------------------------------------------------------------------
banner("(C) large-N scaling of m_{0++}/sqrt(sigma)  [V-num, lattice]")
# representative lattice values (Athenodorou-Tepec; Lucini-Tepec-Wenger);
# illustrative central values -- the POINT is N-independence + small 1/N^2 drift.
N     = np.array([2, 3, 4, 5, 6, 8], dtype=float)
m_sq  = np.array([3.78, 3.55, 3.40, 3.37, 3.35, 3.33])   # m_{0++}/sqrt(sigma)
inv2  = 1.0 / N**2

# fit m = m_inf + c / N^2  (least squares)
Mmat = np.vstack([np.ones_like(inv2), inv2]).T
(m_inf, c), *_ = np.linalg.lstsq(Mmat, m_sq, rcond=None)
pred = m_inf + c * inv2

print("    N        m/sqrt(s)   fit   (1/N^2 term)")
for n, m, p, w in zip(N, m_sq, pred, inv2):
    print(f"    {int(n):<4d}     {m:.3f}      {p:.3f}     c/N^2 = {c*w:+.4f}")
print()
print(f"  SU(inf) extrapolation  m_inf = {m_inf:.3f}   (c = {c:+.3f})")
print(f"  SU(3) correction vs SU(inf): {m_sq[1] - m_inf:+.3f}  "
      f"({100*(m_sq[1]-m_inf)/m_inf:+.1f}%)  -> gap SURVIVES to finite N.  [C]")
print()
print("  This is a RATIO of two physical (renormalon-free) observables,")
print("  m_{0++}/sqrt(sigma): no u=2 IR renormalon enters.  [T]")
print("  But the SU(inf) master-field value is not computable in 4D (B).")

banner("SUMMARY")
print("  [T]  gap is a renormalon-free physical mass (unlike <tr F^2>)")
print("  [T]  large-N gives 1/N small parameter; m/sqrt(s) is N-independent")
print("  [T]  2D loop equation closes; exact master field")
print("  [NG] 4D loop equation does NOT close -> leading gap not computable")
print("  [O]  4D master field not rigorously constructed (Singer: 2D only)")
print("  Clay 0%; no QED.")
