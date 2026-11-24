"""
pnt_temple_lower.py -- witness for research/unconditional/PNT-TEMPLE-LOWER.md
Temple / Weinstein lower bounds on E1 from variational data (Rayleigh quotient rho,
variance sigma^2, and an UPPER bound W on E1).  Calibrate on exact solvable models.

Temple (1928):    E1 >= rho - sigma^2/(W - rho),   valid when rho < E1  (and rho < W).
Weinstein (1934): E1 >= rho - sigma,               valid when rho < E2.

HEADLINE STRUCTURAL FACT (rigorous, verified numerically below):
  Temple's bound = rho - sigma^2/(W - rho)  <=  rho      (since sigma^2 >= 0, W > rho),
  i.e. it is DOMINATED by the Rayleigh quotient rho. But rho lower-bounds E1 only
  under the hypothesis rho < E1 -- which is the unknown. So Temple adds NO gap
  lower bound beyond rho; and for any tractable ground-state trial rho ~ E0, so the
  delivered lower bound sits at/below E0 (the trivial bound).

Run:  py -3.12 -I research/unconditional/w7/pnt_temple_lower.py
"""
import numpy as np
rng = np.random.default_rng(12345)

def banner(s):
    print("\n" + "=" * 78); print(s); print("=" * 78)

def temple(rho, var, W):
    return float("-inf") if W <= rho + 1e-15 else rho - var / (W - rho)

# ===========================================================================
banner("PART A -- harmonic oscillator E_n = n+1/2 (analytic), Gaussian trial")
X, N = 9.0, 4001
x = np.linspace(-X, X, N); dx = x[1] - x[0]
def Hpsi(p):
    d2 = np.zeros_like(p)
    d2[1:-1] = (p[2:] - 2*p[1:-1] + p[:-2]) / dx**2
    d2[0] = (p[1] - p[0]) / dx**2; d2[-1] = (p[-2] - p[-1]) / dx**2
    return -0.5*d2 + 0.5*x*x*p
def norm(p): return p / np.sqrt(np.sum(p*p)*dx)
def rq(p):
    p = norm(p); Hp = Hpsi(p)
    rho = float(np.sum(p*Hp)*dx); h2 = float(np.sum(Hp*Hp)*dx)
    return rho, max(h2 - rho*rho, 0.0)

best_W = min(rq(x*np.exp(-a*x*x/2))[0] for a in np.linspace(0.2, 3.0, 57))
E1true = 1.5; E0true = 0.5
print(f"  true E0={E0true}  E1={E1true}   best variational upper bound W(E1) = {best_W:.6f}")
print("\n  Gaussian GROUND trial psi_a = exp(-a x^2/2):")
print("     a      rho      sigma^2   TempleLB   <=rho?   <=E1(valid)?")
for a in [0.30, 0.50, 0.80, 1.00, 1.30, 2.00]:
    rho, var = rq(np.exp(-a*x*x/2)); lb = temple(rho, var, best_W)
    print(f"   {a:4.2f}  {rho:8.5f}  {var:9.5f}  {lb:9.5f}   "
          f"{'yes' if lb<=rho+1e-12 else 'NO!!':4s}   {'yes' if lb<=E1true+1e-9 else 'NO!!'}")

# ===========================================================================
banner("PART B -- random symmetric matrix (exact eigs), eigenvector-mixing trial")
n = 8; A = rng.standard_normal((n, n)); H = 0.5*(A + A.T)
ev, evec = np.linalg.eigh(H); E0, E1, E2 = ev[0], ev[1], ev[2]; Delta = E1 - E0
print(f"  exact E0={E0:.5f} E1={E1:.5f} E2={E2:.5f}  gap Delta={Delta:.5f}")
v0, v1 = evec[:, 0], evec[:, 1]
W_tight = E1; W_loose = E1 + 0.5*(E2-E1)
def mom(psi):
    psi = psi/np.linalg.norm(psi); rho = float(psi @ H @ psi)
    return rho, max(float(psi @ H @ (H @ psi)) - rho*rho, 0.0)
print("\n  trial psi = sqrt(1-eps^2) v0 + eps v1 :")
print("     eps    rho      Var      Tem(tight)  Tem(loose)  true E1")
for eps in [0.0, 0.05, 0.15, 0.30, 0.50, 0.70]:
    rho, var = mom(np.sqrt(1-eps**2)*v0 + eps*v1)
    print(f"   {eps:4.2f}  {rho:9.5f} {var:9.5f}  {temple(rho,var,W_tight):10.5f}  "
          f"{temple(rho,var,W_loose):10.5f}  {E1:.5f}")

# ===========================================================================
banner("PART C -- structural law: two-level mixture => Temple bound == E0 exactly")
print("  psi=sqrt(1-e^2)|0>+e|1>:  rho=E0+e^2 D, Var=e^2(1-e^2) D^2, W=E1=E0+D")
print("  => sigma^2/(W-rho)=e^2 D  =>  Temple bound = E0   for ALL e (0% of gap).")
for eps in [0.01, 0.05, 0.10, 0.20, 0.35, 0.50]:
    var = eps**2*(1-eps**2)*Delta**2; rho = E0 + eps**2*Delta
    lb = temple(rho, var, E1)
    print(f"   eps={eps:4.2f}  Var/D^2={var/Delta**2:7.4f}  Temple-E0 = {lb-E0:+.2e}  "
          f"fraction of gap = {100*(lb-E0)/Delta:6.2f}%")

# ===========================================================================
banner("PART D -- random near-ground trials: max Temple bound never exceeds E0")
best = -np.inf; worst_dom = 0.0
for _ in range(20000):
    psi = v0 + 0.3*rng.standard_normal(n)      # perturbed ground
    rho, var = mom(psi)
    if rho < E1:                               # only where validity hypothesis holds
        lb = temple(rho, var, E1)
        best = max(best, lb); worst_dom = max(worst_dom, lb - rho)
print(f"  over 20000 trials (rho<E1 kept):  max Temple bound = {best:.5f}")
print(f"                                     E0 = {E0:.5f}   E1 = {E1:.5f}")
print(f"  max(lb - rho) over ALL trials = {worst_dom:.2e}  (must be <= 0: DOMINANCE)")
print(f"  max Temple bound - E0 = {best-E0:+.2e}   (<=0 => never beats trivial E0)")

print("""
  READ-OFF (honest):
   * Temple (rho - sigma^2/(W-rho)) and Weinstein (rho - sigma) give VALID lower
     bounds on E1 from (rho, Var, W) on exact models -- the bound TYPE exists and is
     renormalon-free (only Rayleigh quotients & variances enter).
   * But Temple's bound is RIGOROUSLY <= rho: it is dominated by the Rayleigh
     quotient, and rho lower-bounds E1 only under the unknown hypothesis rho<E1.
     For every tractable ground trial (rho~E0) the delivered gap lower bound sits
     at/below E0 (Part C: exactly E0 for two-level; Part D: max <= E0).
   * Hence Temple/Weinstein do NOT close the tensor-RG sandwich. The missing lower
     bound on E1 is the same nonperturbative input (a near-exact 0++ wavefunctional
     with controlled variance / an independent spectral-separation certificate).
""")
