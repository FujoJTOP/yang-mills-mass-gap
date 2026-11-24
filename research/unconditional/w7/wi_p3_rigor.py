"""
wi_p3_rigor.py -- witness for research/unconditional/WI-P3-RIGOR.md
Effective mass m_O(tau) = -d/dtau log C_O(tau), C_O(tau)=sum_n |c_n|^2 e^{-e_n tau}.

Tests the approach law m(tau) - e1:
  (1) isolated atom  -> EXPONENTIAL approach, rate Delta;  T(tau):= -d log|m-e1|/dtau -> Delta (const)
  (2) continuum edge rho~(s-E)^p -> POWER-LAW approach  m-E ~ (p+1)/tau;  T(tau) ~ 1/tau -> 0
  (3) finite-L: N discrete levels, spacing delta -> m-e1 LOOKS exponential with rate delta,
      and delta -> 0 as N -> infinity (so single finite L cannot discriminate).

Lightweight numpy only.  run:  py -3.12 -I wi_p3_rigor.py
"""
import numpy as np

def m_eff(s, w, tau):
    """Numerically stable: factor out e^{-s0 tau} so large tau does not underflow."""
    s = np.asarray(s, float); w = np.asarray(w, float); tau = np.asarray(tau, float)
    s0 = s.min()
    E = np.exp(-np.outer(tau, s - s0))       # (ntau, ns), exponents <= 0
    num = E @ (w * (s - s0)); den = E @ w
    return s0 + num / den

def Tlocal(tau, dm):
    """T(tau) = -d log|dm|/dtau  (local log-derivative of the approach)."""
    l = np.log(np.abs(dm))
    return -np.gradient(l, tau)

taus = np.linspace(0.5, 30.0, 3000)

print("="*74)
print("(1) ISOLATED ATOM: two levels e1=1, e2=2, equal weights")
print("="*74)
m1 = m_eff([1.0, 2.0], [1.0, 1.0], taus); dm1 = m1 - 1.0
T1 = Tlocal(taus, dm1)
for t in (2.0, 4.0, 8.0):
    i = np.searchsorted(taus, t)
    print(f"  tau={t:4.1f}:  m-e1={dm1[i]:.3e}   T={T1[i]:.4f}   (expect T -> Delta = 1)")

print()
print("="*74)
print("(2) CONTINUUM EDGE: rho ~ (s-2)^1 on [2,3]  (p=1) -> m-E ~ (p+1)/tau")
print("="*74)
sg = np.linspace(2.0, 3.0, 40001); wg = (sg - 2.0); wg = wg / wg.sum()
m2 = m_eff(sg, wg, taus); dm2 = m2 - 2.0
T2 = Tlocal(taus, dm2)
for t in (5.0, 15.0, 30.0):
    i = np.searchsorted(taus, t)
    print(f"  tau={t:5.1f}:  m-E={dm2[i]:.3e} (~(p+1)/tau={2.0/t:.3e})   T={T2[i]:.4f} (~1/tau={1.0/t:.4f})")

print()
print("="*74)
print("(3) FINITE-L CONTINUUM: N discrete levels in [2,3], spacing delta=1/(N-1)")
print("    late-tau effective mass looks EXPONENTIAL with rate ~ delta; delta -> 0 as N grows")
print("="*74)
print(f"  {'N':>5} {'delta':>9} {'T(tau~3/delta)':>16}   (expect ~ delta)")
for N in [11, 21, 51, 101, 201, 401]:
    delta = 1.0 / (N - 1)
    sN = np.linspace(2.0, 3.0, N); wN = np.ones(N)      # (uniform weights; rate result is density-independent)
    tau0 = 3.0 / delta                                   # ~3 level-spacings: only 2 lowest levels survive
    tg = np.array([tau0 * 0.9, tau0, tau0 * 1.1])
    mN = m_eff(sN, wN, tg); TN = Tlocal(tg, mN - 2.0)
    print(f"  {N:5d} {delta:9.5f} {TN[1]:16.5f}")

print()
print("READING:")
print("  atom    -> T(tau) -> const (Delta)          [exponential approach]")
print("  cut     -> T(tau) ~ 1/tau -> 0               [power-law approach]")
print("  finite L-> T ~ delta_N (the level spacing),  delta_N -> 0 as N=L^3 grows.")
print("  => a SINGLE finite L cannot tell atom from cut (both look exponential, rate delta_L);")
print("     the discriminating content is the L-scaling of delta_L (atom: Delta>0; cut: ->0).")
