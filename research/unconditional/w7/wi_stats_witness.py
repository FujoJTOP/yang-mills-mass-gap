"""
wi_stats_witness.py -- WI-STATS statistical-feasibility witness (fast).

See report WI-STATS.md. Toy model on real-tau correlator, T points:
  atom: C(t) = w0 + f2*exp(-m1 t)              (f2 = operator overlap)
  cut : C(t) = w0 + Int_mc^inf exp(-(E-mc)/D)/D * exp(-E t) dE
Operator-reweighting = different f2 (atom) / different D (cut).
Estimator: profile-likelihood m from w0 + f2*exp(-m t) (minimize_scalar).
Discriminant: spread S = std_k(m_hat_k) over K operators.
Run:  py -3.12 -I research/unconditional/w7/wi_stats_witness.py
"""
import numpy as np
from scipy.optimize import minimize_scalar
rng = np.random.default_rng(0)

T = 6
w0 = 1.0
m1 = 0.8
tau = np.arange(T)
atom_f2 = [0.3, 1.0, 5.0]
cut_mc = 1.0
cut_Delta = [0.5, 1.0, 2.0]
E = np.linspace(cut_mc, cut_mc + 60.0, 6001); dE = E[1] - E[0]

def C_atom(f2): return w0 + f2 * np.exp(-m1 * tau)
def C_cut(D):
    rho = np.exp(-(E - cut_mc) / D) / D
    return w0 + np.array([np.sum(rho * np.exp(-E * t)) * dE for t in tau])

C_atoms = [C_atom(f) for f in atom_f2]
C_cuts  = [C_cut(D) for D in cut_Delta]

def fit_m(C, sig):
    Wt = 1.0 / sig**2
    def f(m):
        X = np.column_stack([np.ones(T), np.exp(-m * tau)])
        XtW = X.T * Wt
        try:
            beta = np.linalg.solve(XtW @ X, XtW @ C)
        except np.linalg.LinAlgError:
            return 1e18
        return float(np.sum((C - X @ beta)**2 * Wt))
    return minimize_scalar(f, bounds=(0.3, 2.0), method='bounded',
                           options={'xatol': 1e-4}).x

print("=" * 64)
print("INTRINSIC effective masses (noise-free), true m1=0.8")
print("=" * 64)
aint = [fit_m(C, 1e-7) for C in C_atoms]
cint = [fit_m(C, 1e-7) for C in C_cuts]
for f, m in zip(atom_f2, aint): print(f"  atom f2={f:<4}: m_hat={m:.4f}")
for D, m in zip(cut_Delta, cint): print(f"  cut  D ={D:<4}: m_hat={m:.4f}")
print(f"  atom spread={np.std(aint):.4f} ; cut spread={np.std(cint):.4f}"
      f"  ({100*np.std(cint)/m1:.0f}% of m1)")

print("\n" + "=" * 64)
print("NOISY DISCRIMINATION vs per-point sigma0  (N ~ s_base^2/sigma0^2)")
print("=" * 64)
print(f"{'sigma0':>7} {'atom_S':>8} {'cut_S':>8} {'z':>6} {'N(s_base=1)':>12}")
ntr = 80
for sigma0 in [0.3, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002]:
    aS, cS = [], []
    for _ in range(ntr):
        ma = [fit_m(C + sigma0 * rng.standard_normal(T), sigma0) for C in C_atoms]
        mc = [fit_m(C + sigma0 * rng.standard_normal(T), sigma0) for C in C_cuts]
        aS.append(np.std(ma)); cS.append(np.std(mc))
    aS = np.array(aS); cS = np.array(cS)
    pooled = np.sqrt(0.5 * (aS.var() + cS.var()) + 1e-14)
    z = (cS.mean() - aS.mean()) / pooled
    print(f"{sigma0:>7.4f} {aS.mean():>8.4f} {cS.mean():>8.4f} {z:>6.1f} {int(round(1/sigma0**2)):>12d}")

print("\nz>=3 => per-experiment separable at ~3 sigma. sigma0~noise/sqrt(N_eff).")
print("The late-tau signal f2*exp(-m(T-1)) sets the noise ceiling; full continuation")
print("(real A needs the whole continuum) is strictly harder than this toy.")
