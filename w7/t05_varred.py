"""t05_varred.py -- can variance reduction (control variate / importance sampling /
stratification) cut the strict L=2 Delta certificate cost from N~3e8 to N<=1e6?

Data (real, on disk):
  research/unconditional/w7/o24_L2A.npz   L=2, SU(2), beta=2.5, Metropolis mode A,
      per-sample kappa (K) and F (F), 15000 sweeps x 4 chains, eps=0.9, 10.02 ms/sweep.
Quoted (no re-run): w7/p02_dyndyn.log (constrained xi=0 HMC, 19500x4, 16.7 ms/sweep),
      research/v3/species/R04-L2CERT.md  (Hoeffding arithmetic: K=2/d^2=2222, N~3e8).

Identities used:
  Delta = E[kappa] - Var(F),  kappa = <v,Hess V v>,  F = <v,grad V>.
  E[F] = 0 EXACTLY (divergence/IBP identity on the fiber; o24 log E[F]=+0.0349 +- 0.05).
  => Var(F) = E[F^2] exactly, so  Delta = E[kappa - F^2]  (one-sample pairing, no plug-in mean).

Output: every number below is measured on the above samples unless tagged [T]/[C].
"""
import numpy as np

RNG = np.random.default_rng(0)


def ips_tau(x, maxtau=None):
    """Geyer initial-positive-sequence integrated autocorr time tau (sum of rho_t, t>=1)."""
    x = np.asarray(x, float)
    n = x.size
    x = x - x.mean()
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n].real
    ac /= ac[0]
    if maxtau is None:
        maxtau = n // 4
    s = 0.0
    t = 1
    while t + 1 < maxtau:
        pair = ac[t] + ac[t + 1]
        if pair < 0:
            break
        s += pair
        t += 2
    return 0.5 * s          # tau = 0.5*sum pair  (= (1+2tau_ips-1)/2)


def ips_factor(x):
    """1 + 2 tau  (variance inflation of the mean estimator)."""
    return 1.0 + 2.0 * ips_tau(x)


d = np.load(r"D:\Dev\math\research\unconditional\w7\o24_L2A.npz", allow_pickle=True)
K, F = d["K"], d["F"]                      # (15000, 4)
ns, G = K.shape
P = print

P("=" * 104)
P(f"[T05] variance-reduction audit for Delta = E[kappa]-Var(F)   data=o24_L2A (L=2,SU(2),b=2.5)  {ns}x{G}")
P("=" * 104)
P(f"  kappa: mean={K.mean():+.4f} std={K.std():.4f} max={K.max():.3f}   "
  f"F: mean={F.mean():+.4f} std={F.std():.4f} max|F|={np.abs(F).max():.3f}")
P(f"  Delta(point) = {K.mean()-F.var():+.4f}   (o24 log: +0.3942)")

# ---------------------------------------------------------------- 1. where the variance lives
P("\n[1] variance budget of the ESTIMATOR (per-sample, before 1/N)")
mu = F.mean()
z = (F - mu) ** 2                # the Var(F) channel, two-pass (as used in o13/o24/p02)
w = K - F ** 2                   # paired channel, exact E[F]=0
rows = []
for g in range(G):
    k, f = K[:, g], F[:, g]
    zg = (f - f.mean()) ** 2
    wg = k - f ** 2
    ckk = np.cov(k, zg)[0, 1]
    fK, fZ, fW = ips_factor(k), ips_factor(zg), ips_factor(wg)
    vK = k.var() * fK
    vZ = zg.var() * fZ
    vC = 2 * ckk * np.sqrt(fK * fZ)
    vW = wg.var() * ips_factor(wg)
    rows.append((vK, vZ, vC, vK + vZ + vC, vW, fK, fZ))
    P(f"  ch{g}: 1+2tau kappa={fK:6.1f}  z=(F-mu)^2={fZ:6.1f}  w=k-F^2={ips_factor(wg):6.1f}"
      f" | var_terms: k {vK:.3f}  z {vZ:.3f}  cov {vC:+.3f}  sum {vK+vZ+vC:.3f}"
      f" | paired w {vW:.3f}")
rows = np.array(rows)
P(f"  MEAN over chains:  share(E[kappa] channel)={rows[:,0].mean()/rows[:,3].mean():.3f}"
  f"  share(Var(F) z channel)={rows[:,1].mean()/rows[:,3].mean():.3f}"
  f"  share(cov)={rows[:,2].mean()/rows[:,3].mean():+.3f}")
P(f"  >>>> paired estimator w=kappa-F^2 cuts the per-sample variance by"
  f" x{rows[:,3].mean()/rows[:,4].mean():.2f}  (n=4 chains, spread "
  f"{rows[:,3].mean()/rows[:,4].min():.2f}..{rows[:,3].mean()/rows[:,4].max():.2f})")

# ---------------------------------------------------------------- 2. correlations (discipline: report rho)
P("\n[2] control-variate correlations (measured, same-sample pairs, 4 chains)")
allk, allf = K.ravel(), F.ravel()
zz = (allf - allf.mean()) ** 2
ff2 = allf ** 2
def rho(a, b):
    return float(np.corrcoef(a, b)[0, 1])
cands = {
    "F            (E[F]=0 exact)": (allk, allf),
    "F^2          (E[F^2]=Var F)": (allk, ff2),
    "z=(F-mu)^2   (Var channel)": (allk, zz),
    "kappa        (vs F channel)": (ff2, allk),
}
for nm, (a, b) in cands.items():
    r = rho(a, b)
    P(f"  rho({nm:28s}) = {r:+.4f}   => CR factor 1-rho^2 = {1-r*r:.4f}"
      f"   (kappa side x{(1-r*r)**-1:.3f})")
P(f"  rho(F, |F|) = {rho(allf, np.abs(allf)):+.4f}   rho(F^2,|F|) = {rho(ff2, np.abs(allf)):+.4f}")
P(f"  rho(kappa, |F|) = {rho(allk, np.abs(allf)):+.4f}")

# ---------------------------------------------------------------- 3. stratification (post-hoc, honest split)
P("\n[3] stratification: bins from a PILOT (chains 0,1), variance measured on chains 2,3")
pil = slice(0, 2)
eva = slice(2, 4)
for keyname, key in (("|F|", np.abs(F)), ("kappa", K)):
    kp = key[:, pil].ravel()
    edges = np.quantile(kp, np.linspace(0, 1, 5)[1:-1])
    P(f"  --- stratify by {keyname}, 4 equal-count pilot bins, edges={np.round(edges,3)}")
    idx = np.digitize(key[:, eva].ravel(), edges)
    wE = (K[:, eva].ravel() - F[:, eva].ravel() ** 2)
    v_tot = wE.var()
    ws, ns_ = [], []
    for b in range(4):
        s = wE[idx == b]
        ws.append(s.var(ddof=1)); ns_.append(s.size)
        P(f"      bin{b}: n={s.size:6d}  mean(w)={s.mean():+.4f}  var(w)={s.var(ddof=1):.4f}")
    ws, ns_ = np.array(ws), np.array(ns_)
    v_str = (ns_ / ns_.sum() * ws).sum()          # equal-cost allocation w_i = n_i/N
    P(f"      weighted within-stratum var = {v_str:.4f}  vs total {v_tot:.4f}"
      f"   => VR factor x{v_tot/v_str:.2f}  (pilot overhead x1.5 total cost)")
    # optimal (Neyman-ish) allocation at same total N
    opt = (np.sqrt(ns_ * ws).sum() ** 2) / ns_.sum()
    P(f"      optimal allocation (Neyman, same N) var = {opt:.4f} => x{v_tot/opt:.2f}")

# ---------------------------------------------------------------- 4. strict-bound arithmetic
P("\n[4] what actually sets N: the HOEFFDING bound, not the variance")
DELTA = 0.3942          # o24 point value
EPS = 0.30              # r04 target half-width
DP = 1e-5               # r04-style failure prob (ln(2/dp)/2 = 6.0); see sensitivity below
dstar = 3e-2
Kh = 2.0 / dstar ** 2   # 2222  (kappa <= 2/delta^2 ; O24 heavy-tail mechanism)
P(f"  range bound K = 2/d*^2 = {Kh:.1f}   sigma_1(measured, per sample) = "
  f"{np.sqrt(rows[:,3].mean()):.2f}   ratio K/sigma_1 = {Kh/np.sqrt(rows[:,3].mean()):.1f}")
for dp in (1e-5, 1e-3):
    Nh = Kh ** 2 * np.log(2 / dp) / (2 * EPS ** 2)
    P(f"  (S1) Hoeffding, dp={dp:g}: ln(2/dp)={np.log(2/dp):.2f}  N = {Nh:.3e}   [r04: ~3e8]")
s1 = np.sqrt(rows[:,3].mean())
ln3 = np.log(3 / DP)


def n_bern(sig, eps=EPS, K=Kh, ln=ln3):
    """smallest N with sig*sqrt(2 ln/N) + 3K ln/N <= eps  (empirical Bernstein, Maurer-Pontil)."""
    lo, hi = 1.0, 1e14
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if sig * np.sqrt(2 * ln / mid) + 3 * K * ln / mid <= eps:
            hi = mid
        else:
            lo = mid
    return hi


P(f"  (S2) empirical Bernstein sigma*sqrt(2ln3/N)+3K ln3/N, dp={DP:g} (ln3={ln3:.2f}):")
for lab, sig in (("sigma_1 raw", s1),
                 ("x0.5 var", s1 * 0.5), ("x0.1 var", s1 * 0.1), ("x0.01 var", s1 * 0.01),
                 ("perfect (sigma->0)", 0.0)):
    P(f"        {lab:20s} sigma={sig:8.2f}  -> N = {n_bern(sig):.3e}")
P(f"        floor (sigma=0) = 3K ln3/eps = {3*Kh*ln3/EPS:.3e}   <- the K-term is a 1/N floor")
P(f"  (S2') same with the PAIRED estimator: var reduced x{rows[:,3].mean()/rows[:,4].mean():.2f}"
  f" -> sigma={s1/np.sqrt(rows[:,3].mean()/rows[:,4].mean()):.2f}"
  f" -> N = {n_bern(s1/np.sqrt(rows[:,3].mean()/rows[:,4].mean())):.3e}")

# truncated Hoeffding: bound kappa by T=2/t^2, pay the shell tail P(d<t)<=64 c t^3
P("  (S3) truncated Hoeffding: min(kappa,T) with T=2/t^2 + tail bias K*P(d<t):")
c_meas = float(np.mean(K > Kh))
P(f"        measured P(kappa > K) = {c_meas:.2e} (n={K.size}) => the K=2222 truncation"
  f" is {'SAFE' if c_meas < 0.5/K.size*3 else 'NOT confirmed at this n'}")
for t, cc in ((3e-2, 1.0), (5e-2, 1.0), (1e-1, 1.0), (2e-1, 1.0)):
    T = 2 / t ** 2
    # tail: P(d<t) <= (t/pi)^3 * 64 * cc  (shell volume, slope checked 3.0 in r04 §4.1)
    bnd = 2 * np.pi * 64 * 1.0 * (t / np.pi) ** 3
    print_end = ""
    for need in (0.05, 0.10):
        bib = Kh * bnd if bnd < 1 else np.inf
        if bib >= need:
            P(f"        t={t:.2f} T={T:7.1f} P_tail<={bnd:.2e} bias<={bib:.3f} >= {need}"
              f" -> truncation cannot buy eps={need}")
            break
        e2 = need - bib
        N3 = T ** 2 * np.log(2 / DP) / (2 * e2 ** 2)
        P(f"        t={t:.2f} T={T:7.1f} P_tail<={bnd:.2e} bias={bib:.3f}"
          f" | bias+MC<= {need}: N={N3:.3e}")

# ---------------------------------------------------------------- 5. autocorrelation is the other axis
P("\n[5] the tau axis: information per wall-second (data + quoted logs)")
P(f"  o24 metr.: tau_kappa~{np.mean([ips_tau(K[:,g]) for g in range(G)]):.1f},"
  f" 1+2tau_F~{np.mean([ips_factor(F[:,g]) for g in range(G)]):.1f}, cost 10.02 ms/sweep(4 chains)")
tauK_o24 = np.mean([ips_factor(K[:, g]) for g in range(G)])
tauF_o24 = np.mean([ips_factor(F[:, g]) for g in range(G)])
tauK_p02, tauF_p02, c_p02 = 8.6 * 2 - 1, 1.9 * 2 - 1, 16.7     # p02_dyndyn.log lines 19-23
for nm, (tk, tf, cost) in {"o24 (Metropolis)": (tauK_o24, tauF_o24, 10.02),
                           "p02 (HMC xi=0)  ": (1 + 2 * 8.6, 1 + 2 * 1.9, 16.7)}.items():
    P(f"  {nm}: (1+2tau_kappa)={tk:5.1f} (1+2tau_F)={tf:5.1f} cost={cost:5.2f} ms"
      f" | kappa-cost/sweep={tk*cost/1000:7.3f} s  F-cost/sweep={tf*cost/1000:6.3f} s")
P(f"  HMC/esta: kappa x{(tauK_o24*10.02)/(1+2*8.6)/16.7/1:.2f} cheaper, "
  f"F x{(tauF_o24*10.02)/(1+2*1.9)/16.7:.2f} cheaper per unit info")
P(f"  regeneration break-even: need p_reg > 1/(1+2tau) = {1/tauK_o24:.3f} (kappa),"
  f" {1/tauF_o24:.3f} (F) for a Nummelin atom to pay for itself")

# ---------------------------------------------------------------- 6. verdict arithmetic
P("\n[6] verdict arithmetic (eps=0.30, dp=1e-5)")
Nh = Kh ** 2 * np.log(2 / DP) / (2 * EPS ** 2)
vr = rows[:, 3].mean() / rows[:, 4].mean()
nb = n_bern(s1 / np.sqrt(vr))
P(f"  (S1) pure Hoeffding (what r04 counted)          N = {Nh:.2e}")
P(f"  (S2) empirical Bernstein, raw sigma             N = {n_bern(s1):.2e}")
P(f"  (S2b) empirical Bernstein x paired CV (x{vr:.2f})    N = {nb:.2e}")
P(f"  (S2c) empirical Bernstein, sigma->0 (vector floor) N = {3*Kh*ln3/EPS:.2e}")
P(f"  ratio S1/S2c = {Nh/(3*Kh*ln3/EPS):.1f}   ratio S1/S2b = {Nh/nb:.1f}")
P("=" * 104)
