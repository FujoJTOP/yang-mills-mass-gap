"""QI-BRACKET witness — the residue-level-gap bridge as a certificate with measurable hypotheses.

Bridge (WA-ALIGNMENT eq. 5, T2 §2.4 / GAP-G2-edge L2):
    w1 >= theta := (s2 - Ebar)/(s2 - s1) = 1 - (Ebar - s1)/(s2 - s1) ,
with s1 = inf supp nu_w, s2 = SECOND SUPPORT POINT of nu_w (T3 convention: nu = sum w_n delta_{mu_n},
mu = generator eigenvalue; in 0+1 the two-point function is already a Laplace transform, s = E).

New exact object checked here (QI-1), the DEFICIT IDENTITY
    theta = w1 - Delta/(s2 - s1),    Delta := int_{[s2,inf)} (s - s2) d nu_w  >= 0 ,
i.e. the bridge is the two-level approximation and its slack is exactly the scaled excess first
moment above s2.  Consequence QI-2: a uniform N-state low band gives theta = (3-N)/2.

Sections
  A  deficit identity + identifiability of the degenerate (0/0) point
  B  free field: exact Gaussian two-particle cut, L scan, two-point-ray check
  C  deep double well (1D, exact diagonalisation): theta, certificate interval, thresholds
  D  GEVP / Ritz estimate of s2, t-scan, and the TP2 proxy s2 ~ rho^2 s1
  E  4D beta=2.5 read-only reuse (m1_4d_run_results.json): bare Ebar, fake / impossible certificates
  F  YM-like synthetic measure (atom + cut): tail budget Delta <= (w1-w0)(s2-s1)
  G  beta=0 Haar: degenerate level, N-state band criterion
Deterministic, single process, numpy+scipy only. Writes qi_bracket.log + qi_bracket_results.json.
"""
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh

LOG = []
def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    LOG.append(line)

R = {}
rng = np.random.default_rng(20261124)
JW = "research/v3/species/w7/qi_bracket_results.json"


def theta(s1, s2, E):
    """bridge value (s2-E)/(s2-s1)"""
    return (s2 - E) / (s2 - s1)


def deficit(s, w, s2):
    """Delta = int_{[s2,inf)} (s-s2) d nu  (atoms at s2 contribute 0)"""
    return float(np.sum(np.where(s >= s2, w * (s - s2), 0.0)))


def profile(s, w):
    """(s1, s2, Ebar, w1, theta, Delta) + the two-sided bracket, for a normalised sorted measure"""
    s1, s2 = float(s[0]), float(s[1])
    E = float(w @ s)
    D = deficit(s, w, s2)
    lo = max(0.0, theta(s1, s2, E), 1.0 - E / s2)       # 1-w1 <= E/s2 (Markov on leakage)
    hi = min(1.0, E / s1) if s1 > 0 else 1.0            # w1 <= E/s1 (Markov, needs s1>0)
    return dict(s1=s1, s2=s2, eps=s2 - s1, Ebar=E, w1=float(w[0]),
                theta=theta(s1, s2, E), Delta=D, deficit=D / (s2 - s1), lo=lo, hi=hi,
                neff=2.0 * (E - s1) / (s2 - s1) + 1.0)


# ---------------------------------------------------------------- §A identity
say("=" * 78)
say("(A) QI-1 deficit identity  theta = w1 - Delta/(s2-s1)   [random measures]")
worst = 0.0
worst_ineq = 0.0
for _ in range(500):
    K = int(rng.integers(3, 11))
    s = np.cumsum(rng.uniform(0.2, 1.6, K))
    w = rng.random(K)
    w /= w.sum()
    p = profile(s, w)
    worst = max(worst, abs(p["theta"] - (p["w1"] - p["deficit"])))
    worst_ineq = max(worst_ineq, p["theta"] - p["w1"])     # must be <= 0
say(f"    max |theta - (w1 - Delta/eps)| = {worst:.3e}   over 500 random measures")
say(f"    max (theta - w1)               = {worst_ineq:.3e}   (bridge is sound: <= 0)")
say("    identifiability at the degenerate point (s2 -> s1): two-point family nu_eps = w1 d_{s1} +")
say("    (1-w1) d_{s1+eps} has theta = w1 EXACTLY for every eps, while (s1, s2, Ebar) -> (0,0,0)")
say("    for every w1.  So at eps = 0 the triple carries NO information about w1: the 0/0 is real.")
R["A"] = dict(identity_err=worst, max_overclaim=worst_ineq)


# ---------------------------------------------------------------- §B free field
def free_measure(L, m):
    """free scalar on a ring of L sites, source o = phi_0^2 - <phi^2>.
    C(t) = 2 G(t)^2 with G(t) = (1/L) sum_k e^{-w_k t}/(2 w_k)
    => nu = sum_{k,k'} (2/(4 L^2 w_k w_k')) delta_{w_k+w_k'}.   Exact (Gaussian Wick)."""
    k = np.arange(L)
    om = np.sqrt(m * m + 4.0 * np.sin(np.pi * k / L) ** 2)
    s = (om[:, None] + om[None, :]).ravel()
    w = (2.0 / (4.0 * L * L * om[:, None] * om[None, :])).ravel()
    o = np.argsort(s)
    s, w = s[o], w[o]
    M0 = w.sum()
    w = w / M0
    su, wu = [], []
    for si, wi in zip(s, w):
        if su and abs(si - su[-1]) <= 1e-12 * max(1.0, abs(si)):
            wu[-1] += wi
        else:
            su.append(si)
            wu.append(wi)
    return np.array(su), np.array(wu), float(M0)


say("=" * 78)
say("(B) free field (exact Gaussian, source = bilinear): the 0/0 degeneration, eps -> 0")
say(f"    {'L':>5} {'s1':>9} {'s2':>9} {'eps=s2-s1':>10} {'Ebar-s1':>9} {'theta':>11} {'w1':>10} {'Delta':>9} {'bracket':>16}")
R["B"] = []
for L in (8, 16, 32, 64, 128, 256):
    s, w, M0 = free_measure(L, 0.4)
    p = profile(s, w)
    R["B"].append(dict(L=L, **p, M0=M0))
    say(f"    {L:5d} {p['s1']:9.5f} {p['s2']:9.5f} {p['eps']:10.3e} {p['Ebar']-p['s1']:9.5f}"
        f" {p['theta']:11.3e} {p['w1']:10.3e} {p['Delta']:9.3e}  [{p['lo']:.3f},{p['hi']:.3f}]")
say("    eps ~ 1/L^2, Ebar-s1 = O(1) fixed  =>  theta = 1-(Ebar-s1)/eps -> -infinity; w1 -> 0.")
say("    m-scan at L=64 (massless direction; s1 = 2m is the two-particle threshold):")
R["B_m"] = []
for m in (0.4, 0.1, 0.02):
    s, w, M0 = free_measure(64, m)
    p = profile(s, w)
    R["B_m"].append(dict(m=m, **p, M0=M0))
    say(f"      m={m:5.3f}  s1={p['s1']:.5f}  eps={p['eps']:.3e}  theta={p['theta']:.3e}  w1={p['w1']:.3e}")
say("    two-point ray check (the removable direction): eps=1e-3, s1=0, w1 = 0.01/0.5/0.99")
for w1 in (0.01, 0.5, 0.99):
    eps = 1e-3
    s = np.array([0.0, eps])
    w = np.array([w1, 1 - w1])
    p = profile(s, w)
    say(f"      w1={w1:5.2f}: theta={p['theta']:.6f}  (Ebar={p['Ebar']:.2e}, s2={p['s2']:.2e}, Delta={p['Delta']:.1e})")


# ---------------------------------------------------------------- §C double well
def double_well(g, n=4001, xmax=6.0, k=120):
    """H = -1/2 d^2/dx^2 + (g/2)(x^2-1)^2, source o = x.  FD, tridiagonal, lowest k states."""
    x = np.linspace(-xmax, xmax, n)
    h = x[1] - x[0]
    d = 1.0 / h ** 2 + 0.5 * g * (x * x - 1.0) ** 2
    e = -0.5 / h ** 2 * np.ones(n - 1)
    ev, vec = eigh_tridiagonal(d, e, select="i", select_range=(0, k - 1))
    return x, h, ev, vec


say("=" * 78)
say("(C) deep double well: theta vs w1 (soundness + saturation), s2 = SECOND SUPPORT POINT")
say(f"    {'g':>6} {'s1':>10} {'s2supp':>9} {'Ebar':>9} {'w1':>9} {'theta':>9} {'w1-theta':>9} {'Delta':>9} {'cert.int.':>17}")
R["C"] = []
for g in (0.5, 2.0, 8.0, 16.0, 32.0, 64.0, 100.0):
    x, h, ev, vec = double_well(g)
    psi0 = vec[:, 0]
    Ex = float(np.sum(psi0 * psi0 * x))
    Ex2 = float(np.sum(psi0 * psi0 * x * x))
    Var = Ex2 - Ex * Ex
    s = ev - ev[0]
    mtx = (vec.T @ (x * psi0)).ravel()                  # <n|x|0>
    w = mtx ** 2 / Var
    supp = np.where(w[1:] > 1e-12)[0] + 1               # support indices n>=1
    s_supp = s[supp]
    w_supp = w[supp] / w[supp].sum()
    Ebar = 1.0 / (2.0 * Var)                            # [x,[H,x]] = 1 exactly
    Ebar_spec = float(w_supp @ s_supp)
    p = profile(s_supp, w_supp)
    R["C"].append(dict(g=g, **p, Ebar_spec=Ebar_spec, s2_eig=float(s[2]),
                       varx=Var, nsupp=int(len(supp))))
    say(f"    {g:6.1f} {p['s1']:10.4e} {p['s2']:9.4f} {p['Ebar']:9.6f} {p['w1']:9.6f}"
        f" {p['theta']:9.6f} {p['w1']-p['theta']:9.2e} {p['Delta']:9.2e}  [{p['lo']:.6f},{p['hi']:.6f}]")
say(f"    identity check  max|Ebar - Ebar_spec|/Ebar = "
    f"{max(abs(r['Ebar']-r['Ebar_spec'])/r['Ebar'] for r in R['C']):.2e}"
    "   (double-commutator identity, k=120 truncation)")
say("    threshold form: certified iff s2 >= (Ebar - w0 s1)/(1-w0), i.e. theta >= w0")
for w0 in (0.1, 0.3, 0.5):
    ok = sum(1 for r in R["C"] if r["theta"] >= w0)
    ratio = [(r["Ebar"] - w0 * r["s1"]) / (1 - w0) / r["s2"] for r in R["C"]]
    say(f"    w0={w0:4.2f}: certified {ok}/{len(R['C'])}   required_s2/actual_s2 in"
        f" [{min(ratio):.2e}, {max(ratio):.2e}]  (margin factor {1/max(ratio):.0f}x..{1/min(ratio):.0f}x)")

# ---------------------------------------------------------------- §D GEVP
say("-" * 78)
say("(D) GEVP/Ritz estimate of s2 (basis {x,x^3}) and the TP2 proxy s2 ~ rho^2 s1")
say(f"    {'g':>6} {'s1':>9} {'mu1':>9} {'s2supp':>9} {'mu2':>9} {'theta(mu2)':>10} {'w1':>9} {'overclaim':>10} {'rho2*s1':>9} {'err':>8}")
R["D"] = []
for g in (0.5, 2.0, 8.0, 32.0, 100.0):
    x, h, ev, vec = double_well(g)
    psi0 = vec[:, 0]
    Var = float(np.sum(psi0 ** 2 * x ** 2)) - float(np.sum(psi0 ** 2 * x)) ** 2
    s = ev - ev[0]
    m1 = (vec.T @ (x * psi0)).ravel()
    m2 = (vec.T @ ((x ** 3) * psi0)).ravel()
    w = m1 ** 2 / Var
    supp = np.where(w[1:] > 1e-12)[0] + 1
    s_supp, w_supp = s[supp], w[supp] / w[supp].sum()
    p = profile(s_supp, w_supp)
    best = None
    for t in (0.2, 0.5, 1.0, 2.0):
        Ct = np.zeros((2, 2))
        C0 = np.zeros((2, 2))
        for i, mi in enumerate((m1, m2)):
            for j, mj in enumerate((m1, m2)):
                Ct[i, j] = np.sum(mi * mj * np.exp(-s * t)) - mi[0] * mj[0]
                C0[i, j] = np.sum(mi * mj) - mi[0] * mj[0]
        try:
            lam, _ = eigh(Ct, C0)
        except Exception:
            continue
        if np.any(lam <= 0):
            continue
        mu = np.sort(-np.log(lam) / t)                  # ascending: mu[0] ~ s1, mu[1] ~ s2
        cand = theta(p["s1"], float(mu[1]), p["Ebar"])
        if best is None or cand > best[0]:
            best = (cand, float(mu[0]), float(mu[1]), t)
    th_gev, mu1, mu2, t = best
    proxy = mu2 / mu1 * p["s1"]
    R["D"].append(dict(g=g, **p, mu1=mu1, mu2=mu2, theta_gev=th_gev,
                       overclaim=th_gev - p["w1"], rho2=mu2 / mu1, proxy=proxy,
                       proxy_err=proxy / p["s2"] - 1.0, t=t))
    say(f"    {g:6.1f} {p['s1']:9.4f} {mu1:9.4f} {p['s2']:9.4f} {mu2:9.4f} {th_gev:10.6f}"
        f" {p['w1']:9.6f} {th_gev-p['w1']:+10.2e} {proxy:9.4f} {proxy/p['s2']-1.0:+8.1e}")
say("    reading: mu2 >= s2 (Ritz is an upper bound) => theta(mu2) >= theta(s2): the estimate can only")
say("             push the certificate UP (over-claim).  It can also land BELOW s2 with a bad basis, so")
say("             no GEVP/Ritz/rho^2 estimate controls s2 in the direction the certificate needs (below).")
say("             rho^2*s1 tracks s2 only in the near-harmonic regime; at g=100 it is off by orders.")


# ---------------------------------------------------------------- §E 4D read-only
say("=" * 78)
say("(E) 4D beta=2.5 read-only (m1_4d_run_results.json): bare Ebar and fake / impossible certificates")
with open("research/v3/w7/m1_4d_run_results.json") as f:
    J = json.load(f)
R["E"] = []
for L in (4, 6, 8):
    gb = J[f"L{L}"]["glueball"]
    row = dict(L=L, variants={})
    say(f"    L={L}: fit m1 pinned at {gb['fit_unsmeared']['m1']:.4f}")
    for tag, key in (("unsm", "O_unsmeared"), ("ape2", "O_ape2")):
        C = np.array(gb[key], float)
        C = C / C[0]
        Eb = float(-np.log(C[1]))
        row["variants"][tag] = dict(C=[float(c) for c in C], Ebar_bare=Eb,
                                    theta_if_s2={str(s2): theta(0.0, s2, Eb) for s2 in (0.03, 0.1, 0.3, 1.0)},
                                    s2_required={str(w0): Eb / (1 - w0) for w0 in (0.1, 0.3, 0.5)})
        say(f"         [{tag:4s}] C(t)/C(0)={np.round(C,5).tolist()}")
        say(f"                Ebar = -ln[C(1)/C(0)] = {Eb:+.3e}  bridge(s1=0): " +
            ", ".join(f"s2={k}: {v:.4f}" for k, v in row["variants"][tag]["theta_if_s2"].items()))
    R["E"].append(row)
    d4 = row["variants"]["unsm"]["theta_if_s2"]["0.3"] - row["variants"]["ape2"]["theta_if_s2"]["0.3"]
    say(f"         smoothing spread of theta at s2=0.3 (unsm vs ape2) = {d4:+.2e}"
        "  => O(1) scheme sensitivity of the certificate, not a renormalisation")
say("    L=4: Ebar_bare = +4.8e-3 > 0 (apparently valid): theta >= 0.84 for any s2 >= 0.03, >= 0.98 for")
say("         s2 >= 0.3.  It certifies that >=98% of the source variance sits at s <= 0.005 -- i.e. in the")
say("         contact/zero mode, NOT in an atom.  TRUE but MEANINGLESS (vacuous threshold ~1e-2).")
say("    L=6,8: Ebar_bare < 0  =>  theta > 1 (1.002-1.072): an IMPOSSIBLE value (w1 <= 1).  The bridge")
say("         self-diagnoses: theta>1 iff Ebar<s1, i.e. the hypothesis 's1 is a support point' is violated.")


# ---------------------------------------------------------------- §F synthetic
def synth(m1sq, w1, s2, Lam, alpha=0.0, N=4001):
    s = np.linspace(s2, Lam, N)
    ds = s[1] - s[0]
    tw = np.full(N, ds)
    tw[0] = tw[-1] = 0.5 * ds                            # trapezoid weights
    wt = s ** (-alpha) * tw
    wc = wt / wt.sum() * (1.0 - w1)
    s_all = np.concatenate([[m1sq], s])
    w_all = np.concatenate([[w1], wc])
    return s_all, w_all / w_all.sum()


say("=" * 78)
say("(F) YM-like measure: atom at m1^2=1 (w1=0.30) + two-particle cut from s2=4. Tail budget")
say(f"    {'Lambda':>7} {'alpha':>6} {'Ebar':>8} {'Delta':>8} {'theta':>9} {'>=0.1':>7} {'>=0.3':>7}")
R["F"] = []
for Lam in (4.5, 5.0, 5.5, 6.0, 8.0, 16.0, 64.0):
    for alpha in (0.0, 2.0):
        s, w = synth(1.0, 0.30, 4.0, Lam, alpha)
        p = profile(s, w)
        R["F"].append(dict(Lam=Lam, alpha=alpha, **p))
        say(f"    {Lam:7.1f} {alpha:6.1f} {p['Ebar']:8.4f} {p['Delta']:8.4f} {p['theta']:9.4f}"
            f" {str(p['theta']>=0.1):>7} {str(p['theta']>=0.3):>7}")
say("    budget: theta >= w0  <=>  Delta <= (w1-w0)(s2-s1).  Flat cut: Delta = 0.7(Lambda-s2)/2,")
say("    so w0=0.1 needs Lambda-s2 <= 2(0.2)(3)/0.7 = 1.71 = 0.57 eps ; w0=0.3 needs Lambda=s2 exactly.")
say("    => the bridge needs the source's whole spectral weight above s2 to sit within O(eps) of s2.")


# ---------------------------------------------------------------- §G beta=0 / bands
say("=" * 78)
say("(G) beta=0 Haar: degenerate level => N-state band criterion")
q = rng.normal(size=(400000, 4))
q /= np.linalg.norm(q, axis=1)[:, None]
say(f"    SU(2) Haar: Var((1/2)tr U) = {q[:,0].var():.4f}  (exact 1/4)  => M0 finite at beta=0")
say("    beta=0: T = |1><1| (rank one) => every excitation at s = infinity: nu = delta_infinity,")
say("            no (s1,s2,Ebar) at all => the bridge is UNDEFINED (an inf/inf, not a 0/0).")
say("    regularised: N states spread uniformly over a band of width W => eps=W/(N-1), Ebar-s1=eps(N-1)/2")
say("                 theta = (3-N)/2  (independent of W);  theta>0 needs N<=2")
R["G"] = dict(haar_var=float(q[:, 0].var()), bands=[])
for N in (2, 3, 4, 8, 100, 1000):
    s = np.linspace(1.0, 2.0, N)
    w = np.ones(N) / N
    p = profile(s, w)
    R["G"]["bands"].append(dict(N=N, theta=p["theta"], w1=p["w1"], eps=p["eps"]))
    say(f"      N={N:5d}: eps={p['eps']:.3e}  theta={p['theta']:+.4f}  (w1={p['w1']:.4f})")
say("    => few-level certificate: N states in the low band give theta = (3-N)/2 <= 0 for N>=3.")
say("       free field / beta=0 / any cut: N -> infinity => theta -> -infinity.  No certificate.")


# ---------------------------------------------------------------- §H summary
say("=" * 78)
say("(H) summary")
say("    [T] QI-1  theta = w1 - Delta/(s2-s1),  Delta = int_{[s2,inf)}(s-s2)dnu >= 0   (exact)")
say("    [T] QI-2  uniform N-state band: theta = (3-N)/2  => N<=2 needed for any w0>0")
say("    [T] QI-3  theta>1 (resp. theta<0) iff Ebar<s1 (resp. Ebar>s2): a built-in consistency flag")
say("    [V] free field: eps=s2-s1 ~ 1/L^2 -> 0, w1 -> 0, theta -> -inf  (no positive certificate)")
say("    [V] double well: theta = w1 - (5e-5..1.4e-3), certified at w0=0.1/0.3/0.5 (saturated)")
say("    [V] synthetic YM-like: theta<0 for any cut wider than ~0.6 eps  (tail budget)")
say("    [V] 4D bare read-only: Ebar_bare ~ 0 => theta>=0.98 (meaningless) or theta>1 (impossible)")
say("    [NG] s2 from GEVP/Ritz/rho^2: upper-bound machines, wrong direction; theta over-claims")
with open("research/v3/species/w7/qi_bracket.log", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG) + "\n")
with open(JW, "w", encoding="utf-8") as f:
    json.dump(R, f, indent=1, default=float)
say("    wrote qi_bracket.log + qi_bracket_results.json")
