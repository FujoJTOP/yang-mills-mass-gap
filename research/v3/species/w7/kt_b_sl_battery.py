"""kt_b_sl_battery.py -- KT-RHO2 T1: random Sturm-Liouville battery.

The decisive question: does ANY member of the (physically natural) single-well class
violate rho^2 <= 4 ?   Classes: DD / NN / ND boundary patterns, random potential q>=0,
random density w>0 (string), random stiffness p>0 (mixed), 1D chains up to L=40.
Finite differences, single process, deterministic seeds.
Writes kt_b_battery.json + kt_b_battery.log after every class (crash-safe).
"""
import json, numpy as np
from scipy.linalg import eigh_tridiagonal, eigh

JSON, LOG = "kt_b_battery.json", []
OUT = {}


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


def flush():
    with open(JSON, "w", encoding="utf-8") as f:
        json.dump(OUT, f, indent=1)
    with open("kt_b_battery.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


def sl_sys(N=400, q=None, w=None, p=None, bc="DD"):
    """Symmetric-tridiagonal representation of  -(p u')' + q u = lam w u  on [0,1].

    Dirichlet: u(0)=u(1)=0 (nodes 1..N-1).  Neumann: ghost mirror, 2*p0/h^2 on the
    diagonal and -2*p0/h^2 off it.  ND = Neumann at 0, Dirichlet at 1 (nodes 0..N-1).
    The density w is removed by the symmetric transform u -> sqrt(w) u.
    """
    h = 1.0 / N
    x = np.linspace(0, 1, N + 1)
    pfun = (lambda t: np.ones_like(t)) if p is None else p
    pn = pfun(x)
    ph = 0.5 * (pn[:-1] + pn[1:])            # p at (j+1/2)h, j=0..N-1
    qn = np.zeros(N + 1) if q is None else np.asarray(q(x), float)
    idx = {"DD": np.arange(1, N), "ND": np.arange(0, N), "NN": np.arange(0, N + 1)}[bc]
    m = len(idx)
    diag = np.zeros(m)
    off = np.zeros(max(m - 1, 0))
    for a, j in enumerate(idx):
        if j == 0 and bc in ("NN", "ND"):
            diag[a] = ph[0] / h ** 2 + qn[0]
            off[a] = -ph[0] / h ** 2
            continue
        d = qn[j] + (ph[j - 1] / h ** 2 if j >= 1 else 0.0)
        if j < N:
            d += ph[j] / h ** 2
        elif bc == "NN":                      # mirror at the far end
            d = qn[N] + ph[N - 1] / h ** 2
        diag[a] = d
        if a < m - 1:
            off[a] = -ph[j] / h ** 2
    if w is not None:
        ww = np.asarray(w(x), float)[idx]
        diag = diag / ww
        off = off / np.sqrt(ww[:-1] * ww[1:])
    return diag, off


def sl_rho2(N=400, q=None, w=None, p=None, bc="DD"):
    diag, off = sl_sys(N, q, w, p, bc)
    lam = eigh_tridiagonal(diag, off, eigvals_only=True)
    lam = lam[lam > 1e-9 * max(1.0, lam[-1])]      # drop a possible zero mode
    return lam


def rand_q(rng, kind="smooth", J=6):
    if kind == "smooth":
        a = rng.exponential(1.0, J)
        ph = rng.uniform(0, 2 * np.pi, J)
        f = 0.5 + 1.5 * rng.random()
        return lambda x: sum(a[i] * (1 + np.cos(2 * np.pi * (i + 1) * f * x + ph[i]))
                             for i in range(J)) / J
    if kind == "bumps":
        c = rng.uniform(0.1, 0.9, 3)
        s = 10 ** rng.uniform(-2.5, -1.0, 3)
        hh = rng.exponential(3.0, 3)
        return lambda x: sum(hh[i] * np.exp(-((x - c[i]) / s[i]) ** 2) for i in range(3))
    if kind == "steps":
        k = rng.integers(2, 5)
        cuts = np.sort(rng.uniform(0, 1, k - 1))
        vals = rng.exponential(2.0, k)
        return lambda x: vals[np.searchsorted(cuts, x)]
    raise ValueError(kind)


def _logfield(rng, J, amp):
    a = rng.normal(0, 1.0, J)
    ph = rng.uniform(0, 2 * np.pi, J)
    return lambda x: np.exp(amp * sum(a[i] * np.cos(2 * np.pi * (i + 1) * x + ph[i])
                                     for i in range(J)))


def rand_w(rng, J=5, amp=0.7):
    return _logfield(rng, J, amp)


def rand_p(rng, J=4, amp=0.7):
    return _logfield(rng, J, amp)


P("=" * 100)
P("B0  sanity: the FD operator reproduces the exact BC towers (DD=4, NN=4, ND=9)")
P("=" * 100)
sane = {}
for bc, exact in [("DD", 4.0), ("NN", 4.0), ("ND", 9.0)]:
    lam = sl_rho2(N=1200, bc=bc)
    r = lam[1] / lam[0]
    sane[bc] = dict(fd=float(r), exact=exact)
    P(f"    {bc}: FD rho^2 = {r:.8f}   exact = {exact:.1f}   rel.err = {abs(r-exact)/exact:.2e}")
OUT["B0_sanity"] = sane
flush()

P("=" * 100)
P("B1  random battery: 2000 samples per class, N=400, rho^2 = lam2/lam1")
P("=" * 100)
rng = np.random.default_rng(20261124)
classes = [
    ("DD_q>=0_smooth", "DD", dict(q="smooth")),
    ("DD_q>=0_bumps", "DD", dict(q="bumps")),
    ("DD_q>=0_steps", "DD", dict(q="steps")),
    ("DD_string_w", "DD", dict(w=True)),
    ("DD_mixed_p_q_w", "DD", dict(p=True, q="smooth", w=True)),
    ("NN_q>=0", "NN", dict(q="smooth")),
    ("ND_q>=0", "ND", dict(q="smooth")),
    ("ND_q>=0_bumps", "ND", dict(q="bumps")),
]
NS = 2000
for name, bc, spec in classes:
    rs = []
    for _ in range(NS):
        kw = {}
        if "q" in spec:
            kw["q"] = rand_q(rng, spec["q"])
        if spec.get("w"):
            kw["w"] = rand_w(rng)
        if spec.get("p"):
            kw["p"] = rand_p(rng)
        lam = sl_rho2(N=400, bc=bc, **kw)
        rs.append(lam[1] / lam[0])
    rs = np.array(rs)
    d = dict(bc=bc, n=NS, sup=float(rs.max()), median=float(np.median(rs)),
             lo=float(rs.min()), n_viol=int((rs > 4.0).sum()))
    OUT.setdefault("B1_battery", {})[name] = d
    P(f"  {name:<18s} bc={bc}  n={NS}  sup={rs.max():.8f}  median={np.median(rs):.5f}"
      f"  min={rs.min():.5f}   #(rho2>4)={d['n_viol']}")
    flush()

P("=" * 100)
P("B2  deterministic families (DD: the ratio decreases away from q=0; ND: sup 9 at q=0)")
P("=" * 100)
det = {}
for c in [0.0, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4]:
    lam = sl_rho2(N=800, bc="DD", q=(lambda x, c=c: c * (x - 0.5) ** 2))
    det[f"DD_parab_c{c:g}"] = float(lam[1] / lam[0])
    P(f"    DD q={c:>8g}(x-1/2)^2        rho^2 = {lam[1]/lam[0]:.8f}")
for c in [0.0, 1.0, 1e2, 1e4]:
    lam = sl_rho2(N=800, bc="DD", q=(lambda x, c=c: c * ((x < 0.1) + 1.0 * (x > 0.9))))
    det[f"DD_walls_c{c:g}"] = float(lam[1] / lam[0])
    P(f"    DD q={c:>8g} on the two edges  rho^2 = {lam[1]/lam[0]:.8f}")
for c in [0.0, 1.0, 10.0, 1e2]:
    lam = sl_rho2(N=800, bc="ND", q=(lambda x, c=c: c * (x - 0.5) ** 2))
    det[f"ND_parab_c{c:g}"] = float(lam[1] / lam[0])
    P(f"    ND q={c:>6g}(x-1/2)^2       rho^2 = {lam[1]/lam[0]:.8f}")
# two-step string, cross-checked against the exact secular equation and the dense problem
from scipy.optimize import brentq


def exact_two_step(a, R, k=3):
    """roots of cot(w a) + sqrt(R) cot(sqrt(R) w (1-a)) = 0  (pole crossings filtered)"""
    fr = lambda t: 1.0 / np.tan(t * a) + np.sqrt(R) / np.tan(np.sqrt(R) * t * (1 - a))
    grid = np.linspace(1e-6, 20.0, 400001)
    v = 1.0 / np.tan(grid * a) + np.sqrt(R) / np.tan(np.sqrt(R) * grid * (1 - a))
    out = []
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0 and abs(v[i]) < 1e3 and abs(v[i + 1]) < 1e3:
            out.append(brentq(fr, grid[i], grid[i + 1]))
    return np.array(out[:k])


for rat in [1.0, 10.0, 100.0]:
    wfun = (lambda x, rat=rat: np.where(x < 0.5, 1.0, rat))
    lam = sl_rho2(N=1600, bc="DD", w=wfun)
    ws = exact_two_step(0.5, rat, 3)
    ex = (ws[1] / ws[0]) ** 2 if len(ws) > 1 else np.nan
    det[f"DD_string_step_{rat:g}"] = float(lam[1] / lam[0])
    det[f"DD_string_step_{rat:g}_exact"] = float(ex)
    P(f"    DD string w=1|{rat:g} step   rho^2 = {lam[1]/lam[0]:.8f}"
      f"   (exact secular: {ex:.8f};  omega1={ws[0]:.6f})")
OUT["B2_det"] = det
flush()

P("=" * 100)
P("B2b  adversarial DD families: is `q>=0` enough, or is the real hypothesis `no interior dip`?")
P("=" * 100)
adv = {}
for H in [1e2, 1e4, 1e6, 1e8]:
    for s in [0.5, 0.2, 0.05]:
        qf = (lambda x, H=H, s=s: H * (1.0 * (np.abs(x - 0.5) > s / 2)))
        lam = sl_rho2(N=1200, bc="DD", q=qf)
        adv[f"one_dip_H{H:g}_s{s:g}"] = float(lam[1] / lam[0])
        P(f"    one dip  H={H:<8g} s={s:<5g}   rho^2 = {lam[1]/lam[0]:.8f}")
for H in [1e2, 1e4, 1e6]:
    qf = (lambda x, H=H: H * (1 - np.exp(-((x - 0.3) / 0.02) ** 2)
                              - np.exp(-((x - 0.7) / 0.02) ** 2)))
    lam = sl_rho2(N=1600, bc="DD", q=qf)
    adv[f"two_dips_H{H:g}"] = float(lam[1] / lam[0])
    P(f"    two dips H={H:<8g} (0.3,0.7) rho^2 = {lam[1]/lam[0]:.8f}   "
      f"(merged; raw pair split = {(lam[1]-lam[0])/lam[0]:.2e})")
for n in [1, 2, 3, 5]:
    qf = (lambda x, n=n: 200.0 * (1 - np.cos(2 * np.pi * n * x)) / 2)
    lam = sl_rho2(N=1600, bc="DD", q=qf)
    adv[f"n_dips_{n}"] = float(lam[1] / lam[0])
    P(f"    {n} dips (cos comb)        rho^2 = {lam[1]/lam[0]:.8f}")
OUT["B2b_adversarial_DD"] = adv
flush()

P("=" * 100)
P("B2c  how far can a (physical) non-uniform DENSITY push rho^2 ?  DD, N=600")
P("=" * 100)
dens = {}
for amp in [0.5, 1.0, 2.0, 3.0]:
    rs = []
    for _ in range(400):
        wf = rand_w(rng, J=6, amp=amp)
        lam = sl_rho2(N=600, bc="DD", w=wf)
        rs.append(lam[1] / lam[0])
    rs = np.array(rs)
    dens[f"rand_w_amp{amp:g}"] = dict(sup=float(rs.max()), median=float(np.median(rs)),
                                      viol=int((rs > 4.0).sum()))
    P(f"    random log-density amp={amp:<4g} (contrast <= e^{2*amp:g}={np.exp(2*amp):.0f})"
      f"  sup={rs.max():8.3f}  median={np.median(rs):6.3f}  #(>4)={int((rs>4.0).sum())}/400")
for C in [10.0, 1e2, 1e4]:
    for s in [0.02, 0.1]:
        wf = (lambda x, C=C, s=s: 1.0 + C * np.exp(-((x - 0.5) / s) ** 2))
        lam = sl_rho2(N=1500, bc="DD", w=wf)
        v = lam[1] / lam[0]
        dens[f"lump_C{C:g}_s{s:g}"] = float(v)
        dens[f"lump_C{C:g}_s{s:g}"] = float(v)
        P(f"    heavy lump C={C:<7g} s={s:<5g}            rho^2 = {v:8.4f}")
OUT["B2c_density"] = dens
flush()

P("=" * 100)
P("B3  1D chains L=2..40: rho^2 -> 4 from below (never above)")
P("=" * 100)
ch = {}
for L in [2, 3, 4, 6, 8, 12, 20, 40]:
    d = np.full(L, 2.0)
    e = np.full(L - 1, -1.0)
    mu = eigh_tridiagonal(d, e, eigvals_only=True)      # no vacuum: this matrix is Dirichlet-ended
    if len(mu) < 3:
        ch[f"L{L}"] = dict(rho2=None, note="too few levels")
        P(f"    path L={L:<3d} only {len(mu)} level(s)")
        continue
    ch[f"L{L}"] = dict(rho2=float(mu[1] / mu[0]), r3=float(mu[2] / mu[1]))
    P(f"    path L={L:<3d} rho^2={mu[1]/mu[0]:.10f}  r3={mu[2]/mu[1]:.8f}")
OUT["B3_chains"] = ch
flush()

P("=" * 100)
P("B4  verdict tally")
P("=" * 100)
sup = lambda pre: max(v["sup"] for k, v in OUT["B1_battery"].items() if k.startswith(pre))
viol = lambda pre: int(sum(v["n_viol"] for k, v in OUT["B1_battery"].items() if k.startswith(pre)))
P(f"    sup rho^2 :  DD = {sup('DD'):.8f}   NN = {sup('NN'):.8f}   ND = {sup('ND'):.8f}")
P(f"    violations (rho^2>4):  DD = {viol('DD')} / {5*NS}   NN = {viol('NN')} / {NS}"
  f"   ND = {viol('ND')} / {2*NS}")
OUT["B4_verdict"] = dict(sup_DD=float(sup("DD")), sup_NN=float(sup("NN")), sup_ND=float(sup("ND")),
                         viol_DD=viol("DD"), viol_NN=viol("NN"), viol_ND=viol("ND"))
flush()
P("done. json ->", JSON)
