#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p02_dyndyn.py -- push  Delta_dyn = E[H_vv] - Var(G_v)  to >=5 sigma on L=2, SU(2), beta=2.5.

Model = o27_adjudicate.py Tier-2 model (= o13_ldep.py HMC): constrained (xi=0) HMC,
flat chart X_l = log U_l in R^3, V = S_W + Phi, Phi = -2 log(sin r / r), r in (0,pi).
PLUS the O24-L2-PUSH fundamental-domain correction: |x_l| < pi is the fundamental
domain of SU(2) in THIS chart (U = exp(i x.sigma); shpere |x|=pi -> -I), so any
proposal leaving it is rejected (exact: symmetric proposal + rejection inside the region).

Vectorised over NC chains x 96 plaquettes in ONE numpy process (= ONE job, memory ~50 MB).
Quantities recorded per sweep, all from the SAME configurations:
   Gv = <v, dS_W>          Phi' = <v, dPhi>      Phi'' = v^T Hess_Phi v
   H_vv = v^T Hess_{S_W} v (central FD along v, h=0.05; same discipline as o13/o24)
   Delta_dyn = E[H_vv] - Var(G_v)

Errors: (1) block-jackknife sigma_jk (40 blocks/chain, 4 chains merged);
        (2) autocorrelation-corrected sigma_tau (Geyer initial positive sequence),
            Var(E[H]) = Var(H)(1+2 tau_H)/N ,
            Var(Var G) = Var((G-mu)^2)(1+2 tau_z)/N , z=(G-mu)^2,
            Cov(E[H],Var G) from the symmetrised cross-covariance sum.

[V] self-checks: vectorised S_W/force vs o13_ldep.sw_and_force on random configs,
    analytic force vs central FD, E[G+Phi'] = 0, algebraic identity residual = 0.

Run:  py -3.12 p02_dyndyn.py [ntraj] [nchain]
Labels [E].  No git commit.
"""
import sys, time, json
import numpy as np

np.seterr(all='ignore')

BETA = 2.5
L = 2
NTRAJ = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 4
NSTEP = 8
EPS0 = 0.1402
HFD = 0.05
THERM = 500
RMAX = np.pi
P = lambda *a: print(*a, flush=True)
T0 = time.time()

# ------------------------------------------------------------------ quaternions
def qmul(p, q):
    p0 = p[..., 0:1]; pv = p[..., 1:]
    q0 = q[..., 0:1]; qv = q[..., 1:]
    return np.concatenate([p0 * q0 - (pv * qv).sum(-1, keepdims=True),
                           p0 * qv + q0 * pv - np.cross(pv, qv)], -1)

def qconj(p):
    o = p.copy(); o[..., 1:] *= -1.0; return o

def sincos_series(r, sm, rs):
    f = np.where(sm, 1.0 - r * r / 6.0, np.sin(rs) / rs)
    fp = np.where(sm, -r / 3.0, (rs * np.cos(rs) - np.sin(rs)) / rs ** 2)
    sr = np.where(sm, r, np.sin(rs))
    cr = np.where(sm, 1.0 - 0.5 * r * r, np.cos(rs))
    return f, fp, sr, cr

def x_to_q(x):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    q0 = np.where(sm, 1.0 - 0.5 * r * r, np.cos(rs))
    f = np.where(sm, 1.0 - r * r / 6.0, np.sin(rs) / rs)
    return np.concatenate([q0, f * x], -1)

# ------------------------------------------------------------------ geometry (plaquette lists)
def geom():
    n = L ** 4
    idx = np.arange(n)
    c = np.stack([idx % L, (idx // L) % L, (idx // L ** 2) % L, (idx // L ** 3) % L], 1)
    def shift(mu, s):
        cc = c.copy(); cc[:, mu] = (cc[:, mu] + s) % L
        return cc[:, 0] + L * cc[:, 1] + L ** 2 * cc[:, 2] + L ** 3 * cc[:, 3]
    sp = np.stack([shift(m, +1) for m in range(4)])      # (4,n)
    mu_l, nu_l = [], []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            mu_l.append(np.full(n, mu)); nu_l.append(np.full(n, nu))
    mu_p = np.concatenate(mu_l); nu_p = np.concatenate(nu_l)
    s_p = np.tile(idx, 6)
    return dict(n=n, idx=idx, sp=sp, mu=mu_p, nu=nu_p, s=s_p,
                c2s=sp[mu_p, s_p], c3s=sp[nu_p, s_p])

G = geom()
N = G['n']

# ------------------------------------------------------------------ S_W, force, Phi
def force_corner(xx, m, sg):
    r = np.linalg.norm(xx, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, fp, sr, _ = sincos_series(r, sm, rs)
    mv = m[..., 1:]; m0 = m[..., 0:1]
    xd = (xx * mv).sum(-1, keepdims=True)
    return BETA * (sg * f * mv + (xx / rs) * (m0 * sr + sg * fp * xd))

def sw_and_force(x, want_force=True):
    """x: (C,4,n,3) -> (S_W (C,), dS_W/dx (C,4,n,3))"""
    C = x.shape[0]
    q = x_to_q(x)
    mu, nu, s, c2, c3 = G['mu'], G['nu'], G['s'], G['c2s'], G['c3s']
    A = q[:, mu, s]; B = q[:, nu, c2]; Cc = qconj(q[:, mu, c3]); D = qconj(q[:, nu, s])
    ABC = qmul(qmul(A, B), Cc)
    Pq = qmul(ABC, D)
    tot = BETA * (mu.size - Pq[:, :, 0].sum(1))
    if not want_force:
        return tot, None
    m1 = qmul(qmul(B, Cc), D); m2 = qmul(qmul(Cc, D), A); m3 = qmul(qmul(D, A), B)
    Ff = np.zeros((C, 4 * N, 3))
    rows = np.repeat(np.arange(C)[:, None], mu.size, 1)
    for link, site, sg, mm in ((mu, s, +1, m1), (nu, c2, +1, m2),
                               (mu, c3, -1, m3), (nu, s, -1, ABC)):
        xs = x[:, link, site]
        vals = force_corner(xs, mm, sg)
        np.add.at(Ff, (rows, link * N + site), vals)
    return tot, Ff.reshape(C, 4, N, 3)

def phi_of(x):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, _, _, _ = sincos_series(r, sm, rs)
    return -2.0 * np.log(f)

def grad_phi(x):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, fp, _, _ = sincos_series(r, sm, rs)
    ur = np.where(sm, -1.0 / 3.0, fp / (f * rs))
    return -2.0 * ur * x

def hess_phi_vv(x, w):
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    sm = r < 1e-7
    rs = np.where(sm, 1.0, r)
    f, fp, sr, cr = sincos_series(r, sm, rs)
    fpp = np.where(sm, -1.0 / 3.0 + r * r / 30.0,
                   ((2.0 - rs ** 2) * sr - 2.0 * rs * cr) / rs ** 3)
    u = np.where(sm, 0.0, fp / f)
    ur = np.where(sm, -1.0 / 3.0, fp / (f * rs))
    vv = fpp / f - u * u
    xh = x / rs
    wpar = (w * xh).sum(-1, keepdims=True)
    w2 = (w * w).sum(-1, keepdims=True)
    return -2.0 * (vv * wpar ** 2 + ur * (w2 - wpar ** 2))

# ------------------------------------------------------------------ mode v (k=0, mu=0, sigma_1), ||v||=1
V = np.zeros((4, N, 3)); V[0][:, 0] = 1.0 / L ** 2
V4 = V[None]

def dotv(a):
    return (V4 * a).sum(axis=(1, 2, 3))

# ------------------------------------------------------------------ measurement
def measure(x):
    """returns (C,5): Gv, Phi', Phi'', H_vv, S_W"""
    S0, Fx = sw_and_force(x, True)
    Gv = dotv(Fx)
    Ph = dotv(grad_phi(x))
    Pip = hess_phi_vv(x, V4).sum(axis=(1, 2, 3))
    Sp, _ = sw_and_force(x + HFD * V4, False)
    Sm, _ = sw_and_force(x - HFD * V4, False)
    Hv = (Sp - 2.0 * S0 + Sm) / HFD ** 2
    return np.stack([Gv, Ph, Pip, Hv, S0], 1)

# ------------------------------------------------------------------ [V] checks
def checks():
    import importlib.util
    spec = importlib.util.spec_from_file_location("o13", r"D:\Dev\math\research\v3\w7\o13_ldep.py")
    o13 = importlib.util.module_from_spec(spec); spec.loader.exec_module(o13)
    rng = np.random.default_rng(5)
    g13 = o13.build_geom(L)
    eS = eF = 0.0
    for t in range(3):
        x1 = rng.standard_normal((4, N, 3)) * 0.8
        S1, F1 = o13.sw_and_force(x1, g13, BETA, True)
        S2, F2 = sw_and_force(x1[None], True)
        eS = max(eS, abs(S1 - S2[0]) / abs(S1))
        eF = max(eF, np.abs(F1 - F2[0]).max() / max(np.abs(F1).max(), 1e-30))
    # analytic force vs central FD (vectorised path)
    x = rng.standard_normal((1, 4, N, 3)) * 0.6
    e = rng.standard_normal((1, 4, N, 3)); e /= np.linalg.norm(e)
    _, F = sw_and_force(x, True)
    hh = 1e-6
    Sp, _ = sw_and_force(x + hh * e, False); Sm, _ = sw_and_force(x - hh * e, False)
    eFD = float((abs((Sp - Sm) / (2 * hh) - (F * e).sum()) / abs((F * e).sum())).item())
    P(f"  [V] vectorised vs o13_ldep: rel dev S_W={eS:.2e}  force={eF:.2e}")
    P(f"  [V] analytic force vs central FD (h=1e-6): rel dev={eFD:.2e}")
    return max(eS, eF, eFD)

# ------------------------------------------------------------------ HMC (NC chains in parallel)
def hmc(ntraj, nstep, eps, seed, therm):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((NC, 4, N, 3)) * 0.30
    x -= V4 * dotv(x)[:, None, None, None]
    Vf = lambda z: sw_and_force(z, False)[0] + phi_of(z).sum(axis=(1, 2, 3))
    def force(z):
        f = sw_and_force(z, True)[1] + grad_phi(z)
        return f - V4 * dotv(f)[:, None, None, None]
    S, Fx = sw_and_force(x, True)
    Vx = S + phi_of(x).sum(axis=(1, 2, 3))
    fx = Fx + grad_phi(x); fx -= V4 * dotv(fx)[:, None, None, None]
    acc = np.zeros(NC); natt = 0
    rec = []; trec = []
    t0 = time.time()
    for it in range(ntraj):
        p = rng.standard_normal((NC, 4, N, 3))
        p -= V4 * dotv(p)[:, None, None, None]
        x0 = x.copy(); V0 = Vx; K0 = 0.5 * (p * p).sum(axis=(1, 2, 3))
        p = p + 0.5 * eps * (-fx)
        x = x0 + eps * p
        for k in range(nstep - 1):
            f = force(x)
            p = p - eps * f
            x = x + eps * p
        f = force(x)
        p = p - 0.5 * eps * f
        Sn, Fn = sw_and_force(x, True)
        V1 = Sn + phi_of(x).sum(axis=(1, 2, 3))
        K1 = 0.5 * (p * p).sum(axis=(1, 2, 3))
        dH = (V1 + K1) - (V0 + K0)
        natt += 1
        ok = np.log(rng.random(NC)) < -dH
        # fundamental domain: |x_l| < pi
        rmax = np.linalg.norm(x, axis=-1).max(axis=(1, 2))
        ok &= rmax < RMAX
        idx_bad = ~ok
        if idx_bad.any():
            x[idx_bad] = x0[idx_bad]; V1[idx_bad] = V0[idx_bad]; f[idx_bad] = fx[idx_bad]
            Sn[idx_bad] = S[idx_bad]
            Fn = sw_and_force(x, True)[1]
        acc += ok
        Vx = V1; fx = (Fn + grad_phi(x)); fx -= V4 * dotv(fx)[:, None, None, None]
        S = Sn
        if it >= therm:
            rec.append(measure(x)); trec.append(rmax)
    rec = np.array(rec)                       # (n, NC, 5)
    return rec, acc / max(natt, 1), (np.array(trec) >= RMAX).mean(), time.time() - t0

# ------------------------------------------------------------------ statistics
def tau_geyer(v, ml):
    v = np.asarray(v, float); v = v - v.mean()
    n = v.size; den = (v * v).sum()
    if den <= 0:
        return 0.0, 1.0
    ml = min(ml, n - 2)
    ac = np.array([1.0] + [(v[:-j] * v[j:]).sum() / den for j in range(1, ml + 1)])
    s = 0.0
    for k in range(0, (ml - 1) // 2 + 1):
        gam = ac[2 * k] + ac[2 * k + 1] if 2 * k + 1 <= ml else ac[2 * k]
        if gam <= 0:
            break
        s += gam
    tau = -1.0 + 2.0 * s
    return max(tau, 0.0), 1.0 + 2.0 * max(tau, 0.0)

def xcov_sum(a, b, ml):
    """sum_{j} (gamma_ab(j)+gamma_ba(j)) + gamma_ab(0), Geyer-truncated"""
    a = a - a.mean(); b = b - b.mean(); n = a.size
    vaa = (a * a).sum(); vbb = (b * b).sum()
    if vaa <= 0 or vbb <= 0:
        return 0.0
    ml = min(ml, n // 3)
    g = []
    for j in range(0, ml + 1):
        c = (a[:n - j] * b[j:]).sum() / np.sqrt(vaa * vbb)
        d = (b[:n - j] * a[j:]).sum() / np.sqrt(vaa * vbb)
        g.append(c + d if j > 0 else c)
    g = np.array(g)
    s = 0.0
    for k in range(0, (ml) // 2 + 1):
        gam = g[2 * k] + (g[2 * k + 1] if 2 * k + 1 < g.size else 0.0)
        if gam <= 0:
            break
        s += gam
    return float(s * np.sqrt(vaa * vbb) / n)

def stats_all(rec):
    """rec: (n, NC, 5) -> per-chain & combined estimates with sigma_jk, sigma_tau"""
    n, NCc, _ = rec.shape
    ml = min(600, n // 4)
    out = {}
    names = ['E[H_vv]', "E[Phi'']", 'Var(G_v)', "Var(G_v+Phi')", "Cov(G,Phi')"]
    per = np.zeros((NCc, 12))
    taus = np.zeros((NCc, 6))
    for c in range(NCc):
        Gv, Ph, Pip, Hv = rec[:, c, 0], rec[:, c, 1], rec[:, c, 2], rec[:, c, 3]
        Fv = Gv + Ph
        mH, mPi = Hv.mean(), Pip.mean()
        vG, vF, vP = Gv.var(ddof=1), Fv.var(ddof=1), Ph.var(ddof=1)
        cov = (Gv * Ph).mean() - Gv.mean() * Ph.mean()
        A = mH - vG
        z = (Gv - Gv.mean()) ** 2
        tauH, iH = tau_geyer(Hv, ml)
        tauG, iG = tau_geyer(Gv, ml)
        tauZ, iZ = tau_geyer(z, ml)
        tauF, iF = tau_geyer(Fv, ml)
        covHZ = xcov_sum(Hv, z, ml)
        taus[c] = [tauH, tauG, tauZ, tauF, iH, covHZ]
        per[c] = [mH, mPi, vG, vF, cov, A, iH, iG, iZ, iF, covHZ, (Gv.mean() + Ph.mean())]
    # jackknife: 40 blocks per chain
    nb = 40; k = n // nb
    jk = np.zeros((nb, NCc, 12))
    for c in range(NCc):
        for i in range(nb):
            sel = np.r_[0:i * k, (i + 1) * k:n]
            Gv, Ph, Pip, Hv = rec[sel, c, 0], rec[sel, c, 1], rec[sel, c, 2], rec[sel, c, 3]
            Fv = Gv + Ph
            mH, mPi = Hv.mean(), Pip.mean()
            vG, vF = Gv.var(ddof=1), Fv.var(ddof=1)
            vP = Ph.var(ddof=1)
            cov = (Gv * Ph).mean() - Gv.mean() * Ph.mean()
            jk[i, c] = [mH, mPi, vG, vF, cov, mH - vG, 0, 0, 0, 0, 0, Gv.mean() + Ph.mean()]
    sj = np.sqrt((nb - 1) / nb * ((jk - jk.mean(0)) ** 2).sum(0))       # (NCc,12)
    full_jk = per.copy()
    s_jk = np.sqrt((sj ** 2).sum(0)) / NCc
    # sigma_tau with delta method for A = E[H]-Var(G)
    st = np.zeros((NCc, 12))
    for c in range(NCc):
        iH, iZ, covHZ = taus[c, 4], taus[c, 2], taus[c, 5]
        sH = np.sqrt(rec[:, c, 3].var(ddof=1) * iH / n)
        z = (rec[:, c, 0] - rec[:, c, 0].mean()) ** 2
        sV = np.sqrt(z.var(ddof=1) * iZ / n)
        sA = np.sqrt(max(sH ** 2 + sV ** 2 - 2 * covHZ, 1e-20))
        Fv = rec[:, c, 0] + rec[:, c, 1]
        zF = (Fv - Fv.mean()) ** 2
        sVF = np.sqrt(zF.var(ddof=1) * taus[c, 3] / n)
        sPi = np.sqrt(rec[:, c, 2].var(ddof=1) * tau_geyer(rec[:, c, 2], ml)[1] / n)
        sCov = np.sqrt(2.0) * 0  # not needed
        st[c] = [sH, sPi, sV, sVF, sCov, sA, 0, 0, 0, 0, 0, 0]
    s_tau = np.sqrt((st ** 2).sum(0)) / NCc
    mean = per.mean(0)
    return dict(names=names, per=per, taus=taus, mean=mean, s_jk=s_jk, s_tau=s_tau,
                neff=(n / np.maximum(per[:, 6], 1e-9)).mean())

# ------------------------------------------------------------------ main
if __name__ == "__main__":
    P("=" * 104)
    P(f"[P02-DYNDYN] Delta_dyn = E[H_vv]-Var(G_v) on L={L} SU(2) beta={BETA}; "
      f"{NC} chains x {NTRAJ} sweeps, one job   {time.strftime('%Y-%m-%d %H:%M:%S')}")
    P("model: o13/o27 constrained xi=0 HMC + O24 fundamental-domain reject |x_l|<pi")
    P("=" * 104)
    dev = checks()
    rng = np.random.default_rng(11)
    eps = EPS0
    for _ in range(5):
        rec_t, acc_t, _, _ = hmc(15, NSTEP, eps, 300, 5)
        P(f"  eps tune {eps:.4f} -> acc={acc_t.mean():.3f}")
        if 0.6 < acc_t.mean() < 0.9:
            break
        eps *= np.exp(np.clip((acc_t.mean() - 0.78) / 3.0, -0.6, 0.6))
    t0 = time.time()
    rec, acc, rfrac, wall = hmc(NTRAJ, NSTEP, eps, 3101, min(THERM, NTRAJ // 5))
    P(f"\n  production: n={rec.shape[0]} sweeps x {NC} chains, acc={acc.mean():.3f} "
      f"(per chain {' '.join(f'{a:.3f}' for a in acc)}), wall={wall:.1f}s, "
      f"{(wall / rec.shape[0] * 1000):.1f} ms/sweep(all chains)")
    P(f"  fundamental-domain rejections (all chains, mean) : {float(rfrac):.4f}")
    st = stats_all(rec)
    P(f"\n  {'quantity':<16s} " + "".join(f"{'ch'+str(c):>16s}" for c in range(NC))
      + f"{'COMBINED':>18s}   {'sigma_jk':>9s} {'sigma_tau':>9s}")
    for j, nm in enumerate(st['names']):
        P(f"  {nm:<16s} " + "".join(f"{st['per'][c, j]:16.6f}" for c in range(NC))
          + f"{st['mean'][j]:18.6f}   {st['s_jk'][j]:9.6f} {st['s_tau'][j]:9.6f}")
    P(f"\n  tau per chain (H, G, z=(G-mu)^2, F=G+Phi', 1+2tau_H):")
    for c in range(NC):
        P(f"    ch{c}: tau_H={st['taus'][c,0]:7.1f}  tau_G={st['taus'][c,1]:7.1f}  "
          f"tau_z={st['taus'][c,2]:7.1f}  tau_F={st['taus'][c,3]:7.1f}  "
          f"1+2tau_H={st['taus'][c,4]:6.2f}  Neff_H~{rec.shape[0]/max(st['taus'][c,4],1e-9):6.0f}")
    A_jk, A_tau = st['mean'][5], st['s_jk'][5]
    P(f"\n  Delta_dyn = E[H_vv]-Var(G_v) = {A_jk:+.6f}")
    P(f"    sigma_jk  = {st['s_jk'][5]:.6f}  -> Delta/sigma = {A_jk/st['s_jk'][5]:+.2f} sigma")
    P(f"    sigma_tau = {st['s_tau'][5]:.6f}  -> Delta/sigma = {A_jk/st['s_tau'][5]:+.2f} sigma")
    P(f"  consistent: B = E[H+Phi'']-Var(G+Phi') = {st['mean'][0]+st['mean'][1]-st['mean'][3]:+.6f}"
      f" +- {np.sqrt(st['s_jk'][0]**2+st['s_jk'][1]**2+st['s_jk'][3]**2):.6f} (jk)")
    P(f"  identity |A+W_full-B| = "
      f"{abs(st['mean'][5] + (st['mean'][1]+st['mean'][2]-st['mean'][3]) - (st['mean'][0]+st['mean'][1]-st['mean'][3])):.2e}")
    P(f"  E[G+Phi'] = {st['mean'][11]:+.6f} +- {st['s_jk'][11]:.6f} (must contain 0)")
    P(f"total elapsed {time.time()-T0:.1f}s   [chek dev {dev:.1e}]")
