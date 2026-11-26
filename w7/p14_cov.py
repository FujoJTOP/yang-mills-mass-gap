#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""p14_cov.py -- sign and bound of Cov(G_v, Phi') at L=2, SU(2), beta=2.5.

Reuses o24's sampler VERBATIM (research/unconditional/w7/o24_lat.py: flat chart on the
SU(2) fundamental domain, pair-transfer Metropolis, mode A = one common v per link,
so <v,X> = 0 is conserved to machine precision).  What is new here is only the
*recording*: on the same configurations we split d_v V into its action part and its
Jacobian part,

    G_v  = d_v S_W = (S_W(x+hv)-S_W(x-hv))/2h
    Phi' = d_v Phi = (Phi(x+hv)-Phi(x-hv))/2h
    Phi''= d_v^2 Phi,  H_vv = d_v^2 S_W,  kappa = H_vv+Phi'' = d_v^2 V
and, as a CONTROL, the same bilinear along a fixed direction u _|_ v:
the constraint plane is {<v,X>=0}, so <u,X> is NOT conserved and u is a
"unconstrained" control direction on the identical ensemble.

usage: py -3.12 p14_cov.py <nsweep> <eps> [seed0] [burn] [tag]
"""
import sys, time
import numpy as np

sys.path.insert(0, r"D:\Dev\math\research\unconditional\w7")
import o24_lat as m

BETA, NL, NPQ = m.BETA, m.NL, m.NPQ
H = 1e-3
P = lambda *a: print(*a, flush=True)


def run(v, u, G, ns, eps, rng, burn=800, h=H):
    assert np.allclose(v, v[0])
    x = rng.normal(0, 0.55, (G, NL, 3))
    x -= np.einsum('gla,la->g', x, v)[:, None, None] * v[None]
    w, a = m.quat(x)
    tr = m.ht(w, a, m.C0, m.C1, m.C2, m.C3)
    phS = m.phi_r2((x * x).sum(-1)).sum(-1)
    ar = np.arange(G)
    nacc = np.zeros(G)
    rec = np.zeros((ns, G, 8))
    for sw in range(ns + burn):
        for _ in range(NL):
            l = rng.integers(NL, size=G)
            l2 = (l + rng.integers(1, NL, size=G)) % NL
            wv = rng.uniform(-eps, eps, (G, 3))
            xn = x.copy()
            xn[ar, l] += wv
            xn[ar, l2] -= wv
            r2n = (xn * xn).sum(-1)
            phn = m.phi_r2(r2n).sum(-1)
            wn, an = m.quat(xn)
            rows = np.concatenate([m.ROWS[l], m.ROWS[l2]], axis=1)
            new = m.ht(wn, an, m.C0[rows], m.C1[rows], m.C2[rows], m.C3[rows])
            old = tr[ar[:, None], rows]
            dS = BETA * (old - new).sum(-1) + (phn - phS)
            mk = (dS <= 0) | (rng.random(G) < np.exp(-np.maximum(dS, 0.0)))
            dom = (r2n < m.DOM2).all(-1)
            mk &= dom
            nacc += mk
            if mk.any():
                x[mk] = xn[mk]
                w[mk] = wn[mk]
                a[mk] = an[mk]
                phS[mk] = phn[mk]
                trn = tr.copy()
                trn[ar[:, None], rows] = new
                tr[mk] = trn[mk]
        if sw < burn:
            continue
        j = sw - burn
        sw0 = BETA * (NPQ - tr.sum(-1))
        V0 = sw0 + phS
        vv = v[None]
        Vp, trp, php = m.full_V(x + h * vv)
        Vm, trm, phm = m.full_V(x - h * vv)
        swp = BETA * (NPQ - trp.sum(-1))
        swm = BETA * (NPQ - trm.sum(-1))
        rec[j, :, 0] = (swp - swm) / (2 * h)                 # G_v
        rec[j, :, 1] = (php - phm) / (2 * h)                 # Phi'
        rec[j, :, 2] = (Vp - 2 * V0 + Vm) / h ** 2           # kappa
        rec[j, :, 3] = (php - 2 * phS + phm) / h ** 2        # Phi''
        rec[j, :, 4] = (swp - 2 * sw0 + swm) / h ** 2        # H_vv
        if u is not None:
            uu = u[None]
            Up, trq, phq = m.full_V(x + h * uu)
            Um, tru, phu = m.full_V(x - h * uu)
            rec[j, :, 5] = (BETA * (NPQ - trq.sum(-1)) - BETA * (NPQ - tru.sum(-1))) / (2 * h)
            rec[j, :, 6] = (phq - phu) / (2 * h)
        rec[j, :, 7] = np.einsum('gla,la->g', x, v)          # <v,X>  (must be 0)
    return rec, nacc / ((ns + burn) * NL), x


def acf(y, ml=2000):
    y = y - y.mean()
    vv = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / vv for k in range(1, min(ml, len(y) - 2))])


def geyer(y, ml=2000):
    r = np.concatenate(([1.0], acf(y, ml)))
    if len(r) % 2:
        r = r[:-1]
    s = np.maximum(r[0::2] + r[1::2], 0.0)
    tot = 0.0
    for si in s:
        if si <= 0:
            break
        tot += si
    return max(2 * tot - 1.0, 1.0)


def jk_sigma(f, Y, nb):
    """jackknife sigma of functional f(Y, sel) over nb contiguous blocks"""
    n = Y.shape[0]
    k = n // nb
    full = f(Y, np.arange(n))
    acc = []
    for i in range(nb):
        sel = np.r_[0:i * k, (i + 1) * k:n].astype(int)
        acc.append(f(Y, sel))
    acc = np.array(acc)
    return full, np.sqrt((nb - 1) / nb * ((acc - acc.mean(0)) ** 2).sum(0))


def stat(y, sel):
    """y: (n, 4) columns Gv, PhiP, PhiPP, Hvv -> [Cov, Var G, Var P, Var(G+P), E Phi'', E Hvv, corr]"""
    y = y[sel]
    G, Pp, Pi2, Hv = y[:, 0], y[:, 1], y[:, 2], y[:, 3]
    mG, mP = G.mean(), Pp.mean()
    vG, vP = G.var(), Pp.var()
    cov = (G * Pp).mean() - mG * mP
    return np.array([cov, vG, vP, (G + Pp).var(), Pi2.mean(), Hv.mean(), mG, mP,
                     cov / np.sqrt(vG * vP)])


if __name__ == "__main__":
    ns = int(sys.argv[1])
    eps = float(sys.argv[2])
    seed0 = int(sys.argv[3]) if len(sys.argv) > 3 else 440
    burn = int(sys.argv[4]) if len(sys.argv) > 4 else 800
    tag = sys.argv[5] if len(sys.argv) > 5 else 'cov'
    G = 4
    ctrl = not (len(sys.argv) > 6 and sys.argv[6] == 'nou')
    v = np.zeros((NL, 3)); v[:, 2] = 1.0 / np.sqrt(NL)
    u = np.zeros((NL, 3)); u[:, 0] = 1.0 / np.sqrt(NL)      # u _|_ v, unconstrained
    assert abs(v[0] @ u[0]) < 1e-15
    if not ctrl:
        u = None
    seeds = [seed0 + i for i in range(G)]
    t0 = time.time()
    rec, acc, x = run(v, u, G, ns, eps, np.random.default_rng(seeds), burn=burn)
    el = time.time() - t0
    np.savez(f"p14_{tag}.npz", rec=rec, seeds=np.array(seeds), ns=ns, eps=eps, acc=acc, el=el)
    P(f"[P14] L=2 SU(2) beta={BETA} mode=A  nsw={ns} burn={burn} eps={eps} G={G} "
      f"seeds={seeds} wall={el:.1f}s ({1e3*el/(ns+burn):.2f} ms/sweep all 4 chains)")
    P(f"[P14] acc_link={np.round(acc,3)}  |<v,X>|max={np.abs(rec[:,:,7]).max():.2e} (exact constraint)")
    names = ["Cov", "VarG", "VarPhi'", "Var(G+Phi')", "E[Phi'']", "E[H_vv]", "E[G]", "E[Phi']", "corr"]
    allsel = np.arange(ns)
    for g in range(G):
        s = stat(rec[:, g, :4], allsel)
        P(f"  g={g} seed={seeds[g]} E[kappa]={s[5]+s[4]:+.4f} " +
          "  ".join(f"{n}={s[i]:+.4f}" for i, n in enumerate(names)))
        sc = stat(rec[:, g, [5, 6, 2, 4]], allsel)
        if ctrl:
            P(f"      control u (u_|_v, unconstrained): Cov(G_u,Phi'_u)={sc[0]:+.4f} "
              f"VarGu={sc[1]:.3f} VarP'u={sc[2]:.3f} corr={sc[8]:+.3f}")
    # ---- pooled (no centring: stat subtracts the means, so E[G], E[Phi'] stay as the parity test) --
    Y = rec[:, :, :4].reshape(-1, 4)
    full, sig = jk_sigma(stat, Y, 25 * G)
    Yu = np.zeros((G * ns, 4))
    Yu[:, :2] = rec[:, :, [5, 6]].reshape(-1, 2)
    if ctrl:
        fullu, sigu = jk_sigma(stat, Yu, 25 * G)
    else:
        fullu = sigu = np.zeros(9)
    # ---- autocorrelation errors on the Cov series ----
    zg = (Y[:, 0] * Y[:, 1])
    tG, tP, tZ = geyer(Y[:, 0]), geyer(Y[:, 1]), geyer(zg)
    nev = zg.size / max(tZ, 1.0)
    s_tau = np.sqrt(zg.var() / nev)
    tZu, s_tau_u = 0.0, 0.0
    if ctrl:
        tZu = geyer(Yu[:, 0] * Yu[:, 1])
        s_tau_u = np.sqrt((Yu[:, 0] * Yu[:, 1]).var() / (Yu.shape[0] / max(tZu, 1.0)))
    half = ns // 2
    c1 = stat(rec[:half, :, :4].reshape(-1, 4), np.arange(half * G))[0]
    c2 = stat(rec[half:, :, :4].reshape(-1, 4), np.arange((ns - half) * G))[0]
    per = np.array([stat(rec[:, g, :4], allsel)[0] for g in range(G)])
    P("")
    P(f"[POOLED] n={Y.shape[0]} samples (4 chains x {ns})  tau_G={tG:.1f} tau_Phi'={tP:.1f} "
      f"tau_(G Phi')={tZ:.1f}  Neff={nev:.0f}")
    for i, n in enumerate(names):
        P(f"  {n:<14s} = {full[i]:+10.6f}  +- {sig[i]:.6f} (jk, {25*G} blocks)   "
          f"[tau-series sigma: {'%.6f' % s_tau if i == 0 else '-'}]")
    P(f"  Cov/3-sigma: {full[0]/sig[0]:+.2f} (jk)   {full[0]/s_tau:+.2f} (tau)")
    P(f"  per-chain Cov = {np.round(per,5)}  mean={per.mean():+.6f} "
      f"grp_err={per.std(ddof=1)/np.sqrt(G):.6f} spread={per.std(ddof=1):.6f}")
    P(f"  window stability: Cov(1st half)={c1:+.5f}   Cov(2nd half)={c2:+.5f}")
    P("")
    if ctrl:
        P(f"[CONTROL u_|_v, same ensemble, <u,X> NOT constrained]  tau_(GuP'u)={tZu:.1f}")
        cn = ["Cov_u", "VarGu", "VarP'u", "Var(Gu+P'u)", "E[Phi''_u]", "E[H_u]", "E[G_u]", "E[Phi'_u]", "corr_u"]
        for i, n in enumerate(cn):
            P(f"  {n:<14s} = {fullu[i]:+10.6f}  +- {sigu[i]:.6f}   "
              f"[tau sigma: {'%.6f' % s_tau_u if i == 0 else '-'}]")
        P(f"  Cov_u/sigma: {fullu[0]/sigu[0]:+.2f} (jk)  {fullu[0]/s_tau_u:+.2f} (tau)")
    # ---- bookkeeping in o24 normalisation ----
    Wful = full[4] + full[1] - full[3]
    P("")
    P(f"[BOOKKEEP] W_full = E[Phi''] + VarG - Var(G+Phi') = {Wful:+.4f}   "
      f"= E[Phi''] - VarPhi' - 2Cov = {full[4] - full[2] - 2*full[0]:+.4f}")
    P(f"  cross part -VarPhi'-2Cov = {-(full[2] + 2 * full[0]):+.4f}   "
      f"C-S bound |Cov| <= sigma_G sigma_Phi' = {np.sqrt(full[1]*full[2]):.4f}")
    P(f"total elapsed {time.time()-t0:.1f}s")
