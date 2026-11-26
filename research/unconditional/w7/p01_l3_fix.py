"""P01-L3-REFIX: re-run Delta(L=3) with o24_lat.py's fundamental-domain constraint.

Fixes w.r.t. o25_l3.py (see P01-L3.md sec.1):
  (D1) o25 had NO |x_l| < 2*pi restriction -> chain escapes into repeated coverings.
       Here: hard reject any proposal with |x_l| >= 2*pi (DOM2 = (2pi)^2), rejection counted.
  (D2) o25's proposal was  x[l] += d ; x[:,2] -= d_z/NL  (collective z shift on ALL links)
       while dS only re-evaluated the touched link/plaquettes -> the collective part of the
       energy was silently dropped (stale-energy sampler = free drift of the z-centre, the
       analogue of o24 (F2)).  Replaced by o24's PAIR-TRANSFER:  x_l += w, x_l2 -= w
       (w ~ U(-eps,eps)^3), which conserves <v,X>=0 exactly with zero cross-talk, so the
       maintained plaquette traces stay exact -> dS is the EXACT energy difference.

Model identical to o24/o25:  V = beta sum_P (1 - 1/2 ReTr U_P) + sum_l Phi(r_l),
Phi = -log J = -2 log|sin(r/2)/(r/2)|,  U_l = exp(i x_l.sigma/2),  mode A (k=0, |v|=1).
Observable identical: kappa = d_v^2 V, F = d_v V, Delta = E[kappa]-Var(F) = G''(0).

One job, one process: G chains are a leading numpy axis (NOT G processes).

usage: py -3.12 p01_l3_fix.py <L> <nsweep> <eps> <G> <seed0> <burn> <tag>
"""
import sys, time, math
import numpy as np
import o11_physical as o

BETA = 2.5
H = 1e-3
TWO_PI = 2.0 * np.pi
DOM2 = TWO_PI ** 2


def quat(x):
    r = np.sqrt((x * x).sum(-1))
    w = np.cos(0.5 * r)
    s = np.where(r > 1e-14, np.sin(0.5 * r) / np.where(r > 1e-14, r, 1.0), 0.5)
    return w, x * s[..., None]


def phi_r2(r2):
    r = np.sqrt(r2)
    q = np.where(r > 1e-14, np.sin(0.5 * r) / np.where(r > 1e-14, 0.5 * r, 1.0), 1.0)
    return -2.0 * np.log(np.abs(q))


def htq(w, a, idx):
    """(1/2)Re Tr(U0 U1 U2^dag U3^dag) from quaternions.

    w,a: (G, NL) ; idx: broadcastable to (G, n, 4) link ids -> returns (G, n)."""
    ar = np.arange(w.shape[0])[:, None]
    idx = np.asarray(idx)
    if idx.ndim == 2:
        idx = idx[None]                    # (n,4) -> (1,n,4), broadcast over chains
    i0, i1, i2, i3 = idx[..., 0], idx[..., 1], idx[..., 2], idx[..., 3]
    w0, a0 = w[ar, i0], a[ar, i0]
    w1, a1 = w[ar, i1], a[ar, i1]
    w2, b2 = w[ar, i2], a[ar, i2]
    w3, b3 = w[ar, i3], a[ar, i3]
    a2, a3 = -b2, -b3
    tw = w0 * w1 - (a0 * a1).sum(-1)
    ta = w0[..., None] * a1 + w1[..., None] * a0 - np.cross(a0, a1)
    sw = w2 * w3 - (a2 * a3).sum(-1)
    sa = w2[..., None] * a3 + w3[..., None] * a2 - np.cross(a2, a3)
    return tw * sw - (ta * sa).sum(-1)


def parts(x, tab, G, NL):
    """(S_W, Phi) with a leading chain axis."""
    w, a = quat(x)
    tr = htq(w, a, tab)
    return BETA * (len(tab) - tr.sum(-1)), phi_r2((x * x).sum(-1)).sum(-1), tr


def geyer_tau(x, ml=400):
    y = x - x.mean()
    vv = (y * y).sum()
    if vv <= 0:
        return 1.0
    ml = min(ml, len(y) - 3)
    r = np.array([1.0] + [float((y[:-k] * y[k:]).sum() / vv) for k in range(1, ml + 1)])
    if len(r) % 2:
        r = r[:-1]
    s = np.maximum(r[0::2] + r[1::2], 0.0)
    tot = 0.0
    for si in s:
        if si <= 0:
            break
        tot += si
    return max(2.0 * tot - 1.0, 1.0)


def jk_blocks(K, F, nblk=20):
    n = len(K)
    idx = np.array_split(np.arange(n), nblk)

    def dl(sel):
        return K[sel].mean() - F[sel].var()

    d0 = dl(slice(None))
    ds = np.array([dl(np.concatenate([idx[j] for j in range(nblk) if j != i])) for i in range(nblk)])
    sig = math.sqrt((nblk - 1.0) / nblk * ((ds - ds.mean()) ** 2).sum())
    return d0, sig


def run(L, ns, eps, G, seed0, burn):
    tab, near = o.lat_tables(L)
    NL = 4 * L ** 4
    NPQ = len(tab)
    tabI = tab.astype(np.intp)
    ROWS = np.array([near[l] for l in range(NL)], dtype=np.intp)     # (NL, 6)
    v = np.zeros((NL, 3))
    v[:, 2] = 1.0 / math.sqrt(NL)
    rng = np.random.default_rng([seed0 + i for i in range(G)])
    ar = np.arange(G)

    x = rng.normal(0, 0.55, (G, NL, 3))
    x -= np.einsum('gla,la->g', x, v)[:, None, None] * v
    w, a = quat(x)                                   # maintained (G,NL), (G,NL,3)
    tr = htq(w, a, tabI)                             # maintained (G,NPQ)
    phS = phi_r2((x * x).sum(-1)).sum(-1)

    # --- validation of the quaternion energy vs o11_physical (same check o25 reported) ---
    dv = float(np.abs(BETA * (NPQ - tr.sum(-1)) + phS - np.array(
        [o.Vtot(x[g], tab, BETA) for g in range(G)])).max())
    trv = float(np.abs(tr[0] - o.tr_plaq(o.expi(x[0]), tab, np.arange(NPQ))).max())
    print(f"[P01] validate: |dV|={dv:.2e}  |dtr|={trv:.2e}  (L={L}, {NL} links, {NPQ} plaq)",
          flush=True)

    nacc = np.zeros(G)
    ndom = 0
    Ks = np.empty((ns, G)); Fs = np.empty((ns, G))
    KA = np.empty((ns, G)); KB = np.empty((ns, G))
    FA = np.empty((ns, G)); FB = np.empty((ns, G))
    cache_dev = []
    # precompute: which plaquette row (if any) is shared by the pair (l,l2)
    sh = (ROWS[:, None, :, None] == ROWS[None, :, None, :]).any(-1)          # (NL,NL,6)
    SHROW = np.where(sh.any(-1), sh.argmax(-1), -1).astype(np.intp)          # (NL,NL)
    aoff = ar * NL
    ar2 = ar[:, None]
    arbase = aoff[:, None, None]
    wf = w.reshape(-1)
    af = a.reshape(-1, 3)
    ph2pi = 2.0 * np.pi
    t0 = time.time()
    for sw in range(ns + burn):
        for _ in range(NL):
            l = rng.integers(NL, size=G)
            l2 = (l + rng.integers(1, NL, size=G)) % NL
            wv = rng.uniform(-eps, eps, (G, 3))
            rows = np.concatenate([ROWS[l], ROWS[l2]], axis=1)               # (G,12)
            xo1, xo2 = x[ar, l], x[ar, l2]
            xn1, xn2 = xo1 + wv, xo2 - wv
            r1 = (xn1 * xn1).sum(-1); r2 = (xn2 * xn2).sum(-1)
            q1 = np.sqrt(r1); q2 = np.sqrt(r2)
            dom = (q1 < TWO_PI) & (q2 < TWO_PI)
            wl, al = np.cos(0.5 * q1), xn1 * (0.5 * np.sinc(q1 / ph2pi))[:, None]
            w2, a2 = np.cos(0.5 * q2), xn2 * (0.5 * np.sinc(q2 / ph2pi))[:, None]
            i1, i2 = aoff + l, aoff + l2
            wo1, ao1 = wf[i1].copy(), af[i1].copy()
            wo2, ao2 = wf[i2].copy(), af[i2].copy()
            wf[i1], af[i1] = wl, al
            wf[i2], af[i2] = w2, a2
            ic = tabI[rows]                                                  # (G,12,4)
            icf = ic + arbase
            wc = wf[icf]; ac = af[icf]                                       # (G,12,4)/(G,12,4,3)
            c0w, c1w, c2w, c3w = wc[..., 0], wc[..., 1], wc[..., 2], wc[..., 3]
            c0a, c1a = ac[..., 0, :], ac[..., 1, :]
            c2a, c3a = -ac[..., 2, :], -ac[..., 3, :]
            tw = c0w * c1w - (c0a * c1a).sum(-1)
            ta = c0w[..., None] * c1a + c1w[..., None] * c0a - np.cross(c0a, c1a)
            sw_ = c2w * c3w - (c2a * c3a).sum(-1)
            sa = c2w[..., None] * c3a + c3w[..., None] * c2a - np.cross(c2a, c3a)
            new = tw * sw_ - (ta * sa).sum(-1)                               # (G,12)
            old = tr[ar2, rows]
            vals = old - new
            sg = SHROW[l, l2]
            corr = np.where(sg >= 0, vals[ar, np.maximum(sg, 0)], 0.0)
            dE = vals.sum(1) - corr
            # phi = -2 log|sin(q/2)/(q/2)| = -2 log|sinc(q/2pi)|
            phi4 = -2.0 * np.log(np.abs(np.sinc(np.stack([q1, q2], 1) / ph2pi)))
            phn = (phS - phi_r2((xo1 * xo1).sum(-1)) - phi_r2((xo2 * xo2).sum(-1))
                   + phi4[:, 0] + phi4[:, 1])
            dS = BETA * dE + (phn - phS)
            m = (dS <= 0) | (rng.random(G) < np.exp(-np.maximum(dS, 0.0)))
            ndom += int((~dom).sum())
            m &= dom
            nacc += m
            gs = np.nonzero(m)[0]
            if len(gs):
                x[gs, l[gs]] = xn1[gs]
                x[gs, l2[gs]] = xn2[gs]
                phS[gs] = phn[gs]
                tr[gs[:, None], rows[gs]] = new[gs]
            rej = np.nonzero(~m)[0]
            if len(rej):
                j1, j2 = aoff[rej] + l[rej], aoff[rej] + l2[rej]
                wf[j1], af[j1] = wo1[rej], ao1[rej]
                wf[j2], af[j2] = wo2[rej], ao2[rej]
        if sw < burn:
            continue
        j = sw - burn
        if j % 1000 == 0:                     # integrity: cached energy vs full recompute
            Swf, Phf, trf = parts(x, tabI, G, NL)
            c = float(np.abs((BETA * (NPQ - tr.sum(-1)) + phS) - (Swf + Phf)).max())
            cache_dev.append(c)
            print(f"   [check] sweep={sw} cache-full={c:.2e} |dtr|={np.abs(trf-tr).max():.2e} "
                  f"max_r={np.linalg.norm(x, axis=-1).max():.3f}", flush=True)
        xp = x + H * v[None]
        xm = x - H * v[None]
        Swp, Php, _ = parts(xp, tabI, G, NL)
        Swm, Phm, _ = parts(xm, tabI, G, NL)
        Sw0 = BETA * (NPQ - tr.sum(-1))
        Ph0 = phS
        h2 = H * H
        kA = (Swp - 2 * Sw0 + Swm) / h2
        kB = (Php - 2 * Ph0 + Phm) / h2
        FA[j] = (Swp - Swm) / (2 * H)
        FB[j] = (Php - Phm) / (2 * H)
        KA[j] = kA; KB[j] = kB
        Ks[j] = kA + kB
        Fs[j] = FA[j] + FB[j]
    el = time.time() - t0
    return dict(Ks=Ks, Fs=Fs, KA=KA, KB=KB, FA=FA, FB=FB,
                acc=nacc / ((ns + burn) * NL), domrej=ndom / (G * (ns + burn) * NL),
                el=el, NL=NL, NPQ=NPQ, x=x, cache_dev=max(cache_dev) if cache_dev else -1.0)


if __name__ == "__main__":
    L = int(sys.argv[1]); ns = int(sys.argv[2]); eps = float(sys.argv[3])
    G = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    seed0 = int(sys.argv[5]) if len(sys.argv) > 5 else 500
    burn = int(sys.argv[6]) if len(sys.argv) > 6 else 800
    tag = sys.argv[7] if len(sys.argv) > 7 else "FX"
    t0 = time.time()
    r = run(L, ns, eps, G, seed0, burn)
    Ks, Fs = r['Ks'], r['Fs']
    rmax = np.linalg.norm(r['x'], axis=-1).max()
    print(f"[P01] L={L} tag={tag} ns={ns} burn={burn} eps={eps} G={G} "
          f"seeds={[seed0+i for i in range(G)]} NL={r['NL']} NPQ={r['NPQ']}", flush=True)
    print(f"[P01] wall={r['el']:.1f}s ({1e3*r['el']/(ns+burn):.2f} ms/sweep for {G} chains) "
          f"acc_link={np.round(r['acc'],3)} domain-reject={r['domrej']:.4f} "
          f"final r: mean={np.linalg.norm(r['x'],axis=-1).mean():.3f} max={rmax:.3f} "
          f"(2pi={2*np.pi:.3f})", flush=True)
    ds, sgs = [], []
    for g in range(G):
        K, F = Ks[:, g], Fs[:, g]
        d0, sg = jk_blocks(K, F, 20)
        tk = geyer_tau(K); tf = geyer_tau(F)
        semK = math.sqrt(K.var(ddof=1) * (1 + 2 * tk) / ns)
        semV = math.sqrt(2 * F.var(ddof=1) ** 2 * (1 + 2 * tf) / ns)
        st = math.hypot(semK, semV)
        ds.append(d0); sgs.append(sg)
        print(f"   g={g} seed={seed0+g} acc={r['acc'][g]:.3f} E[kappa]={K.mean():+.4f} "
              f"Var(F)={F.var():.4f} E[F]={F.mean():+.4f} Delta={d0:+.4f} "
              f"sig_jk={sg:.4f} sig_tau={st:.4f} tauK={tk:.1f} tauF={tf:.1f} "
              f"Neff={ns/max(tk,tf):.0f} E[H]={r['KA'][:,g].mean():+.4f} "
              f"E[Phipp]={r['KB'][:,g].mean():+.4f} Var(Fw)={r['FA'][:,g].var():.4f} "
              f"Var(Fphi)={r['FB'][:,g].var():.4f}", flush=True)
    ds = np.array(ds); sgs = np.array(sgs)
    Dm = ds.mean()
    sjk = math.sqrt((sgs ** 2).sum()) / G
    # tau-error: rebuild per chain then combine
    sts = []
    for g in range(G):
        K, F = Ks[:, g], Fs[:, g]
        tk = geyer_tau(K); tf = geyer_tau(F)
        semK = math.sqrt(K.var(ddof=1) * (1 + 2 * tk) / ns)
        semV = math.sqrt(2 * F.var(ddof=1) ** 2 * (1 + 2 * tf) / ns)
        sts.append(math.hypot(semK, semV))
    sts = np.array(sts)
    stau = math.sqrt((sts ** 2).sum()) / G
    sdb = ds.std(ddof=1) / math.sqrt(G)
    print(f"[P01] COMBINED G={G}: E[kappa]={Ks.mean():+.4f} Var(F)={Fs.var():.4f} "
          f"E[F]={Fs.mean():+.4f} E[H]={r['KA'].mean():+.4f} E[Phi'']={r['KB'].mean():+.4f} "
          f"Delta={Dm:+.4f} sig_jk={sjk:.4f} ({Dm/sjk:+.2f}s) sig_tau={stau:.4f} "
          f"({Dm/stau:+.2f}s) sd_between={sdb:.4f} Var(F)/chain std="
          f"{np.std([Fs[:,g].var() for g in range(G)],ddof=1):.4f}", flush=True)
    print(f"[P01] L2=L^4={L**4} Delta_norm=L^2*Delta={L**2*Dm:+.4f} "
          f"total_wall={time.time()-t0:.1f}s", flush=True)
    np.savez(f"p01_{tag}.npz", K=Ks, F=Fs, KA=r['KA'], KB=r['KB'], FA=r['FA'], FB=r['FB'],
             acc=r['acc'], domrej=r['domrej'], ns=ns, eps=eps, burn=burn, L=L,
             seeds=np.array([seed0 + i for i in range(G)]))
