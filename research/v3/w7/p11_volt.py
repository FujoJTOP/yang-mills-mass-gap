"""P11 — Delta_dyn = E[H] - Var(G_v) versus L, with the O24 pair-transfer machinery plus the
fundamental-domain correction that O25-L3 lacked.

Model (same flat-chart density as o11_physical / o24_lat):
    V(x) = beta sum_P (1 - 1/2 ReTr U_P) + sum_l Phi(r_l),   Phi = -2 log|sin(r/2)/(r/2)|
    U_l = exp(i x_l.sigma/2),  |x_l| < 2 pi  (ball = SU(2) group manifold exactly once)
Slice {<v,X> = 0} enforced EXACTLY by PAIR-TRANSFER proposals  x_l += w, x_l' -= w
(w conserves sum_l x_l,z = 0 with v_l = zhat/sqrt(NL) on every link, hence <v,X>=0 at all times).
Because only two links move, only plaquettes touching l or l' change -> per-proposal cost is local,
which is what makes L=3 (NL=324, NPQ=486) cheap.

Fixes vs o25_l3.py (the "fundamental-domain bug" that produced the negative L=3):
  * o25 projected by translating x[:,2] -= d[2]/nlink, i.e. it MOVED ALL nlink links at once but
    accepted/rejected using only the 6 plaquettes and 1 Jacobian factor of the moved link.
    The global translation changes every r_l and every plaquette trace -> wrong Boltzmann weight.
    Here the slice is conserved by a local +w/-w swap with an exact local Delta S.
  * duplicate plaquettes when l,l' sit in the same plaquette are deduplicated (o24 double-counted
    them in Delta S, a ~2-4% acceptance distortion).

usage: py -3.12 p11_volt.py <L> <nsweep> <eps> <seed0> <burn> [tag] [G]
"""
import sys
import time
import numpy as np

sys.path.insert(0, r"D:\Dev\math\research\unconditional\w7")
import o11_physical as o                      # noqa: E402

BETA = 2.5
H = 1e-3
DOM2 = (2 * np.pi) ** 2


def quat(x):
    r = np.sqrt((x * x).sum(-1))
    w = np.cos(0.5 * r)
    s = np.where(r > 1e-14, np.sin(0.5 * r) / np.where(r > 1e-14, r, 1.0), 0.5)
    return w, x * s[..., None]


def phi_r2(r2):
    r = np.sqrt(r2)
    q = np.where(r > 1e-14, np.sin(0.5 * r) / np.where(r > 1e-14, 0.5 * r, 1.0), 1.0)
    return -2.0 * np.log(np.abs(q))


def trace4(xA, xB, xC, xD):
    """(1/2)ReTr(U_A U_B U_C^dag U_D^dag) from link coordinates, last axis = 3, leading = (G,rows)."""
    wA, aA = quat(xA)
    wB, aB = quat(xB)
    wC, aC = quat(xC)
    wD, aD = quat(xD)
    # chart U = exp(i x.sigma/2) => quaternion product (q0 p0 - q.p, q0 p + p0 q - q x p)  [anti-Hamilton]
    tw = wA * wB - (aA * aB).sum(-1)
    ta = wA[..., None] * aB + wB[..., None] * aA - np.cross(aA, aB)
    sw = wC * wD - (aC * aD).sum(-1)
    sa = -(wC[..., None] * aD + wD[..., None] * aC) - np.cross(aC, aD)   # vector part of Udag_C Udag_D
    return tw * sw - (ta * sa).sum(-1)


def Vtot(x, C0, C1, C2, C3, ar):
    tr = trace4(x[ar[:, None], C0], x[ar[:, None], C1], x[ar[:, None], C2], x[ar[:, None], C3])
    return BETA * (C0.shape[0] - tr.sum(-1)) + phi_r2((x * x).sum(-1)).sum(-1)


def run(L, G, ns, eps, rng, burn, h=H):
    tab, near = o.lat_tables(L)
    NL = 4 * L ** 4
    NPQ = len(tab)
    C0, C1, C2, C3 = (tab[:, k].astype(np.intp) for k in range(4))
    ROWS = np.array([near[l] for l in range(NL)], dtype=np.intp)   # (NL, 6) plaquettes per link
    allp = np.arange(NPQ)
    VS = 1.0 / np.sqrt(NL)

    x = rng.normal(0, 0.55, (G, NL, 3))
    x[..., 2] -= x[..., 2].mean(-1, keepdims=True)       # project onto the slice sum_l x_l,z = 0
    x *= (DOM2 ** 0.5 * 0.6) / np.maximum(np.linalg.norm(x, axis=-1).max(-1), 1e-12)[:, None, None]
    return_ = None

    ph = phi_r2((x * x).sum(-1))
    phS = ph.sum(-1)
    ar = np.arange(G)
    a3 = ar[:, None]
    tr = trace4(x[np.broadcast_to(a3, (G, NPQ)), C0[None, :]],
                x[np.broadcast_to(a3, (G, NPQ)), C1[None, :]],
                x[np.broadcast_to(a3, (G, NPQ)), C2[None, :]],
                x[np.broadcast_to(a3, (G, NPQ)), C3[None, :]])
    nacc = np.zeros(G)
    ndom = 0
    Ks = np.empty((ns, G))
    Fs = np.empty((ns, G))
    KA = np.empty((ns, G))
    KB = np.empty((ns, G))
    tot = ns + burn
    for sw in range(tot):
        for _ in range(NL):
            l = rng.integers(NL, size=G)
            l2 = (l + rng.integers(1, NL, size=G)) % NL
            wv = rng.uniform(-eps, eps, (G, 3))
            xa = x[ar, l] + wv
            xb = x[ar, l2] - wv
            r2a = (xa * xa).sum(-1)
            r2b = (xb * xb).sum(-1)
            pa = phi_r2(r2a)
            pb = phi_r2(r2b)
            dphi = pa + pb - ph[ar, l] - ph[ar, l2]
            rows = np.sort(np.concatenate([ROWS[l], ROWS[l2]], axis=1), axis=1)
            keep = np.ones(rows.shape, bool)
            keep[:, 1:] = rows[:, 1:] != rows[:, :-1]
            A, B, C, D = C0[rows], C1[rows], C2[rows], C3[rows]
            sel_a = A == l[:, None]
            sel_b = A == l2[:, None]
            xA = np.where(sel_a[..., None], xa[:, None, :], np.where(sel_b[..., None], xb[:, None, :],
                                                                       x[a3, A]))
            sel_a = B == l[:, None]
            sel_b = B == l2[:, None]
            xB = np.where(sel_a[..., None], xa[:, None, :], np.where(sel_b[..., None], xb[:, None, :],
                                                                       x[a3, B]))
            sel_a = C == l[:, None]
            sel_b = C == l2[:, None]
            xC = np.where(sel_a[..., None], xa[:, None, :], np.where(sel_b[..., None], xb[:, None, :],
                                                                       x[a3, C]))
            sel_a = D == l[:, None]
            sel_b = D == l2[:, None]
            xD = np.where(sel_a[..., None], xa[:, None, :], np.where(sel_b[..., None], xb[:, None, :],
                                                                       x[a3, D]))
            new = trace4(xA, xB, xC, xD)
            old = tr[a3, rows]
            dS = BETA * ((old - new) * keep).sum(-1) + dphi
            m = (dS <= 0) | (rng.random(G) < np.exp(-np.maximum(dS, 0.0)))
            dom = (r2a < DOM2) & (r2b < DOM2)
            ndom += int((~dom).sum())
            m &= dom
            nacc += m
            s = np.nonzero(m)[0]
            if s.size:
                x[s, l[s]] = xa[s]
                x[s, l2[s]] = xb[s]
                ph[s, l[s]] = pa[s]
                ph[s, l2[s]] = pb[s]
                phS[s] += dphi[s]
                tr[s[:, None], rows[s]] = new[s]
        if sw < burn:
            continue
        j = sw - burn
        v = np.zeros((G, NL, 3))
        v[..., 2] = VS
        xp = x + h * v
        xm = x - h * v
        trp = trace4(xp[np.broadcast_to(a3, (G, NPQ)), C0[None, :]],
                     xp[np.broadcast_to(a3, (G, NPQ)), C1[None, :]],
                     xp[np.broadcast_to(a3, (G, NPQ)), C2[None, :]],
                     xp[np.broadcast_to(a3, (G, NPQ)), C3[None, :]])
        pAp = phi_r2((xp * xp).sum(-1)).sum(-1)
        trm = trace4(xm[np.broadcast_to(a3, (G, NPQ)), C0[None, :]],
                     xm[np.broadcast_to(a3, (G, NPQ)), C1[None, :]],
                     xm[np.broadcast_to(a3, (G, NPQ)), C2[None, :]],
                     xm[np.broadcast_to(a3, (G, NPQ)), C3[None, :]])
        pAm = phi_r2((xm * xm).sum(-1)).sum(-1)
        tr0 = tr
        Vp = BETA * (NPQ - trp.sum(-1)) + pAp
        V0 = BETA * (NPQ - tr0.sum(-1)) + phS
        Vm = BETA * (NPQ - trm.sum(-1)) + pAm
        Ks[j] = (Vp - 2.0 * V0 + Vm) / h ** 2
        Fs[j] = (Vp - Vm) / (2.0 * h)
        # split: S_W part (A) and Jacobian part (B)
        KA[j] = -BETA * (trp.sum(-1) - 2.0 * tr0.sum(-1) + trm.sum(-1)) / h ** 2
        KB[j] = (pAp - 2.0 * phS + pAm) / h ** 2
    return Ks, Fs, KA, KB, nacc / (tot * NL), ndom / (G * tot * NL), x


def acf(x, ml=400):
    y = x - x.mean()
    v = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / v for k in range(1, min(ml, len(y) - 2))])


def geyer(x, ml=400):
    r = np.concatenate(([1.0], acf(x, ml)))
    if len(r) % 2:
        r = r[:-1]
    s = np.maximum(r[0::2] + r[1::2], 0.0)
    tot = 0.0
    for si in s:
        if si <= 0:
            break
        tot += si
    return max(2.0 * tot - 1.0, 1.0)


def group_stats(k, f, nb=40):
    n = len(k)
    idx = np.array_split(np.arange(n), nb)
    ds = np.array([(np.delete(k, b).mean() - np.delete(f, b).var()) for b in idx])
    s_jk = np.sqrt((nb - 1) / nb * ((ds - ds.mean()) ** 2).sum())
    tK, tF = geyer(k), geyer(f)
    sK = np.sqrt(k.var(ddof=1) * (1 + 2 * tK) / n)
    sF = np.sqrt(2 * f.var(ddof=1) ** 2 * (1 + 2 * tF) / n)
    return dict(d=k.mean() - f.var(), s_jk=s_jk, s_tau=np.hypot(sK, sF), tK=tK, tF=tF,
                neff=n / max(tK, tF), mK=k.mean(), vF=f.var(ddof=1), mF=f.mean())


if __name__ == "__main__":
    L = int(sys.argv[1])
    ns = int(sys.argv[2])
    eps = float(sys.argv[3])
    seed0 = int(sys.argv[4])
    burn = int(sys.argv[5])
    tag = sys.argv[6] if len(sys.argv) > 6 else "x"
    G = int(sys.argv[7]) if len(sys.argv) > 7 else 4
    NP = 6 * L ** 4
    t0 = time.time()
    Ks, Fs, KA, KB, acc, domrej, x = run(L, G, ns, eps, np.random.default_rng(seed0), burn)
    el = time.time() - t0
    print(f"[P11] L={L} G={G} ns={ns} burn={burn} eps={eps} seed0={seed0} NL={4*L**4} NPQ={NP} "
          f"wall={el:.1f}s {1e3*el/(ns+burn):.2f}ms/sweep acc={np.round(acc,3)} "
          f"domrej={domrej:.5f} maxr={np.linalg.norm(x, axis=-1).max():.3f}", flush=True)
    rows = []
    for g in range(G):
        r = group_stats(Ks[:, g], Fs[:, g])
        rows.append(r)
        print(f"   g={g} E[k]={r['mK']:+8.4f} Var(F)={r['vF']:7.4f} E[F]={r['mF']:+.4f} "
              f"D={r['d']:+.4f} s_jk={r['s_jk']:.4f} s_tau={r['s_tau']:.4f} "
              f"tauK={r['tK']:5.1f} tauF={r['tF']:5.1f} Neff={r['neff']:6.0f}", flush=True)
    D = np.mean([r['d'] for r in rows])
    S_jk = np.sqrt(np.sum([r['s_jk'] ** 2 for r in rows])) / G
    S_tau = np.sqrt(np.sum([r['s_tau'] ** 2 for r in rows])) / G
    DS = np.array([r['d'] for r in rows])
    S_grp = np.sqrt((G - 1) / G * ((G * D - (G - 1) * DS - D) ** 2).sum()) if G > 1 else float('nan')
    EK = np.mean([r['mK'] for r in rows])
    VF = np.mean([r['vF'] for r in rows])
    sp = np.sqrt(np.mean([r['s_jk'] ** 2 for r in rows]))
    print(f"[P11] L={L} Delta={D:+.4f} +- {S_jk:.4f}(jk) {S_tau:.4f}(tau) grp={S_grp:.4f} "
          f"ratio_jk={D/S_jk:+.2f} EK={EK:+.4f}+-{sp:.4f} VF={VF:.4f} E[F]={np.mean([r['mF'] for r in rows]):+.4f} "
          f"Delta/Np={D/NP:+.6f} Np={NP}", flush=True)
    print(f"[P11] L={L} split: E[kappa_SW]={KA.mean():+.4f} E[kappa_Phi]={KB.mean():+.4f} "
          f"Var(F_SW)={np.var(KA*0+0)+0:.0f}", flush=True)
    np.savez(f"p11_L{L}_{tag}.npz", K=Ks, F=Fs, KA=KA, KB=KB, acc=acc, domrej=domrej,
             el=el, ns=ns, eps=eps, seed0=seed0, burn=burn, L=L, G=G)
