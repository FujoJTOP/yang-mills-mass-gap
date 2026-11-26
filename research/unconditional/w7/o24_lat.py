"""O24 — push Delta(L=2) at beta=2.5 to >=5 sigma.

Four INDEPENDENT Metropolis chains, batched into one numpy process (one job, four chains as a
leading array axis).  Model: flat-coordinate density on the SU(2) fundamental domain

    V(x) = beta sum_P (1 - 1/2 ReTr U_P) + sum_l Phi(r_l),   Phi = -log J = -2 log|sin(r/2)/(r/2)|
    U_l = exp(i x_l.sigma/2),   |x_l| < 2*pi  (ball = group manifold exactly once)
    slice {<v,X> = 0} enforced EXACTLY.  mode A (v_l identical on all links) uses PAIR-TRANSFER
    proposals x_l += w, x_l' -= w, which conserve <v,X> with no global cross-talk, so the
    maintained plaquette traces stay exact.

usage: py -3.12 o24_lat.py <mode:A> <nsweep> <eps> [epsg-ignored] [seed0] [burn] [tag]
"""
import sys
import time
import numpy as np
import o11_physical as o

BETA = 2.5
NL = 64
DOM2 = (2 * np.pi) ** 2
tab, near = o.lat_tables(2)
NPQ = len(tab)
ROWS = np.array([near[l] for l in range(NL)], dtype=np.intp)
C0, C1, C2, C3 = (tab[:, k].astype(np.intp) for k in range(4))


def quat(x):
    r = np.sqrt((x * x).sum(-1))
    w = np.cos(0.5 * r)
    s = np.where(r > 1e-14, np.sin(0.5 * r) / np.where(r > 1e-14, r, 1.0), 0.5)
    return w, x * s[..., None]


def phi_r2(r2):
    r = np.sqrt(r2)
    q = np.where(r > 1e-14, np.sin(0.5 * r) / np.where(r > 1e-14, 0.5 * r, 1.0), 1.0)
    return -2.0 * np.log(np.abs(q))


def _g(w, a, idx):
    idx = np.broadcast_to(np.asarray(idx, np.intp), (w.shape[0],) + np.shape(idx)[-1:])
    r = np.arange(w.shape[0])[:, None]
    return w[r, idx], a[r, idx]


def ht(w, a, i0, i1, i2, i3):
    """(1/2)Re Tr(U0 U1 U2^dag U3^dag) from quaternions, fully vectorised."""
    w0, a0 = _g(w, a, i0)
    w1, a1 = _g(w, a, i1)
    w2, b2 = _g(w, a, i2)
    w3, b3 = _g(w, a, i3)
    a2, a3 = -b2, -b3
    tw = w0 * w1 - (a0 * a1).sum(-1)
    ta = w0[..., None] * a1 + w1[..., None] * a0 - np.cross(a0, a1)
    sw = w2 * w3 - (a2 * a3).sum(-1)
    sa = w2[..., None] * a3 + w3[..., None] * a2 - np.cross(a2, a3)
    return tw * sw - (ta * sa).sum(-1)


def full_V(x):
    w, a = quat(x)
    tr = ht(w, a, C0, C1, C2, C3)
    ph = phi_r2((x * x).sum(-1)).sum(-1)
    return BETA * (NPQ - tr.sum(-1)) + ph, tr, ph


def run(v, G, ns, eps, rng, h=1e-3, burn=800):
    assert np.allclose(v, v[0]), "pair-transfer assumes one common v per link"
    x = rng.normal(0, 0.55, (G, NL, 3))
    x -= np.einsum('gla,la->g', x, v)[:, None, None] * v[None]
    w, a = quat(x)
    tr = ht(w, a, C0, C1, C2, C3)
    phS = phi_r2((x * x).sum(-1)).sum(-1)
    ar = np.arange(G)
    nacc = np.zeros(G)
    ndom = 0
    Ks = np.empty((ns, G))
    Fs = np.empty((ns, G))
    for sw in range(ns + burn):
        for _ in range(NL):
            l = rng.integers(NL, size=G)
            l2 = (l + rng.integers(1, NL, size=G)) % NL
            wv = rng.uniform(-eps, eps, (G, 3))
            xn = x.copy()
            xn[ar, l] += wv
            xn[ar, l2] -= wv
            r2n = (xn * xn).sum(-1)
            phn = phi_r2(r2n).sum(-1)
            wn, an = quat(xn)
            rows = np.concatenate([ROWS[l], ROWS[l2]], axis=1)
            new = ht(wn, an, C0[rows], C1[rows], C2[rows], C3[rows])
            old = tr[ar[:, None], rows]
            dS = BETA * (old - new).sum(-1) + (phn - phS)
            m = (dS <= 0) | (rng.random(G) < np.exp(-np.maximum(dS, 0.0)))
            dom = (r2n < DOM2).all(-1)
            ndom += int((~dom).sum())
            m &= dom
            nacc += m
            if m.any():
                x[m] = xn[m]
                w[m] = wn[m]
                a[m] = an[m]
                phS[m] = phn[m]
                trn = tr.copy()
                trn[ar[:, None], rows] = new
                tr[m] = trn[m]
        if sw < burn:
            continue
        j = sw - burn
        V0 = BETA * (NPQ - tr.sum(-1)) + phS
        Vp = full_V(x + h * v[None])[0]
        Vm = full_V(x - h * v[None])[0]
        Ks[j] = (Vp - 2 * V0 + Vm) / h ** 2
        Fs[j] = (Vp - Vm) / (2 * h)
    return Ks, Fs, nacc / ((ns + burn) * NL), ndom / (G * (ns + burn) * NL), x


def acf(x, ml=600):
    y = x - x.mean()
    vv = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / vv for k in range(1, min(ml, len(y) - 2))])


def geyer_tau(x, ml=600):
    r = np.concatenate(([1.0], acf(x, ml)))
    if len(r) % 2:
        r = r[:-1]
    s = np.maximum(r[0::2] + r[1::2], 0.0)
    tot = 0.0
    for si in s:
        if si <= 0:
            break
        tot += si
    return max(2 * tot - 1.0, 1.0)


if __name__ == "__main__":
    mode = sys.argv[1]
    ns = int(sys.argv[2])
    eps = float(sys.argv[3])
    seed0 = int(sys.argv[5]) if len(sys.argv) > 5 else 400
    burn = int(sys.argv[6]) if len(sys.argv) > 6 else 800
    tag = sys.argv[7] if len(sys.argv) > 7 else 'run'
    G = 4
    v = np.zeros((NL, 3))
    if mode == 'A':
        v[:, 2] = 1.0 / np.sqrt(NL)
    else:
        v[0, 2] = 1.0
    seeds = [seed0 + i for i in range(G)]
    t0 = time.time()
    Ks, Fs, acc, domrej, x = run(v, G, ns, eps, np.random.default_rng(seeds), burn=burn)
    el = time.time() - t0
    np.savez(f"o24_{tag}.npz", K=Ks, F=Fs, seeds=np.array(seeds), mode=mode, ns=ns, eps=eps,
             acc=acc, domrej=domrej, el=el)
    r = np.linalg.norm(x, axis=-1)
    print(f"[O24] mode={mode} tag={tag} nsw={ns} burn={burn} eps={eps} G={G} seeds={seeds} "
          f"wall={el:.1f}s ({1e3*el/(ns+burn):.2f} ms/sweep for {G} chains)")
    print(f"[O24] acc_link={np.round(acc,3)} domain-reject={domrej:.4f} "
          f"final r: mean={r.mean():.3f} max={r.max():.3f} (2pi={2*np.pi:.3f})")
    for g in range(G):
        k, f = Ks[:, g], Fs[:, g]
        tk, tf = geyer_tau(k), geyer_tau(f)
        print(f"   g={g} seed={seeds[g]} E[kappa]={k.mean():+.4f} Var(F)={f.var():.4f} "
              f"E[F]={f.mean():+.4f} Delta={k.mean()-f.var():+.4f} tauK={tk:.1f} tauF={tf:.1f} "
              f"Neff={ns/max(tk,tf):.0f}")
