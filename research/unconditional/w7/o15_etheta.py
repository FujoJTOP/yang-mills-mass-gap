"""
O15 — E[H_vv|0](N,theta) on ONE plaquette: exact structure, zero set, sign at theta=0.

Model = betastar_scaling.py (same Cartan-transverse slice, same J, same domain), but
B(theta) = exp(i*theta*T3), so N=2, theta=pi reproduces H1A-EXACT.md §2.1's B = i*sigma_3 exactly.
  X = sum_{a off-Cartan} xi_a T_a,  U = exp(iX),  S = beta(1 - ReTr(U B)/N),  measure e^{-S} J dX,
  J = prod_{i<j}[sin(mu_ij/2)/(mu_ij/2)]^2,  domain |mu_i-mu_j|<=2pi,  v = T3 = diag(1/2,-1/2,0,..)
  H_vv = d^2_t S(X+tv)|_{t=0} = -(beta/N) d^2_t ReTr(exp(i(X+tv)) B)

Exact identity (O15-THETAPI.md §1, verified to machine precision here):
  B(theta) = c P + (1-P) + i s V,  c=cos(theta/2), s=sin(theta/2), P=diag(1,1,0..), V=2T3
  H_vv(X;theta) = (beta/N)[ A0 + (1-c) A1 + s A2 ]
    A0 = -d^2_t ReTr e^{i(X+tv)},  A1 = d^2_t ReTr(e^{i(X+tv)}P),  A2 = d^2_t ImTr(e^{i(X+tv)}V)
  E[H|0] = (beta/N)[ <A0> + (1-c)<A1> + s<A2> ],  <> = slice average, weight e^{-S(theta)}J
  X -> DXD (D=diag(1,-1,1..)) keeps slice/domain/J and flips V  ==>  <A2> = 0, E[H] even in theta.

numpy only. Finite differences in t, h=1e-3 (O(h^2)); self-normalised weighted means, block errors.
"""
import numpy as np

from betastar_scaling import su_gens, centre_B

np.set_printoptions(precision=8, suppress=True)
Hd = 1e-3
IU = {}


def slices(N):
    G = su_gens(N)
    off = G[:N * (N - 1)]
    Tv = G[N * (N - 1)]
    P = np.diag(np.array([1.0, 1.0] + [0.0] * (N - 2), complex))
    IU[N] = np.triu_indices(N, 1)
    return off, Tv, P, 2 * Tv


def build(N, nwant, seed=11, chunk=50000, A=2 * np.pi):
    """beta-free pool on the slice; per-sample (J, qI, qP, qV, A0, A1, A2)."""
    rng = np.random.default_rng(seed)
    off, Tv, P, V = slices(N)
    iu = IU[N]
    out = {k: [] for k in ('J', 'qI', 'qP', 'qV', 'A0', 'A1', 'A2')}
    got = tried = 0
    while got < nwant:
        xi = rng.uniform(-A, A, size=(chunk, N * (N - 1)))
        tried += chunk
        X = np.zeros((chunk, N, N), complex)
        for a, T in enumerate(off):
            X += xi[:, a, None, None] * T[None]
        mu = np.linalg.eigvalsh(X)
        d = (mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]]
        keep = np.all(np.abs(d) <= 2 * np.pi, axis=1)
        if not keep.any():
            continue
        X, d = X[keep], d[keep]
        hf = d / 2
        J = np.prod((np.sin(hf) / hf) ** 2, axis=1)
        tr = {}
        for tag, sg in ((-1, -1.0), (0, 0.0), (+1, +1.0)):
            mu2, V2 = np.linalg.eigh(X + sg * Hd * Tv)
            U = (V2 * np.exp(1j * mu2)[:, None, :]) @ V2.conj().transpose(0, 2, 1)
            tr[tag] = (np.real(np.trace(U, axis1=1, axis2=2)),
                       np.real(np.trace(U @ P, axis1=1, axis2=2)),
                       np.imag(np.trace(U @ V, axis1=1, axis2=2)))
        A0 = -(tr[+1][0] - 2 * tr[0][0] + tr[-1][0]) / Hd ** 2
        A1 = (tr[+1][1] - 2 * tr[0][1] + tr[-1][1]) / Hd ** 2
        A2 = (tr[+1][2] - 2 * tr[0][2] + tr[-1][2]) / Hd ** 2
        m = min(nwant - got, len(J))
        for k, v in (('J', J), ('qI', tr[0][0]), ('qP', tr[0][1]), ('qV', tr[0][2]),
                     ('A0', A0), ('A1', A1), ('A2', A2)):
            out[k].append(v[:m])
        got += m
    r = {k: np.concatenate(v) for k, v in out.items()}
    r['acc'] = got / tried
    return r


def wmean(y, w, nblk=10):
    idx = np.array_split(np.arange(len(y)), nblk)
    m = np.array([np.sum(y[i] * w[i]) / np.sum(w[i]) for i in idx])
    return np.sum(y * w) / np.sum(w), m.std(ddof=1) / np.sqrt(len(m))


def stats(pool, N, theta, beta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    qB = c * pool['qP'] + (pool['qI'] - pool['qP']) - s * pool['qV']
    w = np.exp(-beta * (1 - qB / N)) * pool['J']
    a0, _ = wmean(pool['A0'], w)
    a1, _ = wmean(pool['A1'], w)
    a2, _ = wmean(pool['A2'], w)
    eh, dh = wmean((pool['A0'] + (1 - c) * pool['A1'] + s * pool['A2']) * beta / N, w)
    return dict(A0=a0, A1=a1, A2=a2, EH=eh, dH=dh)


def identity_check():
    print("[V1] pointwise identity  H = (beta/N)[A0 + (1-c)A1 + s A2]  (200 random slice pts, 4 thetas)")
    rng = np.random.default_rng(5)
    for N in (2, 3, 4):
        off, Tv, P, V = slices(N)
        xi = rng.uniform(-1, 1, size=(200, len(off)))
        X = np.zeros((200, N, N), complex)
        for a, T in enumerate(off):
            X += xi[:, a, None, None] * T[None]
        worst = 0.0
        for th in (0.0, np.pi / 3, np.pi, 4 * np.pi / 3):
            c, s = np.cos(th / 2), np.sin(th / 2)
            B = c * P + (np.eye(N) - P) + 1j * s * V
            q = {}
            for tag, sg in ((-1, -1.0), (0, 0.0), (+1, +1.0)):
                mu2, V2 = np.linalg.eigh(X + sg * Hd * Tv)
                U = (V2 * np.exp(1j * mu2)[:, None, :]) @ V2.conj().transpose(0, 2, 1)
                q[tag] = (np.real(np.einsum('nij,ji->n', U, B)),
                          np.real(np.trace(U, axis1=1, axis2=2)),
                          np.real(np.trace(U @ P, axis1=1, axis2=2)),
                          np.imag(np.trace(U @ V, axis1=1, axis2=2)))
            Hx = -(q[+1][0] - 2 * q[0][0] + q[-1][0]) / Hd ** 2
            A0 = -(q[+1][1] - 2 * q[0][1] + q[-1][1]) / Hd ** 2
            A1 = (q[+1][2] - 2 * q[0][2] + q[-1][2]) / Hd ** 2
            A2 = (q[+1][3] - 2 * q[0][3] + q[-1][3]) / Hd ** 2
            form = A0 + (1 - c) * A1 + s * A2
            worst = max(worst, np.max(np.abs(Hx - form)) / np.max(np.abs(Hx)))
        print(f"   N={N}: max rel |H_exact - formula| = {worst:.2e}")


def main():
    print("=" * 80, flush=True)
    identity_check()
    print("=" * 80)
    print("[V2] pools (Cartan-transverse slice, de Wit-'t Hooft J, domain |mu_ij|<=2pi)", flush=True)
    pools, n = {}, {2: 200000, 3: 200000, 4: 15000}
    for N in (2, 3, 4):
        pools[N] = build(N, n[N])
        a2, d2 = wmean(pools[N]['A2'], pools[N]['J'])
        print(f"   N={N}: {n[N]} valid, acc {pools[N]['acc']:.1e}, <A2>_Haar = {a2:+.2e} +- {d2:.1e}", flush=True)

    print("=" * 80)
    print("[V3] H1A cross-check N=2 (targets: E[H|0](theta=0,beta=2.5)=0.553061 ; E[H|0](theta=pi)=0 exact)")
    for beta in (0.5, 2.5, 6.0):
        r0, rp = stats(pools[2], 2, 0.0, beta), stats(pools[2], 2, np.pi, beta)
        print(f"   beta={beta:5.2f}: E[H](theta=0)={r0['EH']:+.6f}+-{r0['dH']:.1e}  "
              f"E[H](theta=pi)={rp['EH']:+.2e}+-{rp['dH']:.1e}  <A0>={r0['A0']:+.6f} <A1>={r0['A1']:+.6f}", flush=True)

    print("=" * 80)
    print("[V4] E[H|0](N,theta), family B=exp(i theta T3); theta=pi is H1A's i*sigma_3 for N=2")
    ths = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi, 5 * np.pi / 4, 3 * np.pi / 2, 7 * np.pi / 4, 2 * np.pi]
    print("   theta/pi:  " + "".join(f"{t / np.pi:8.2f}" for t in ths))
    for N in (2, 3, 4):
        for beta in (2.5, 6.0):
            row = [stats(pools[N], N, t, beta) for t in ths]
            print(f"   N={N} b={beta:3.1f}: " + "".join(f"{r['EH']:+8.4f}" for r in row)
                  + f"    (min {min(r['EH'] for r in row):+.4f})", flush=True)

    print("=" * 80)
    print("[V5] theta=0: E[H|0]=(beta/N)<A0>, <A0> carries the whole sign; pointwise A0 stats")
    for N in (2, 3, 4):
        p = pools[N]
        r = stats(p, N, 0.0, 2.5)
        print(f"   N={N}: <A0>_Haar={wmean(p['A0'], p['J'])[0]:+.5f} <A0>_(beta=2.5)={r['A0']:+.5f} "
              f"E[H](beta=2.5)={r['EH']:+.5f}+-{r['dH']:.1e} | pointwise min A0={p['A0'].min():+.4f}, "
              f"frac(A0<0)={np.mean(p['A0'] < 0):.3f}", flush=True)
        for beta in (0.5, 5.0, 10.0):
            print(f"        beta={beta:4.1f}: <A0>={stats(p, N, 0.0, beta)['A0']:+.5f} "
                  f"E[H]={stats(p, N, beta=beta, theta=0.0)['EH']:+.5f}", flush=True)

    print("=" * 80)
    print("[V6] Z_N centre element (BETASTAR-N §4, a DIFFERENT B): N=3, beta=6 target E[H]=+0.306")
    for N, nw in ((2, 40000), (3, 40000)):
        B = centre_B(N)
        off, Tv, _, _ = slices(N)
        iu = IU[N]
        rng = np.random.default_rng(3)
        hs, ws, got = [], [], 0
        while got < nw:
            xi = rng.uniform(-2 * np.pi, 2 * np.pi, size=(40000, N * (N - 1)))
            X = np.zeros((40000, N, N), complex)
            for a, T in enumerate(off):
                X += xi[:, a, None, None] * T[None]
            mu = np.linalg.eigvalsh(X)
            d = (mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]]
            keep = np.all(np.abs(d) <= 2 * np.pi, axis=1)
            if not keep.any():
                continue
            X, d = X[keep], d[keep]
            hf = d / 2
            J = np.prod((np.sin(hf) / hf) ** 2, axis=1)

            def fB(Xs):
                mu2, V2 = np.linalg.eigh(Xs)
                U = (V2 * np.exp(1j * mu2)[:, None, :]) @ V2.conj().transpose(0, 2, 1)
                return np.real(np.trace(U @ B, axis1=1, axis2=2)) / N

            W = np.exp(-6.0 * (1 - fB(X))) * J
            hh = -6.0 / N * (fB(X + Hd * Tv) - 2 * fB(X) + fB(X - Hd * Tv)) / Hd ** 2
            m = min(nw - got, len(J))
            hs.append(hh[:m]); ws.append(W[:m]); got += m
        hh, W = np.concatenate(hs), np.concatenate(ws)
        e, d = wmean(hh, W)
        print(f"   N={N}: Z_N centre, beta=6 -> E[H|0] = {e:+.4f} +- {d:.4f}", flush=True)


if __name__ == '__main__':
    main()
