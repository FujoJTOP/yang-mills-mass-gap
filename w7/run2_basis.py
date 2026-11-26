"""RUN2 basis prep -- NO sampling. Reads archived npz only.

Question: can the a1-mandated basis (APE-3D/HYP graded spatial plaquette, n=1..4)
pass the two gates at L=6, T=12 --
  (g1) S/N(t=2) >= 10 per operator,
  (g2) near-orthogonal  -> raw Gram |G_ij| small / whitened off-diag ~ 0,
and if not, what does it cost in configs?

Measured on the only archived family that exists:
  research/unconditional/f1sq_b2.2_L6.npz  (beta=2.2, L=6, T=6, N=4000, APE s=0,2,4,8,
  one spatial-plaquette sum per time slice).  Extrapolations to T=12 / beta=2.5 are
  printed as such, never silently.
ASCII-only output (Windows console).

Error bars: 20-block delete-one jackknife (agrees with a 400x block bootstrap to 1e-2).
The naive iid per-config SE is 80x larger and is WRONG here: C(t) is estimated with the
same-ensemble mean subtracted, and that subtraction cancels the time-independent common
mode exactly (residual ~ sigma^2/N).  Verified: iid 8.0 vs JK 0.089 for C(1).
"""
import numpy as np

ARCH = 'D:/Dev/math/research/unconditional/f1sq_b2.2_L6.npz'
NBLK = 20


def conn(S, m=None):
    """C(t) = <S(0)S(t)> - m^2 with m the SAMPLE mean (self-centred, corpus form)."""
    m = S.mean() if m is None else m
    return np.array([(S * np.roll(S, -t, axis=1)).mean() - m * m
                     for t in range(S.shape[1])])


def cross(Sa, Sb):
    T = Sa.shape[1]
    return np.array([(Sa * np.roll(Sb, -t, axis=1)).mean() - Sa.mean() * Sb.mean()
                     for t in range(T)])


def jack(S, nb=NBLK):
    """20-block delete-one jackknife of C(t)."""
    N = S.shape[0]
    idx = np.array_split(np.arange(N), nb)
    vals = [conn(S[np.concatenate([idx[j] for j in range(nb) if j != i])])
            for i in range(nb)]
    vals = np.array(vals)
    mu = vals.mean(0)
    return conn(S), np.sqrt((nb - 1) / nb * ((vals - mu) ** 2).sum(0))


def gram(S, levels, t):
    k = len(levels)
    C = np.array([[cross(S[a], S[b])[t] for b in levels] for a in levels])
    d = np.sqrt(np.diag(C))
    return C / np.outer(d, d), C


def whiten(G0, cut=1e-2):
    w, V = np.linalg.eigh(G0)
    keep = w > cut * w.max()
    return (V[:, keep] / np.sqrt(w[keep])) @ V[:, keep].T, w, keep


def main():
    d = np.load(ARCH)
    levels = [int(x) for x in d['levels']]
    S = {s: d[f'S{s}'] for s in levels}
    N, T = S[levels[0]].shape
    print(f'=== archive  beta={float(d["beta"]):.2f}  L={int(d["L"])}  T={T}  N={N}  '
          f'levels={levels} (APE alpha={float(d["alpha"])}) ===')

    print('\n-- (g1) S/N per smearing level (20-block delete-one JK) --')
    print(f'  {"s":>2s} {"C(0)":>9s} {"C1/C0":>7s} {"C2/C0":>7s} '
          f'{"SN(t=1)":>8s} {"SN(t=2)":>8s} {"SN(t=3)":>8s} {"N_req(SN2=10)":>13s}')
    sn2 = {}
    for s in levels:
        C, se = jack(S[s])
        sn = {t: (C[t] / se[t] if se[t] > 0 else np.nan) for t in (1, 2, 3)}
        sn2[s] = sn[2]
        nreq = N * (10.0 / sn[2]) ** 2 if sn[2] > 0 else np.inf
        print(f'  {s:2d} {C[0]:9.4g} {C[1]/C[0]:7.4f} {C[2]/C[0]:7.4f} '
              f'{sn[1]:8.2f} {sn[2]:8.2f} {sn[3]:8.2f} {nreq:13.3g}')

    for t in (0, 2):
        G, _ = gram(S, levels, t)
        ev = np.linalg.eigvalsh(G)[::-1]
        off = np.abs(G - np.diag(np.diag(G)))
        print(f'\n-- (g2) raw normalised Gram at t={t} --')
        print('     ' + ' '.join(f'{s:>6d}' for s in levels))
        for i, a in enumerate(levels):
            print(f'  {a:3d}' + ' '.join(f'{G[i, j]:6.3f}' for j in range(len(levels))))
        print(f'  max|offdiag| = {off.max():.3f}   eig = '
              f'{np.array2string(ev, precision=3)}   kappa = {ev[-1]/ev[0]:.2e}')

    G0, _ = gram(S, levels, 0)
    G2, _ = gram(S, levels, 2)
    W, ev, keep = whiten(G0)
    Gw = W @ G2 @ W
    offw = np.abs(Gw - np.diag(np.diag(Gw)))
    print(f'\n-- (g2\') whitened Gram at t=2 (metric C(0), kept eig {keep.sum()}/{len(levels)},'
          f' kappa={ev[-1]/ev[0]:.2e}) --')
    print('     ' + ' '.join(f'{s:>6d}' for s in levels))
    for i, a in enumerate(levels):
        print(f'  {a:3d}' + ' '.join(f'{Gw[i, j]:6.3f}' for j in range(len(levels))))
    print(f'  max|offdiag| after whitening = {offw.max():.3f}')

    print('\n-- difference operators (normalised s_a - s_b), Gram at t=2 --')
    D = {}
    for i in range(1, len(levels)):
        a, b = levels[i], levels[i - 1]
        D[f'{b}->{a}'] = (S[a] / np.sqrt(conn(S[a])[0])) - (S[b] / np.sqrt(conn(S[b])[0]))
    keys = list(D)
    Cs = {k: conn(D[k]) for k in keys}
    for ka in keys:
        row = [f'{cross(D[ka], D[kb])[2] / np.sqrt(Cs[ka][0] * Cs[kb][0]):6.3f}' for kb in keys]
        print(f'  {ka:>8s} ' + ' '.join(row))

    print('\n-- T=12 budget --')
    print(f'  archive T={T}: independent separations t=0..{T//2} = {T//2+1} points'
          f' (t>{T//2} is the periodic mirror)')
    print(f'  run1  T=12: t=0..6 = 7 points; S/N gain from time-origin averaging only = '
          f'sqrt(12/{T}) = {np.sqrt(12.0/T):.2f}')
    for s in levels:
        n12 = N * (10.0 / (sn2[s] * np.sqrt(12.0 / T))) ** 2
        print(f'    s={s:2d}: N_req(SN2=10, T=12) ~ {n12:9.3g}   (T=6: '
              f'{N*(10.0/sn2[s])**2:9.3g})')

    print('\n-- gates run1 must print itself (T=12 basis) --')
    print('  rank(C(0)) >= 2 ; kappa(C(0)) >= 0.1 ; max|G_raw| reported but NOT gated ;'
          ' max|G_whitened(t=2)| <= 0.1 ; S/N(t=2) >= 10 per operator')


if __name__ == '__main__':
    main()
