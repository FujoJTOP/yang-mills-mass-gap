"""
SU(2) beta=2.5 correlation length from the connected 0++-like plaquette
correlator C(r) along one axis, raw and APE-smeared.  Used only to state the
xi<L regime for the softest mode k=2*pi/L (estimate, not a proof).
"""
import sys
import argparse
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W  # noqa: E402

adj = W.adj


def plaq_field(U):
    """site field P(x) = mean over the 6 plane types of 1/2 Re tr U_p."""
    L = U.shape[1]
    tot = np.zeros((L, L, L, L))
    n = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))
            tot += 0.5 * np.real(np.trace(Up, axis1=-2, axis2=-1))
            n += 1
    return tot / n


def corr(P, axis=0):
    L = P.shape[axis]
    X = P - P.mean()
    return np.array([np.mean(X * np.roll(X, -r, axis=axis)) for r in range(L)])


def fit_m(C, L, tmin=1, tmax=None):
    if tmax is None:
        tmax = L // 2
    ts = np.arange(tmin, tmax + 1)
    y = C[ts]
    ms = np.linspace(0.05, 4.0, 6000)
    Z = np.exp(-ms[:, None] * ts[None, :]) + np.exp(-ms[:, None] * (L - ts)[None, :])
    A = (Z @ y) / np.sum(Z * Z, axis=1)
    r = np.sum((y[None, :] - A[:, None] * Z) ** 2, axis=1)
    k = int(np.argmin(r))
    if k == 0 or k == len(ms) - 1:
        return float('nan')
    return float(ms[k])


def run(L, beta, N, therm, gap, seed, nsmear=(0, 8), alpha=0.5):
    rng = np.random.default_rng(seed)
    cfgs = W.sample2(L, beta, rng, therm, N, gap)
    levels = list(nsmear)
    C = {n: np.zeros(L) for n in levels}
    plaq = []
    for U in cfgs:
        plaq.append(float(plaq_field(U).mean()))
        prev = 0
        Us = U
        for n in levels:
            if n > prev:
                Us = W.ape_all2(Us, alpha=alpha, nsteps=n - prev)
                prev = n
            C[n] += corr(plaq_field(Us), axis=0)
    nc = len(cfgs)
    out = {}
    print(f"# SU2 L={L} beta={beta} N={nc} plaq={np.mean(plaq):.4f}")
    for n in levels:
        Cc = C[n] / nc
        m = fit_m(Cc, L)
        m1 = None
        if Cc[1] != 0 and Cc[0] != 0 and Cc[2] != 0 and min(Cc[0], Cc[1], Cc[2]) > 0:
            m1 = float(np.arccosh((Cc[0] + Cc[2]) / (2 * Cc[1])))
        print(f"  n_smear={n}: C/C0={np.round(Cc/Cc[0],4)}  m_fit={m:.3f} "
              f"xi={1.0/m:.2f}  m_eff(1)={m1 if m1 is None else round(m1,3)}")
        out[n] = dict(m=m, C=Cc)
    np.savez(f'corr_su2_L{L}_b{beta}.npz', L=L, beta=beta, plaq=np.array(plaq),
             **{f'C{n}': out[n]['C'] for n in levels},
             **{f'm{n}': out[n]['m'] for n in levels})
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=8)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--N', type=int, default=400)
    ap.add_argument('--therm', type=int, default=400)
    ap.add_argument('--gap', type=int, default=4)
    ap.add_argument('--seed', type=int, default=321)
    a = ap.parse_args()
    run(a.L, a.beta, a.N, a.therm, a.gap, a.seed)
