"""t23_monopole.py — T2/T3 witness: 4D defect-gas numbers from EXISTING configs.

Purpose (task T23-TEMPLATE, T3): decide whether the 4D deficiency of the
defect-gas class (C2) is "not dilute enough" (quantitative) or "no gas"
(qualitative).  We measure, in the SAME maximal-centre gauge (MCG) used for the
corpus vortex numbers (probe-center/REPORT.md):

  * Z(3) vortex density per dual plaquette, and the cluster count
    (percolation indicator);
  * abelian-projected monopole density per lattice 3-cube (dual site), for the
    three U(1) subgroups (colour pairs 01/02/12) of the SU(3) link;
  * the strong-coupling expansion parameter u(beta)=I2(beta)/I1(beta) and the
    naive Kotecky-Preiss-type sum 6u (6 plaquettes/site).

Densities are reported in lattice units AND, at beta=6.0, in physical units
(a^-1 = 2.0 GeV, i.e. a = 0.10 fm -- standard SU(3) Wilson, [R]); the
dimensionless diluteness n*xi^3 uses xi = 1/m_{0++} = 1/1.710 GeV = 0.116 fm.

Everything is a gauge-fixed DIAGNOSTIC, not a certificate (Singer/Elitzur:
no gauge-invariant monopole density).  numpy only (+numba via mcg); one
process; L=6; 20 snapshots per beta.

Run:  py -3.12 -I t23_monopole.py
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/probe-center')
import mcg  # noqa: E402

RUNS = 'D:/Dev/math/research/lattice-cpb/runs'
FILES = {2.6: 'L6_b2.600_s51.npz', 3.0: 'L6_b3.000_s53.npz', 6.0: 'L6_b6.000_s32.npz'}
NSNAP = 20
NSW, MAXIT = 30, 6


def random_su3(L, rng):
    g = np.zeros((L, L, L, L, 3, 3), np.complex128)
    a = rng.normal(size=(L, L, L, L, 3, 3)) + 1j * rng.normal(size=(L, L, L, L, 3, 3))
    q, r = np.linalg.qr(a)
    d = np.einsum('...ii->...i', r)
    d = d / np.abs(d)
    q = q * d[..., None, :]
    det = np.linalg.det(q)
    q = q / (det ** (1.0 / 3.0))[..., None, None]
    g[...] = q
    return g


def gauge_fix(U, L, rng):
    best = None
    for kind in ('id', 'rand'):
        if kind == 'id':
            g = np.zeros((L, L, L, L, 3, 3), np.complex128)
            for i in range(3):
                g[..., i, i] = 1.0
        else:
            g = random_su3(L, rng)
        for _ in range(NSW):
            mcg.mcg_sweep(U, g, L, MAXIT)
        V = mcg.apply_gauge(U, g, L)
        R = mcg.center_functional(V, L)
        if best is None or R > best[0]:
            best = (R, V)
    return best


def plaq_flux_abelian(theta, L):
    """f[mu,nu](x) = theta_mu(x)+theta_nu(x+mu)-theta_mu(x+nu)-theta_nu(x),
    reduced to (-pi,pi]; array (4,4,L,L,L,L), only mu<nu used."""
    f = np.zeros((4, 4, L, L, L, L))
    for mu in range(4):
        for nu in range(mu + 1, 4):
            raw = (theta[mu]
                   + np.roll(theta[nu], -1, axis=mu)
                   - np.roll(theta[mu], -1, axis=nu)
                   - theta[nu])
            f[mu, nu] = np.angle(np.exp(1j * raw))
    return f


def monopole_cubes(f, L):
    """Q(x; a<b<c) = round( outward flux / 2pi ) for the 4 triples per site."""
    tot_abs = 0.0
    n_nonzero = 0
    n_cube = 0
    for a, b, c in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        s = (f[a, b]
             - np.roll(f[a, b], -1, axis=c)
             + f[b, c]
             - np.roll(f[b, c], -1, axis=a)
             + f[a, c]
             - np.roll(f[a, c], -1, axis=b))
        Q = np.rint(s / (2.0 * np.pi))
        tot_abs += np.abs(Q).sum()
        n_nonzero += int((Q != 0).sum())
        n_cube += Q.size
    return tot_abs / n_cube, n_nonzero / n_cube


def main():
    t0 = time.time()
    rng = np.random.default_rng(20261125)
    out = {}
    for beta, fn in FILES.items():
        d = np.load(os.path.join(RUNS, fn), allow_pickle=False)
        snaps = d['snapshots']
        L = int(d['L'])
        ns = min(NSNAP, snaps.shape[0])
        rows = []
        for i in range(ns):
            U = snaps[i].copy()
            R, V = gauge_fix(U, L, rng)
            # --- Z(3) vortex density (same pipeline as probe-center) ---
            n, _amb = mcg.center_project(V, L)
            flux = mcg.plaquette_flux(n, L)
            plaqs, labels, sizes = mcg.vortex_clusters(flux, L)
            dens = len(plaqs) / (6.0 * L ** 4)
            ncl = len(sizes)
            smax = float(np.max(sizes)) / (6.0 * L ** 4) if len(sizes) else 0.0
            # --- abelian-projected monopoles, three U(1) subgroups ---
            mon = []
            for (i_, j_) in ((0, 1), (0, 2), (1, 2)):
                theta = np.angle(V[..., i_, i_] * np.conj(V[..., j_, j_]))
                f = plaq_flux_abelian(theta, L)
                m_abs, m_frac = monopole_cubes(f, L)
                mon.append((m_abs, m_frac))
            mon = np.array(mon)  # (3, 2)
            rows.append(dict(R=R, dens=dens, ncl=ncl, smax=smax,
                             mon_abs=mon[:, 0].mean(), mon_frac=mon[:, 1].mean(),
                             mon_abs_by=mon[:, 0].tolist(),
                             mon_frac_by=mon[:, 1].tolist()))
            if (i + 1) % 5 == 0:
                print(f"  beta={beta} [{i+1}/{ns}] R={R:.3f} dens={dens:.3f} "
                      f"ncl={ncl} |mono|={mon[:,0].mean():.4f} "
                      f"frac={mon[:,1].mean():.4f} t={time.time()-t0:.0f}s",
                      flush=True)
        arr = {k: np.array([r[k] for r in rows])
               for k in ('R', 'dens', 'ncl', 'smax', 'mon_abs', 'mon_frac')}
        out[f'{beta}'] = dict(
            beta=beta, L=L, nsnap=ns,
            R=float(arr['R'].mean()), R_sem=float(arr['R'].std(ddof=1) / np.sqrt(ns)),
            vdens=float(arr['dens'].mean()), vdens_sem=float(arr['dens'].std(ddof=1) / np.sqrt(ns)),
            nclust=float(arr['ncl'].mean()),
            vmaxfrac=float(arr['smax'].mean()),
            mon_abs=float(arr['mon_abs'].mean()),
            mon_abs_sem=float(arr['mon_abs'].std(ddof=1) / np.sqrt(ns)),
            mon_frac=float(arr['mon_frac'].mean()),
            mon_frac_sem=float(arr['mon_frac'].std(ddof=1) / np.sqrt(ns)),
        )
        o = out[f'{beta}']
        print(f"[beta={beta}] vdens={o['vdens']:.4f}+-{o['vdens_sem']:.4f} "
              f"nclust={o['nclust']:.2f} spanfrac={o['vmaxfrac']:.3f} | "
              f"mono/abs={o['mon_abs']:.4f}+-{o['mon_abs_sem']:.4f} "
              f"mono/frac={o['mon_frac']:.4f}+-{o['mon_frac_sem']:.4f}", flush=True)

    # --- controls: (a) cold links U=exp(i eps X), (b) random links ---
    # calibrates the UV sensitivity of the gauge-fixed abelian/vortex projection.
    Lc = 6

    def measure(U):
        R, V = gauge_fix(U, Lc, rng)
        n, _ = mcg.center_project(V, Lc)
        flux = mcg.plaquette_flux(n, Lc)
        plaqs, labels, sizes = mcg.vortex_clusters(flux, Lc)
        mm = []
        for (i_, j_) in ((0, 1), (0, 2), (1, 2)):
            theta = np.angle(V[..., i_, i_] * np.conj(V[..., j_, j_]))
            f = plaq_flux_abelian(theta, Lc)
            mm.append(monopole_cubes(f, Lc))
        return len(plaqs) / (6.0 * Lc ** 4), float(np.mean([x[0] for x in mm])), R

    ctrl = {}
    for eps in (0.05, 0.5):
        vs, ms, Rs = [], [], []
        for _ in range(4):
            X = rng.normal(size=(4, Lc, Lc, Lc, Lc, 3, 3)) \
                + 1j * rng.normal(size=(4, Lc, Lc, Lc, Lc, 3, 3))
            X = 0.5 * (X + np.conj(np.swapaxes(X, -1, -2)))
            X = X - np.trace(X, axis1=-2, axis2=-1)[..., None, None] / 3.0
            w, Q = np.linalg.eigh(X)
            U = np.einsum('...ik,...k,...jk->...ij', Q, np.exp(1j * eps * w), np.conj(Q))
            det = np.linalg.det(U)
            U = U / (det ** (1.0 / 3.0))[..., None, None]
            v, m, R = measure(U)
            vs.append(v)
            ms.append(m)
            Rs.append(R)
        ctrl[f'cold_eps{eps}'] = dict(vdens=float(np.mean(vs)), mon_abs=float(np.mean(ms)),
                                     R=float(np.mean(Rs)))
        print(f"[cold eps={eps}] vdens={np.mean(vs):.4f} mono/abs={np.mean(ms):.4f} "
              f"R={np.mean(Rs):.2f}", flush=True)
    vs, ms, Rs = [], [], []
    for _ in range(5):
        U = np.stack([random_su3(Lc, rng) for _ in range(4)], axis=0)
        v, m, R = measure(U)
        vs.append(v)
        ms.append(m)
        Rs.append(R)
    ctrl['random'] = dict(vdens=float(np.mean(vs)), mon_abs=float(np.mean(ms)),
                          R=float(np.mean(Rs)))
    print(f"[random] vdens={np.mean(vs):.4f} mono/abs={np.mean(ms):.4f} "
          f"(UV noise floor)", flush=True)
    out['controls'] = ctrl

    # --- strong-coupling parameter u(beta)=I2/I1 (SU(2) character ratio) ---
    from scipy.special import iv
    ub = {}
    for b in (1.0 / 48.0, 0.5, 2.5, 6.0):
        ub[f'{b:.6f}'] = float(iv(2, b) / iv(1, b))
    out['u'] = ub
    print("[u] u(beta)=I2/I1: " + ", ".join(f"{k}:{v:.5f}" for k, v in ub.items()))

    # --- physical units at beta=6.0 (a^-1 = 2.0 GeV) and diluteness ---
    hbar_c = 0.1973  # GeV*fm
    ainv = 2.0  # GeV  (beta=6.0 SU(3) Wilson, [R])
    a_fm = hbar_c / ainv
    xi = hbar_c / 1.710  # fm, = 1/m_0++
    o6 = out['6.0']
    # monopole number density per 3-cube -> per fm^3 (cube volume = a^3)
    n_mon_fm3 = o6['mon_abs'] / a_fm ** 3
    n_dil = n_mon_fm3 * xi ** 3
    # vortex sheet area per 4-volume: (6*vdens) plaquettes per site, each area a^2
    area_dens_fm2 = 6.0 * o6['vdens'] / a_fm ** 2
    area_dil = area_dens_fm2 * xi ** 2
    out['phys'] = dict(ainv_GeV=ainv, a_fm=a_fm, xi_fm=xi,
                       n_mon_fm3=n_mon_fm3, n_mon_xi3=n_dil,
                       area_dens_fm2=area_dens_fm2, area_xi2=area_dil,
                       sigma_fm2=0.497, sqrt_sigma_GeV=0.440, m0_GeV=1.710,
                       Lam_MSbar_GeV=0.250)
    print("[phys] a={a_fm:.4f} fm  xi={xi_fm:.4f} fm".format(**out['phys']))
    print(f"[phys] n_mon = {n_mon_fm3:.3f} /fm^3 ; n_mon*xi^3 = {n_dil:.4f}")
    print(f"[phys] vortex area density = {area_dens_fm2:.2f} /fm^2 ; "
          f"area*xi^2 = {area_dil:.4f}")
    print(f"[phys] sigma = "f"{out['phys']['sigma_fm2']:.3f} /fm^2 ; "
          f"sqrt(sigma)/Lambda = {0.440/0.250:.3f} ; sigma/Lambda^2 = {(0.440/0.250)**2:.3f}")

    with open('t23_monopole.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f"saved t23_monopole.json  ({time.time()-t0:.0f}s)")


if __name__ == '__main__':
    main()
