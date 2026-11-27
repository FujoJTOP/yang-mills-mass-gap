#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ht3_classang.py -- HT3: measure Var(2 arccos a0), a0 = (1/2) tr U_l, on the
beta=2.5 SU(2) flat-chart configurations already on disk.

ZERO new MC.  Inputs (read-only):
  research/v3/w7/gevp_4d_L{4,6,8}_field.npy   SU(2) U (4,L,L,L,L,2,2) complex128,
                                              saved by gevp_4d.py: BETA=2.5, ALPHA=0.5,
                                              make_stream from m1_4d_run.py (SU(2) heat bath),
                                              therm=800, gap=10, cold, seeds 101/202/303,
                                              no gauge fixing; observe() smears a COPY
                                              (w4_blocked.ape_all2 starts with U.copy()).
  Everything else is deterministic quadrature / closed form.

Sections
  [A] exact Haar reference for the class-angle functional (quadrature + closed form)
  [B] what the Elitzur-allowed 1D closures would give for the SAME functional
  [C] measured class-angle statistics on the three stored beta=2.5 configs
  [D] histogram / CDF test against Haar
  [E] verdict numbers
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
W7 = os.path.join(ROOT, 'research', 'v3', 'w7')
OUT = []


def P(s=''):
    print(s, flush=True)
    OUT.append(str(s))


def hdr(t):
    P()
    P('=' * 78)
    P(t)
    P('=' * 78)


# ---------------------------------------------------------------- quadrature
def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * w


PHI, WPHI = gl(4000, 0.0, 2.0 * np.pi)
HAARW = np.sin(0.5 * PHI) ** 2 * WPHI                # sin^2(phi/2) dphi  (Haar class measure)


def stats(phi, w):
    w = w / w.sum()
    a0 = np.cos(0.5 * phi)
    m1 = (w * phi).sum()
    m2 = (w * phi ** 2).sum()
    am = (w * a0).sum()
    av = (w * a0 ** 2).sum() - am ** 2
    a4 = (w * a0 ** 4).sum()
    return dict(Ephi=m1, Ephi2=m2, Var=m2 - m1 ** 2, third=m2 / 3.0,
                Ea0=am, Vara0=av, Ea04=a4,
                big=(w * (np.abs(a0) > 0.9)).sum())


hdr('[A] exact Haar reference, a0 = (1/2)tr U = cos(phi/2), phi = 2 arccos a0 in [0,2pi)')
H = stats(PHI, HAARW)
P(f'  quadrature 4000 pt GL:  E[phi]   = {H["Ephi"]:.9f}   (closed form pi = {np.pi:.9f})')
P(f'                         E[phi^2] = {H["Ephi2"]:.9f}   (4pi^2/3 - 2 = {4*np.pi**2/3-2:.9f})')
P(f'                         Var(phi) = {H["Var"]:.9f}   (pi^2/3 - 2   = {np.pi**2/3-2:.9f})')
P(f'                         (1/3)E[phi^2] = {H["third"]:.9f}   ((4pi^2-6)/9 = {(4*np.pi**2-6)/9:.9f})')
th = np.arccos(0.9)
pbig = 2.0 * (th / 2.0 - np.sin(2 * th) / 4.0) / (np.pi / 2.0)
P(f'                         E[a0] = {H["Ea0"]:+.2e}, Var(a0) = {H["Vara0"]:.9f} (=1/4), '
  f'E[a0^4] = {H["Ea04"]:.9f} (=1/8)')
P(f'                         P(|a0|>0.9) = {H["big"]:.6f}   (closed form {pbig:.6f})')
P('  => THE IDENTITY: 3.719824 = (1/3)E[(2 arccos a0)^2] = (1/3)E[|xi|^2], NOT Var(2 arccos a0).')
P('     Var(2 arccos a0)|Haar = pi^2/3 - 2 = 1.289868134, because E[2 arccos a0] = pi.')
P('     So the two differ by the mean square of the mean: 1.289868 + pi^2 = 11.159473 = 3 x 3.719824.')

hdr('[B] same functional under the two Elitzur-allowed 1D closures (what a NON-Haar edge would give)')
P('  w(phi) ~ sin^2(phi/2) cosh(kappa cos(phi/2))   [even in a0 -> allowed by center symmetry]')
for kap in (2.5, 3.927, 6.0, 12.0):
    w = HAARW * np.cosh(kap * np.cos(0.5 * PHI))
    s = stats(PHI, w)
    P(f'   kappa={kap:6.3f}:  Var(phi) = {s["Var"]:8.4f}   (1/3)E[phi^2] = {s["third"]:8.4f}   '
      f'Var(a0) = {s["Vara0"]:7.4f}   P(|a0|>0.9) = {s["big"]:.4f}')
P('  => a mass-like bias pushes phi toward 0: Var(phi) DROPS below 1.2899 at fixed mean-square.')
P('     Reaching Var(phi) ~ 4.0-6.6 needs mass near BOTH phi~0 and phi~2pi (a0 ~ +1 and -1):')
for kap in (4.0, 8.0):
    for kap2 in (4.0, 12.0):
        w = HAARW * (np.cosh(kap * np.cos(0.5 * PHI)) + np.cosh(kap2 * np.cos(0.5 * (PHI - 2 * np.pi))))
        s = stats(PHI, w)
        P(f'   bimodal kappa={kap}/{kap2}: Var(phi) = {s["Var"]:8.4f}  (1/3)E[phi^2] = {s["third"]:8.4f}')

hdr('[C] measured: beta=2.5 SU(2) flat-chart fields on disk (zero new MC)')


def plaquette(U):
    tot, n = 0.0, 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Pm = U[mu] @ np.roll(U[nu], -1, axis=mu) @ np.conj(np.swapaxes(np.roll(U[mu], -1, axis=nu), -1, -2)) \
                 @ np.conj(np.swapaxes(U[nu], -1, -2))
            tot += float(np.mean(0.5 * np.real(np.trace(Pm, axis1=-2, axis2=-1))))
            n += 1
    return tot / n


rows = []
pool_a0 = []
for L in (4, 6, 8):
    f = os.path.join(W7, f'gevp_4d_L{L}_field.npy')
    if not os.path.exists(f):
        P(f'  TOOL-FAIL: missing {f}')
        continue
    U = np.load(f)
    unit = float(np.abs(U @ np.conj(np.swapaxes(U, -1, -2)) - np.eye(2)).max())
    a0 = 0.5 * np.real(np.trace(U, axis1=-2, axis2=-1)).ravel()
    aim = float(np.abs(np.imag(np.trace(U, axis1=-2, axis2=-1))).max())
    phi = 2.0 * np.arccos(np.clip(a0, -1.0, 1.0))
    pool_a0.append(a0)
    # jackknife over x0 time slices (correlation-aware error on Var(phi))
    n0 = U.shape[1]
    vs = []
    for t in range(n0):
        m = np.ones(n0, bool)
        m[t] = False
        ph = phi.reshape(4, n0, -1)[:, m, :].ravel()
        vs.append(ph.var())
    vs = np.array(vs)
    jk = float(np.sqrt((n0 - 1) / n0 * ((vs - vs.mean()) ** 2).sum()))
    rows.append(dict(L=L, n=a0.size, plaq=plaquette(U), unit=unit, aim=aim,
                     Ea0=a0.mean(), Va0=a0.var(), Ea04=np.mean(a0 ** 4),
                     big=np.mean(np.abs(a0) > 0.9),
                     Ephi=phi.mean(), Ephi2=np.mean(phi ** 2),
                     Vphi=phi.var(), third=np.mean(phi ** 2) / 3.0, jk=jk, phi=phi))
    P(f'  L={L}: nlink={a0.size:6d}  plaq_meas={rows[-1]["plaq"]:.4f}  '
      f'max|UU^dag-1|={unit:.2e}  max|Im tr|={aim:.2e}')
    P(f'         E[a0]={rows[-1]["Ea0"]:+.5f}  Var(a0)={rows[-1]["Va0"]:.5f} (Haar 0.25)  '
      f'E[a0^4]={rows[-1]["Ea04"]:.5f} (Haar {H["Ea04"]:.5f})  P(|a0|>.9)={rows[-1]["big"]:.4f} '
      f'(Haar {H["big"]:.4f})')
    P(f'         E[phi]={rows[-1]["Ephi"]:.4f} (pi)  Var(phi)={rows[-1]["Vphi"]:.4f} +- {jk:.4f}  '
      f'(1/3)E[phi^2]={rows[-1]["third"]:.4f}')

P()
P('  reference plaquette at beta=2.5 (GY3, 0.5<1/2 tr U_p>): 0.6526  -- compare the plaq_meas column')
tot = sum(r['n'] for r in rows)
Vall = np.concatenate([r['phi'] for r in rows])
pooled_V = Vall.var()
pooled_E2 = np.mean(Vall ** 2)
se_naive = pooled_V * np.sqrt(2.0 / Vall.size)
se_jk = float(np.sqrt(sum(r['jk'] ** 2 for r in rows)))
per = np.array([r['Vphi'] for r in rows])
se_cfg = float(per.std(ddof=1) / np.sqrt(len(per)))
P(f'  POOLED: nlink={Vall.size}  Var(phi)={pooled_V:.5f}  (1/3)E[phi^2]={pooled_E2/3:.5f}  '
  f'E[phi]={Vall.mean():.5f}')
P(f'          SE: naive iid {se_naive:.5f}; per-config x0-jackknife {se_jk:.5f}; '
  f'config-to-config {se_cfg:.5f} (per-config Var = {", ".join(f"{v:.4f}" for v in per)})')
se = max(se_jk, se_cfg)
P(f'          ADOPTED SE = {se:.5f}')
P(f'          z vs Haar Var(phi)=1.289868134 : {(pooled_V-H["Var"])/se:+.2f}')
P(f'          z vs Haar (1/3)E[phi^2]=3.719824178 : {(pooled_E2/3-H["third"])/se:+.2f}')

hdr('[D] distribution test: empirical CDF / histogram vs Haar (not just two moments)')
edges = np.linspace(0.0, 2 * np.pi, 11)
cnt, _ = np.histogram(Vall, bins=edges)
p = np.array([(HAARW[(PHI >= edges[i]) & (PHI < edges[i + 1])].sum() / HAARW.sum())
              for i in range(10)])
exp = p * Vall.size
chi2 = float(((cnt - exp) ** 2 / exp).sum())
P(f'{"bin [phi]":>16s} {"meas":>7s} {"Haar":>9s} {"dev/sqrt":>8s}')
for i in range(10):
    P(f'  [{edges[i]:5.2f},{edges[i+1]:5.2f}) {cnt[i]:7d} {exp[i]:9.1f} {((cnt[i]-exp[i])/np.sqrt(exp[i])):+8.2f}')
P(f'  chi2 (10 bins, no free parameter) = {chi2:.2f} / dof 10 = {chi2/10:.3f}')
P('  (the Poisson yardstick above is NOT the right error: 3 configs, links correlated.')
fr = np.array([[np.mean((r['phi'] >= edges[i]) & (r['phi'] < edges[i + 1])) for i in range(10)]
               for r in rows])
spread = fr.std(axis=0, ddof=1) / np.sqrt(len(rows))
P('  per-bin frequency: config-to-config SE vs Poisson SE')
P('   bin        Haar   meas    cfg-SE   pois-SE   dev/pois  dev/cfg')
for i in range(10):
    P(f'   {i+1:2d}  {p[i]:9.5f} {cnt[i]/Vall.size:8.5f} {spread[i]:8.5f} '
      f'{np.sqrt(p[i]*(1-p[i])/Vall.size):8.5f} {(cnt[i]-exp[i])/np.sqrt(exp[i]):+9.2f} '
      f'{(cnt[i]/Vall.size-p[i])/spread[i]:+8.2f}')
# KS distance
a0all = np.concatenate(pool_a0)
grid = np.linspace(-1, 1, 4001)
cdf_h = np.array([(HAARW * (np.cos(0.5 * PHI) <= x)).sum() / HAARW.sum() for x in grid])
cdf_e = np.array([np.mean(a0all <= x) for x in grid])
ks = float(np.abs(cdf_h - cdf_e).max())
P(f'  KS distance (a0, measured vs Haar) = {ks:.5f}   (iid-scale 1/sqrt(n) = {1/np.sqrt(a0all.size):.5f})')

hdr('[E] verdict arithmetic')
P(f'  measured Var(2 arccos a0) = {pooled_V:.4f} +- {se:.4f}')
P(f'  branch "flat true" (stated threshold 3.72):  distance {3.719824-pooled_V:+.4f} = '
  f'{(3.719824-pooled_V)/se:+.1f} sigma')
P(f'  branch "saturation" (4.0 - 6.6):            distance {4.0-pooled_V:+.4f} = '
  f'{(4.0-pooled_V)/se:+.1f} sigma (lower edge)')
P(f'  exact Haar prediction for the SAME functional: {H["Var"]:.6f} '
  f'=> deviation {pooled_V-H["Var"]:+.4f} = {(pooled_V-H["Var"])/se:+.1f} sigma')
P(f'  the number the brief calls Var is (1/3)E[phi^2] = {pooled_E2/3:.4f} '
  f'(Haar identity 3.719824, deviation {(pooled_E2/3-H["third"])/se:+.2f} sigma)')
P(f'  closure brackets (Sec B, kappa=2.5 / 3.927): Var(phi) = '
  f'{stats(PHI, HAARW*np.cosh(2.5*np.cos(0.5*PHI)))["Var"]:.3f} / '
  f'{stats(PHI, HAARW*np.cosh(3.927*np.cos(0.5*PHI)))["Var"]:.3f}  '
  f'vs measured {pooled_V:.3f}: excluded' )

with open(os.path.join(HERE, 'ht3_classang.log'), 'w') as fh:
    fh.write('\n'.join(OUT) + '\n')
print('\n[written] research/v3/species/ht3_classang.log')
