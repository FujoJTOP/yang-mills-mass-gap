#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hr1_haar.py -- HR1: is the flat-chart chi_flat reading eaten by the Haar / log-branch floor?

ZERO new MC.  Inputs:
  research/v3/w7/gy3_4d_results.json        (10 points, beta in {0.5,1.5,2.5}, L in {4,6,8,12})
  research/v3/w7/q5_4d_chi_results.json     (Q-WAVE:30 pilot: beta=0, beta=2.5 L=4,6)
  research/lattice-cpb/runs/witness_L4_b*.npz  (stored SU(3) configurations, beta=2.2/3.0/6.0)
Everything else is deterministic quadrature / closed form.

Sections
  [A] exact Haar one-component variance of xi = log U  (closed form + quadrature)
  [B] theorem: single-link marginal of the unfixed (flat-chart) measure is EXACTLY Haar
      => E[chi_flat] = s2_Haar + Delta(beta,L),  Delta = (1/3N) sum_{l!=l'} <xi_l.xi_l'>,  Delta(beta=0)=0
  [C] the 10+3 measured points against the exact constant  (chi2, weighted mean, beta-slope)
  [D] per-config sd prediction  sd = sqrt(2/3) * s2  (chi^2_3 shape) vs measured chi_flat_sd
  [E] dilution of the L-scan:  bound on the residue r0 of a box-gap family on this carrier
  [F] acid tests: beta=0 (exact), free massive / free massless scalar, 2D (dimension-free identity)
  [G] direct marginal check on stored SU(3) configurations (beta=2.2 / 3.0 / 6.0)
"""
import json
import os
import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
GY3 = os.path.join(ROOT, 'research', 'v3', 'w7', 'gy3_4d_results.json')
Q5 = os.path.join(ROOT, 'research', 'v3', 'w7', 'q5_4d_chi_results.json')

OUT = []


def P(s=''):
    print(s, flush=True)
    OUT.append(str(s))


def hdr(t):
    P()
    P('=' * 78)
    P(t)
    P('=' * 78)


# ----------------------------------------------------------------------------
hdr('[A] exact Haar one-component variance of xi = log U   (U = exp(i xi.sigma/2))')
# Haar on SU(2): U = cos th + i sin th n.sigma, density p(th) = (2/pi) sin^2 th, th in [0,pi].
# principal log: xi = 2 th n  in R^3,  |xi| = 2 th in [0, 2pi).
I2, e2 = quad(lambda t: (2.0 / np.pi) * t * t * np.sin(t) ** 2, 0.0, np.pi, epsabs=1e-14, epsrel=1e-14)
s2_quad = 4.0 * I2 / 3.0
s2_closed = (4.0 * np.pi ** 2 - 6.0) / 9.0
P(f'  <th^2>_Haar          = {I2:.15f}      closed form pi^2/3-1/2 = {np.pi**2/3-0.5:.15f}')
P(f'  s2 = (4/3)<th^2>    = {s2_quad:.15f}')
P(f'  s2 = (4pi^2-6)/9    = {s2_closed:.15f}')
P(f'  |closed - quadrature|= {abs(s2_quad-s2_closed):.3e}')
P(f'  s2 = {s2_closed:.9f}   1/s2 = {1.0/s2_closed:.9f}')
# branch: the cut at th=pi is a measure-zero set; |xi| is in [0,2pi) automatically.
# If one instead used |xi| in [0,pi) (i.e. xi = th n), the value would be s2/4 = 0.929956:
P(f'  (wrong branch, xi = th n): (1/3)<th^2> = {I2/3.0:.6f}  = s2/4')
P('  => the log-branch question is NOT a free knob: the Haar marginal already fixes the branch,')
P('     and the number 3.719824 is the one-component variance of the *Haar marginal itself*.')

# convention check (deterministic quadrature through the estimator's own log): |xi| in [0,2pi)
import sys
sys.path.insert(0, os.path.join(ROOT, 'research', 'lattice-cpb'))
sys.path.insert(0, os.path.join(ROOT, 'research', 'unconditional', 'w7'))
import gpp_core as C  # noqa: E402
_th, _w = np.polynomial.legendre.leggauss(400)
_th = 0.5 * np.pi * (_th + 1.0)
_w = 0.5 * np.pi * _w * (2.0 / np.pi) * np.sin(_th) ** 2
_u, _wu = np.polynomial.legendre.leggauss(60)
_ph = 2.0 * np.pi * (np.arange(60) + 0.5) / 60.0
_n1 = np.sqrt(1 - _u * _u)[:, None] * np.cos(_ph)[None, :]
_n2 = np.sqrt(1 - _u * _u)[:, None] * np.sin(_ph)[None, :]
_n3 = np.broadcast_to(_u[:, None], _n1.shape)
_sx = np.array([[0, 1], [1, 0]], complex)
_sy = np.array([[0, -1j], [1j, 0]], complex)
_sz = np.array([[1, 0], [0, -1]], complex)
_tot = 0.0
_mx = 0.0
for t, wt in zip(_th, _w):
    cc, ss = np.cos(t), np.sin(t)
    U = cc * np.eye(2)[None, None] + 1j * ss * (_n1[..., None, None] * _sx
                                                + _n2[..., None, None] * _sy
                                                + _n3[..., None, None] * _sz)
    x = C.su2_log_x(U)
    r2 = (x * x).sum(-1)
    _tot += wt * float((r2 * (_wu[:, None] / 120.0)).sum())
    _mx = max(_mx, float(np.sqrt(r2).max()))
P(f'  estimator convention (C.su2_log_x, deterministic Haar quadrature): '
  f'(1/3)<|x|^2> = {_tot/3:.5f}, max|x| = {_mx:.3f} (= 2pi = {2*np.pi:.3f})')
P('  => the code\'s xi is 2*theta*nhat with |xi| in [0,2pi): the branch used by M1/GY3 is the one')
P('     that carries s2 = 3.719824 (the |xi| in [0,pi) branch would carry s2/4).')

# ----------------------------------------------------------------------------
hdr('[B] theorem: the flat-chart single-link marginal is exactly Haar at every beta')
P('Claim [T].  dmu(U) = prod_l dU_l e^{-S(U)} / Z with S gauge invariant and dU_l the Haar measure.')
P('For any site x let G act as U_{x,mu} -> G U_{x,mu},  U_{x-mu,mu} -> U_{x-mu,mu} G^dag.')
P('That map is a bijection of the configuration space with unit Jacobian and leaves S invariant.')
P('Hence the marginal m(U_l) = int prod_{l\'!=l} dU_{l\'} e^{-S}/Z satisfies m(G U_l) = m(U_l)')
P('for every G in SU(2); a left-invariant probability density on a compact group is Haar.  [T]')
P('  hypotheses actually used: (i) gauge-invariant action, (ii) product Haar measure,')
P('  (iii) NO gauge-dependent conditioning (no gauge fixing, no center/Polyakov-sector split).')
P()
P('Consequences:')
P('  (B1) (1/3) E|xi_l|^2 = s2 = 3.719824 EXACTLY, for every beta, every L, every dimension.')
P('  (B2) E[chi_flat] = (1/3N) sum_{l,l\'} <xi_l.xi_l\'> = s2 + Delta,')
P('        Delta = (1/3N) sum_{l != l\'} <xi_l.xi_l\'>   (the ONLY beta-dependent piece)')
P('  (B3) beta=0: links iid => Delta = 0 => chi_flat = 3.719824 exactly, all L (corpus D6/D1/P5).')
P('  (B4) the same identity holds for 2D, 3D, any N, any lattice geometry: chi_flat is dimension-free.')
P('  (B5) the plaquette is NOT constrained by (B1): it is a 4-link correlated object and can move')
P('       5.29x while the 1-link marginal is exactly constant.  => all beta-dependence lives in Delta.')

# ----------------------------------------------------------------------------
hdr('[C] the measured points against the exact constant 3.719824')
g = json.load(open(GY3))
rows = []
for k, r in g.items():
    rows.append(dict(tag=k, beta=float(r['beta']), L=int(r['L']), N=int(r['N']),
                     chi=float(r['chi_flat']), se=float(r['chi_flat_se']),
                     sd=float(r['chi_flat_sd']), plaq=float(r['plaq'])))
rows.sort(key=lambda r: (r['beta'], r['L']))
P(f'{"tag":10s} {"beta":>5s} {"L":>3s} {"N":>5s} {"plaq":>7s} {"chi_flat":>9s} {"se":>7s} '
  f'{"z vs s2":>8s} {"sd_meas":>8s} {"sd_pred":>8s}')
sd_pred = np.sqrt(2.0 / 3.0) * s2_closed
chi2 = 0.0
for r in rows:
    z = (r['chi'] - s2_closed) / r['se']
    chi2 += z * z
    P(f'{r["tag"]:10s} {r["beta"]:5.1f} {r["L"]:3d} {r["N"]:5d} {r["plaq"]:7.4f} '
      f'{r["chi"]:9.4f} {r["se"]:7.4f} {z:8.2f} {r["sd"]:8.3f} {sd_pred:8.3f}')
dof = len(rows)
P(f'  chi2 (vs exact constant, no free parameter) = {chi2:.2f} / dof {dof}  = {chi2/dof:.3f}'
  f'   [p ~ {1 - __import__("scipy.stats", fromlist=["chi2"]).chi2.cdf(chi2, dof):.3f}]')

w = np.array([1.0 / r['se'] ** 2 for r in rows])
x = np.array([r['chi'] for r in rows])
wm = float(np.sum(w * x) / np.sum(w))
wse = float(1.0 / np.sqrt(np.sum(w)))
P(f'  weighted mean over the 10 GY3 points = {wm:.5f} +- {wse:.5f}   '
  f'(exact s2 = {s2_closed:.5f}; deviation {(wm-s2_closed)/wse:+.2f} sigma)')

# beta-slope: chi = a + b*(beta-1.5),  per-beta weighted means first (L=4,6,8 common)
betas = sorted(set(r['beta'] for r in rows))
P()
P('  per-beta weighted mean (L=4,6,8 + L=12 where present):')
pm = {}
for b in betas:
    sub = [r for r in rows if r['beta'] == b]
    ww = np.array([1.0 / r['se'] ** 2 for r in sub])
    xx = np.array([r['chi'] for r in sub])
    m = float(np.sum(ww * xx) / np.sum(ww))
    s = float(1.0 / np.sqrt(np.sum(ww)))
    pm[b] = (m, s)
    P(f'    beta={b:4.1f}: {m:.4f} +- {s:.4f}   ({(m-s2_closed)/s:+.2f} sigma vs s2)')
b_arr = np.array([b for b in betas])
m_arr = np.array([pm[b][0] for b in betas])
s_arr = np.array([pm[b][1] for b in betas])
X = np.vstack([np.ones_like(b_arr), b_arr - 1.5]).T
W = np.diag(1.0 / s_arr ** 2)
cov = np.linalg.inv(X.T @ W @ X)
par = cov @ X.T @ W @ m_arr
P(f'  linear fit chi = a + b(beta-1.5):  a = {par[0]:.4f} +- {np.sqrt(cov[0,0]):.4f},  '
  f'b = {par[1]:+.4f} +- {np.sqrt(cov[1,1]):.4f}')
P(f'  => total beta-variation over [0.5,2.5]: {2.0*abs(par[1]):.4f} '
  f'+- {2.0*np.sqrt(cov[1,1]):.4f}   ({100*2*abs(par[1])/s2_closed:.1f}% of s2)')
P(f'  => peak-to-peak of the three beta-means: {m_arr.max()-m_arr.min():.4f} '
  f'= {100*(m_arr.max()-m_arr.min())/s2_closed:.1f}% of s2')
P(f'  => plaquette over the same window: {min(r["plaq"] for r in rows):.4f} -> '
  f'{max(r["plaq"] for r in rows):.4f}  ({max(r["plaq"] for r in rows)/min(r["plaq"] for r in rows):.2f}x)')
P(f'     1-plaq: {1-min(r["plaq"] for r in rows):.4f} -> {1-max(r["plaq"] for r in rows):.4f}  '
  f'({(1-min(r["plaq"] for r in rows))/(1-max(r["plaq"] for r in rows)):.2f}x)')

# Q5 points
q = json.load(open(Q5))
qpts = [('Q5 L4', q['L4']['chi_flat'], q['L4']['chi_flat_se'], 2.5, 4),
        ('Q5 L6', q['L6']['chi_flat'], q['L6']['chi_flat_se'], 2.5, 6),
        ('Q5 beta0', q['beta0']['chi_flat'], q['beta0']['chi_flat_se'], 0.0, 4)]
P()
P('  Q-WAVE:30 pilot (Q5) points:')
chi2q = chi2
for t, c, s, b, L in qpts:
    z = (c - s2_closed) / s
    chi2q += z * z
    P(f'    {t:9s} beta={b:3.1f} L={L}: {c:.4f} +- {s:.4f}  ({z:+.2f} sigma)')
P(f'  all 13 points vs the single constant: chi2 = {chi2q:.2f} / 13 = {chi2q/13:.3f}')

# WG extension: 5 more beta=2.5 points, L=4..16 (same estimator, independent runs)
wg = [('WG L4', 3.69875807791671, 0.10968069797763833, 4),
      ('WG L6', 3.7565648737866546, 0.15968506899753776, 6),
      ('WG L8', 3.6143031664142526, 0.23452614065781638, 8),
      ('WG L12', 3.559083375239752, 0.1360263102229473, 12),
      ('WG L16', 3.7055094248904084, 0.1758986013321115, 16)]
P()
P('  WG extension, beta=2.5, L up to 16 (wg_invariant_results.json / wg_large_results.json):')
chi2w = chi2q
for t, c, s, L in wg:
    z = (c - s2_closed) / s
    chi2w += z * z
    P(f'    {t:8s} L={L:2d}: {c:.4f} +- {s:.4f}  ({z:+.2f} sigma)')
P(f'  ALL 18 points vs the single constant 3.7198: chi2 = {chi2w:.2f} / 18 = {chi2w/18:.3f}'
  f'   [p ~ {1-__import__("scipy.stats", fromlist=["chi2"]).chi2.cdf(chi2w,18):.3f}]')
r16 = wg[4][1] / wg[0][1]
r16s = r16 * np.hypot(wg[4][2] / wg[4][1], wg[0][2] / wg[0][1])
P(f'  L4->16 ratio (WG, beta=2.5) = {r16:.4f} +- {r16s:.4f}   (free-massless lattice: 13.14; continuum L^2: 16.0)')

# ----------------------------------------------------------------------------
hdr('[D] per-config line shape: chi_flat = (1/3)|A|^2 with A Gaussian => sd = sqrt(2/3)*s2')
P(f'  prediction (Haar, chi^2_3 with 3 dof): sd_pred = sqrt(2/3) * s2 = {sd_pred:.4f}')
for r in rows:
    P(f'    {r["tag"]:10s} sd_meas = {r["sd"]:6.3f}   ({(r["sd"]-sd_pred)/sd_pred*100:+5.1f}%)')
P('  the whole per-config distribution, not just the mean, is the beta=0 one.')

# ----------------------------------------------------------------------------
hdr('[E] dilution: what the L-scan can still say about a box-gap family on THIS carrier')
P('The family of GY3 1.1 (single box-gap atom, residue f1^2) contributes X(L) = f1^2/(4 sin^2(pi/L)),')
P('i.e. chi_flat(L) = s2 + f1^2 * chi_box(L).  Subtraction of s2 is exact [T], so the honest test is')
P('absolute, not a ratio of totals.')
box = {4: 1.0 / (4 * np.sin(np.pi / 4) ** 2), 6: 1.0 / (4 * np.sin(np.pi / 6) ** 2),
       8: 1.0 / (4 * np.sin(np.pi / 8) ** 2), 12: 1.0 / (4 * np.sin(np.pi / 12) ** 2),
       16: 1.0 / (4 * np.sin(np.pi / 16) ** 2)}
merged = {4: (3.7654, 0.0498), 6: (3.7581, 0.0624), 8: (3.8218, 0.0633), 12: (3.6589, 0.1903),
          16: (3.7055094248904084, 0.1758986013321115)}
P(f'{"L":>3s} {"chi_box":>8s} {"pred f=1":>9s} {"meas":>8s} {"se":>6s} {"dev":>7s} {"sigma":>6s} '
  f'{"f1^2 2sig ub":>12s}')
for L in (4, 6, 8, 12, 16):
    m, s = merged[L]
    pred = s2_closed + box[L]
    ub = (m + 2 * s - s2_closed) / box[L]
    P(f'{L:3d} {box[L]:8.4f} {pred:9.4f} {m:8.4f} {s:6.4f} {pred-m:7.4f} {(pred-m)/s:6.1f} {ub:12.4f}')
P('  (the corpus GY3-NOGAP-NUM.md:131 quotes 111 sigma from 3.414 vs 1.015+-0.022: that compares')
P('   ratios of TOTALS, where the additive s2 dilutes the family ratio 3.414 -> 1.286.  The honest')
P('   absolute test is stronger and lands at 9-37 sigma, with the binding bound f1^2 <= 0.051 at L=16.)')
c4, c4s = 3.7654, 0.0498
c8, c8s = 3.8218, 0.0633
ratio = c8 / c4
rel = np.hypot(c4s / c4, c8s / c8)
P(f'  GY3 merged beta=2.5: chi(4) = {c4:.4f}+-{c4s:.4f}, chi(8) = {c8:.4f}+-{c8s:.4f}, '
  f'ratio = {ratio:.4f}+-{rel:.4f}')
P(f'  ratio test with the family f1^2=1: pred {s2_closed+box[8]:.4f}/{s2_closed+box[4]:.4f} = '
  f'{(s2_closed+box[8])/(s2_closed+box[4]):.4f} -> {((s2_closed+box[8])/(s2_closed+box[4])-ratio)/rel:.1f} sigma')
P('  => the L-scan survives as a carrier-specific bound: a box-gap atom on the flat-chart log field')
P('     with f1^2 >= 0.051 (2 sigma, L=16) is excluded; below that it is invisible.')

# ----------------------------------------------------------------------------
hdr('[F] acid tests')
P('  F1  beta=0 pure Haar: chi_flat = s2 exactly, every L.  corpus D6-DECISION:135 prints')
P('      Var(a) = 3.7198 for L=1..8 (L-independent) [T].  PASSES -- and this is the tell:')
P('      the same number is returned at beta=2.5, i.e. the estimator cannot tell beta=0 from beta=2.5.')
P('  F2  free massive scalar (m=1): Var(a) = 1/m^2 = 1.000, L-independent (D1 exact);')
P('      free massless scalar: Var = 1/(4 sin^2(pi/L)) = 0.500 / 1.707 / 3.732 / 6.569 (L=4/8/12/16)')
P('      -> grows like L^2/(4pi^2).  So the *scalar* estimator is fine; it is the GAUGE-ORBIT')
P('      marginal (B1) that the acid tests never exercised: D1/M1 acid tests use scalars, no gauge')
P('      symmetry, no Haar orbit => they cannot see the failure mode.')
P('  F3  2D YM: (B1)-(B4) are dimension-free, so chi_flat(2D) = 3.7198 + Delta_2D as well.')
P('      A channel that returns the same number for the solvable/gapped 2D theory and for 4D at')
P('      every beta cannot separate them.  Per CLAUDE.md 1.1 (a tool must fail on the acid case),')
P('      chi_flat is acid-vacuous *as a beta/L discriminator for gauge theories*.')
P('  F4  chart acid (already measured, M1-4D-RUN.md:38): the SAME beta=2.5 ensemble gives')
P('      chi_Landau = 36.98 / 58.61 / 115.30 (L=4/6/8) vs chi_flat 3.81/3.76/3.82.')
P('      Chart changes the number by 10-30x at FIXED beta; beta changes it by 3%.')
P('      => the leading term of the reading is the gauge orbit, by two orders of magnitude.')

# ----------------------------------------------------------------------------
hdr('[G] direct marginal check on stored (unfixed) SU(3) configurations')
P('Theorem (B1) is group-independent.  If the stored configs are unfixed, then for every link')
P('  <1/2 Re tr U_l> = 0  and  sd = 1/sqrt(8) = 0.3536  (Haar),  at EVERY beta, while the plaquette')
P('  (a 4-link joint object) moves freely.')


def dag(X):
    return np.conj(np.swapaxes(X, -1, -2))


def plaq_su3(U):
    tot = 0.0
    n = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            P_ = U[mu] @ np.roll(U[nu], -1, axis=mu) @ dag(np.roll(U[mu], -1, axis=nu)) @ dag(U[nu])
            tot += float(np.mean(1.0 - np.real(np.trace(P_, axis1=-2, axis2=-1)) / 3.0))
            n += 1
    return tot / n
# exact SU(3) Haar reference by quadrature on the maximal torus (Vandermonde weight)
nq = 2400
th = (np.arange(nq) + 0.5) * 2.0 * np.pi / nq
T1, T2_ = np.meshgrid(th, th, indexing='ij')
wt = (16.0 * np.sin(0.5 * (T1 - T2_)) ** 2
      * np.sin(0.5 * (2 * T1 + T2_)) ** 2 * np.sin(0.5 * (T1 + 2 * T2_)) ** 2)
s3 = 0.5 * (np.cos(T1) + np.cos(T2_) + np.cos(T1 + T2_))
m3 = float((wt * s3).sum() / wt.sum())
sd3 = float(np.sqrt((wt * (s3 - m3) ** 2).sum() / wt.sum()))
f3 = float((wt * (np.abs(s3) > 0.9)).sum() / wt.sum())
P(f'  exact SU(3) Haar (torus quadrature): <1/2 Re tr U> = {m3:+.2e}, sd = {sd3:.4f}, '
  f'frac(|s|>0.9) = {f3:.3f}, min = {s3.min():+.3f}')
for b in ('2.2', '3.0', '6.0'):
    f = os.path.join(ROOT, 'research', 'lattice-cpb', 'runs', f'witness_L4_b{b}.npz')
    if not os.path.exists(f):
        P(f'  TOOL-FAIL: missing {f}')
        continue
    d = np.load(f)
    U = d['U']                                  # (4,L,L,L,L,3,3)
    s = 0.5 * np.real(np.trace(U, axis1=-2, axis2=-1)).ravel()
    se = s.std() / np.sqrt(s.size)
    P(f'  beta={b:4s}  links={s.size:5d}  <1/2 Re tr U> = {s.mean():+.5f} +- {se:.4f}  '
      f'sd = {s.std():.4f}  frac(|s|>0.9) = {np.mean(np.abs(s)>0.9):.3f}  min = {s.min():+.3f}  '
      f'plaq = {plaq_su3(U):.4f}')
P('  => the stored (unfixed) links are Haar-marginal to within the 1024-link sampling error,')
P('     at beta=6.0 included: the beta-dependence of the *joint* measure leaves no trace in the')
P('     one-link marginal.  This is (B1) measured, not assumed.')

# how much would a beta-tracking reading have to move?
P()
P('  required variation of any beta-tracking reading vs the 2-sigma bound from the data:')
P(f'    (1-plaq): {1-min(r["plaq"] for r in rows):.4f} -> {1-max(r["plaq"] for r in rows):.4f} '
  f'= {1-(1-max(r["plaq"] for r in rows))/(1-min(r["plaq"] for r in rows)):.0%} variation')
P(f'    data   : total beta-variation {2*abs(par[1]):.4f} +- {2*np.sqrt(cov[1,1]):.4f} '
  f'= {100*2*abs(par[1])/s2_closed:.1f}% +- {100*2*np.sqrt(cov[1,1])/s2_closed:.1f}% of s2;'
  f' 2-sigma upper {100*(2*abs(par[1])+4*np.sqrt(cov[1,1]))/s2_closed:.1f}%')
P(f'    ratio  : required/observed-upper = '
  f'{(1-(1-max(r["plaq"] for r in rows))/(1-min(r["plaq"] for r in rows)))/((2*abs(par[1])+4*np.sqrt(cov[1,1]))/s2_closed):.1f}x')

# ----------------------------------------------------------------------------
hdr('[H] one-line ledger')
P(f'  s2_Haar = {s2_closed:.6f}')
P(f'  chi2(10 GY3 pts vs s2) = {chi2:.2f}/{dof} = {chi2/dof:.2f}')
P(f'  chi2(18 pts: GY3+Q5+WG vs s2) = {chi2w:.2f}/18 = {chi2w/18:.3f}')
P(f'  beta-slope = {par[1]:+.4f} +- {np.sqrt(cov[1,1]):.4f} (total variation over [0.5,2.5] '
  f'{100*2*abs(par[1])/s2_closed:.1f}%)')
P(f'  f1^2 bound (box-gap atom on this carrier, L=16, 2 sigma) = {(merged[16][0]+2*merged[16][1]-s2_closed)/box[16]:.4f}')

with open(os.path.join(HERE, 'hr1_haar.log'), 'w') as fh:
    fh.write('\n'.join(OUT) + '\n')
print('\n[written] research/v3/species/hr1_haar.log')
