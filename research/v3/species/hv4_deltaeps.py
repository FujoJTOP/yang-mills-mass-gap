"""HV4-DELTAEPS -- measure Delta/eps, the whole distance between (G-4) and the target T.

Deterministic only: reads the archived gh2 S(t) arrays / gh2_lt*.json, forms the
connected correlator C(t), and does exact rational algebra on the corpus test
families.  NO new MC, no sampling, no random numbers.

Run:  py -3.12 -I hv4_deltaeps.py
"""
import os
import sys
import json
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = os.path.join(HERE, 'w7')
sys.path.insert(0, W7)
from run3_apebase import correlators, A_at  # noqa: E402

M_FIT = 1.546          # corpus convention for w1 (GH2-LTART / HT5)
NB = 20                # delete-one jackknife blocks
OUT = []


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def jk(Cs, fn):
    v = np.array([fn(c) for c in Cs], float)
    m = float(np.nanmean(v))
    e = float(np.sqrt((NB - 1) / NB * np.nansum((v - m) ** 2)))
    return m, e


# ---------------------------------------------------------------- [1] archives
P('=' * 84)
P('[1] archived gh2 ensembles: C(t)/C(0), Ebar floor, w1   (zero new MC)')
P('=' * 84)
FILES = ['gh2_b2.5_Ls6_Lt12', 'gh2_b2.5_Ls8_Lt12',
         'gh2_b2.5_Ls6_Lt20', 'gh2_b2.5_Ls8_Lt20']
ARCH = {}
for f in FILES:
    d = np.load(os.path.join(W7, f + '.npz'), allow_pickle=True)
    Lt = int(d['Lt']); Ls = int(d['Ls']); S = d['S']
    Sm = S[:, 0, :]
    N = Sm.shape[0]
    C = correlators(Sm, Lt)
    idx = np.array_split(np.arange(N), NB)
    Cs = [correlators(Sm[np.concatenate([idx[j] for j in range(NB) if j != b]), :], Lt)
          for b in range(NB)]
    r = C / C[0]
    r1, r1e = jk(Cs, lambda c: c[1] / c[0])
    elo, eloe = jk(Cs, lambda c: -np.log(max(c[1] / c[0], 1e-300)))
    w1, w1e = jk(Cs, lambda c: A_at(c, Lt, M_FIT) / c[0])
    m1, m2 = r[1], r[2]
    puv = 1.0 - m1 * m1 / (m2 if m2 > 0 else np.nan)
    ARCH[(Ls, Lt)] = dict(Ls=Ls, Lt=Lt, N=N, C0=C[0], r=[float(x) for x in r[:7]],
                          r1=r1, r1e=r1e, elo=elo, eloe=eloe, w1=w1, w1e=w1e,
                          m1=m1, m2=m2, puv=float(puv))
    P(f'{f}: Ls={Ls} Lt={Lt} N={N} C(0)={C[0]:.4f}')
    P('    C(t)/C(0) t=1..6: ' + ' '.join(f'{r[t]:.6f}' for t in range(1, 7)))
    P(f'    Ebar_lo = -ln(C1/C0) = {elo:.5f} +- {eloe:.5f}   (1-C1/C0 = {1 - r[1]:.5f})')
    P(f'    w1 (m={M_FIT}) = {w1:.5f} +- {w1e:.5f}')
    P(f'    Cauchy-Schwarz cap on the weight at s=inf: p_UV <= 1 - m1^2/m2 = {puv:.5f}')

# ------------------------------------------------- replica cross-check (ltA/ltB)
P('-' * 84)
P('[2] replica check: npz == gh2_ltB.json bit-for-bit; gh2_ltA.json = 2nd replica')
P('-' * 84)
REP = {}
for jf in ['gh2_ltA.json', 'gh2_ltB.json']:
    J = json.load(open(os.path.join(W7, jf)))
    for rec in J:
        Ls, Lt = int(rec['Ls']), int(rec['Lt'])
        rr1, rr2 = rec['r']['1'], rec['r']['2']
        puv = 1 - rr1 * rr1 / rr2
        REP[(jf, Ls, Lt)] = dict(r1=rr1, r2=rr2, elo=-np.log(rr1), puv=puv)
        a = ARCH[(Ls, Lt)]
        same = abs(rr1 - a['r1']) < 1e-12
        P(f'  {jf} Ls={Ls} Lt={Lt}: r1={rr1:.6f} r2={rr2:.6f} Ebar_lo={-np.log(rr1):.5f}'
          f'  p_UV<={puv:.5f}   npz r1={a["r1"]:.6f} '
          f'{"IDENTICAL" if same else "independent replica"}')

# --------------------------------------- [3] Jensen lower bound on Ebar (rigour)
P('-' * 84)
P('[3] Jensen: Ebar = <s> = <-ln lam> >= -ln<lam> = -ln(C1/C0)   [exact; needs only nu>=0]')
P('    (lam_n = exp(-s_n) in (0,1]; -ln convex; nothing else used)')
P('-' * 84)
for (Ls, Lt), a in sorted(ARCH.items()):
    P(f'  Ls={Ls} Lt={Lt}: 1-C1/C0 = {1 - a["r1"]:.5f}  <=  Ebar_lo = {a["elo"]:.5f}'
      f'   (Jensen gain x{a["elo"] / (1 - a["r1"]):.2f})')

# ------------------------------------ [4] two-moment extremal family (Markov-Krein)
P('-' * 84)
P('[4] what the two archived moments pin: nu on lam in (0,1], m1=<lam>=C1/C0,')
P('    m2=<lam^2>=C2/C0.  Cap: any weight p below lam=d satisfies p <= 1-m1^2/m2')
P('    (Cauchy-Schwarz), attained at the exact 2-point solution (d=0).  Ebar is')
P('    bounded BELOW by the 2-point solution with the upper atom at lam=1, and is')
P('    UNBOUNDED ABOVE.  => delta1 = Ebar-s1 has a floor and no ceiling.')
P('-' * 84)
for (Ls, Lt), a in sorted(ARCH.items()):
    m1, m2 = a['m1'], a['m2']
    P(f'  Ls={Ls} Lt={Lt}: m1={m1:.6f} m2={m2:.6f}')
    if not (m2 > m1 * m1):
        P('    m2 <= m1^2 -> no 2-point solution (single atom / noise)  skip')
        continue
    var = m2 - m1 * m1
    p0 = 1 - m1 * m1 / m2
    P(f'    d=0 exact: p={p0:.4f} at s=inf, 1-p={1 - p0:.4f} at s={-np.log(m2 / m1):.4f}'
      f'  -> Ebar=inf')
    for lam2 in (1.0, 0.8, 0.5, 0.25, 0.10):
        d = m1 - lam2
        if d >= 0:
            continue
        p = d * d / (var + d * d)
        lam1 = (m1 - (1 - p) * lam2) / p
        if lam1 <= 0 or lam1 > 1:
            P(f'    lam2={lam2}: solution leaves (0,1]  skip')
            continue
        ebar = p * (-np.log(lam1)) + (1 - p) * (-np.log(lam2))
        c1 = p * lam1 + (1 - p) * lam2
        c2 = p * lam1 ** 2 + (1 - p) * lam2 ** 2
        P(f'    lam2={lam2}: p={p:.4f} lam1={lam1:.4f} -> Ebar={ebar:.5f}'
          f'  [m1 err {c1 - m1:+.1e} ; m2 err {c2 - m2:+.1e}]')

# -------------------------------- [5] exact algebra on the corpus families A-E
P('-' * 84)
P('[5] exact rational check: (G-4)_w0  <=>  Delta <= (w1-w0)eps  <=> Delta/eps <= w1-w0')
P('    Landmark 1: Delta/eps = w1-w0 exactly at eps = delta1/(1-w0)  (the G-4 edge)')
P('    Landmark 2: Delta/eps = 0       exactly at eps = delta1/(1-w1)  (eps_max)')
P('-' * 84)


def fam2(name, atoms):
    Wsum = sum(w for w, s in atoms)
    at = [(w / Wsum, s) for w, s in atoms]
    s1 = min(s for w, s in at)
    w1 = sum(w for w, s in at if s == s1)
    rest = [(w, s) for w, s in at if s > s1]
    ebar = sum(w * s for w, s in at)
    s2 = min(s for w, s in rest)
    eps = s2 - s1
    Delta = sum(w * (s - s2) for w, s in rest)
    return dict(name=name, w1=w1, s1=s1, s2=s2, eps=eps, ebar=ebar,
                d1=ebar - s1, Delta=Delta)


CASES = [
    fam2('A  0.6d1+0.4d10', [(F(6, 10), F(1)), (F(4, 10), F(10))]),
    fam2('B  0.01d1+0.99d2', [(F(1, 100), F(1)), (F(99, 100), F(2))]),
    fam2('D  0.9d1+0.1d2', [(F(9, 10), F(1)), (F(1, 10), F(2))]),
    fam2('E  0.9d1+.05d2+.05d100', [(F(9, 10), F(1)), (F(5, 100), F(2)),
                                   (F(5, 100), F(100))]),
    fam2('two-point 0.6d1+0.4d4', [(F(6, 10), F(1)), (F(4, 10), F(4))]),
]
for c in CASES:
    P(f'  {c["name"]}: w1={c["w1"]} s1={c["s1"]} s2={c["s2"]} eps={c["eps"]}'
      f' Ebar={c["ebar"]} d1={c["d1"]} Delta={c["Delta"]}')
    for w0 in (F(1, 2), F(1, 64), F(3, 10)):
        lhs = c['d1'] <= (1 - w0) * c['eps']
        dr = c['Delta'] / c['eps']
        ok = (lhs == (dr <= c['w1'] - w0))
        P(f'     w0={w0}: (G-4)={lhs}  Delta/eps={dr}  w1-w0={c["w1"] - w0}'
          f'  [(G-4)<=>Delta/eps<=w1-w0 : {ok}]')

# continuous-tail family C (exact by hand): (3/5)d1 + (2/5)Unif[4,9]
w1c, s1c, s2c = F(3, 5), F(1), F(4)
ec = s2c - s1c
d1c = (1 - w1c) * (F(13, 2) - s1c)          # tail mean 6.5
Dc = (1 - w1c) * (F(13, 2) - s2c)           # tail spread beyond s2 = 2.5
ebc = d1c + s1c
P(f'  C  (3/5)d1+(2/5)Unif[4,9]: w1={w1c} s1={s1c} s2={s2c} eps={ec} Ebar={ebc}'
  f' d1={d1c} Delta={Dc}')
for w0 in (F(1, 2), F(1, 64), F(3, 10)):
    lhs = d1c <= (1 - w0) * ec
    P(f'     w0={w0}: (G-4)={lhs}  Delta/eps={Dc / ec}  w1-w0={w1c - w0}'
      f'  [(G-4)<=>Delta/eps<=w1-w0 : {lhs == (Dc / ec <= w1c - w0)}]')
# threshold in w1 (family with (1-w1) on the Unif tail): at w1=8/11, Delta/eps = w1-1/2
th = F(8, 11)
P(f'     threshold: at w1={th}, Delta/eps={F(5, 2) * (1 - th) / 3} and w1-1/2={th - F(1, 2)}'
  f'  -> equal: {F(5, 2) * (1 - th) / 3 == th - F(1, 2)}')

# ------------------------- [6] what (G-4) demands of s2, from the archived numbers
P('-' * 84)
P('[6] (G-4)_w0 <=> eps >= delta1/(1-w0) => s2 >= s1 + (Ebar-s1)/(1-w0)')
P('    necessary condition (uses only the rigorous floor Ebar >= Ebar_lo)')
P('-' * 84)
for (Ls, Lt) in sorted(ARCH):
    a = ARCH[(Ls, Lt)]
    for s1 in (1.546, 1.71):
        d1lo = a['elo'] - s1
        P(f'  Ls={Ls} Lt={Lt} s1={s1}: Ebar_lo={a["elo"]:.4f} -> delta1_lo={d1lo:.4f}'
          f'  w1={a["w1"]:.4f}')
        for w0 in (F(1, 64), F(1, 5), F(3, 10)):
            need = s1 + float(d1lo / (1 - float(w0)))
            P(f'      w0={float(w0):.4f} ({w0}): needs s2 >= {need:.4f}'
              f'  (= {need / s1:.3f} s1)')
        eps = s1
        dr = float(d1lo) / eps - (1 - a['w1'])
        P(f'      most favourable eps=s1={eps:.3f}: Delta/eps >= {dr:.4f}'
          f'  => (G-4) at that eps needs w0 <= {a["w1"] - dr:.4f}')

# ----------------------------------------- [7] the window algebra, exactly (C)
P('-' * 84)
P('[7] eps-window at fixed (w1,s1,Ebar): Delta/eps = delta1/eps-(1-w1) strictly')
P('    decreasing in eps; 0 at eps_max = delta1/(1-w1); = w1-w0 at the (G-4) edge.')
P('    family C (w1=3/5, delta1=11/5, eps_max=11/2):')
P('-' * 84)
emax = d1c / (1 - w1c)
for eps in (F(1), F(3), F(9, 2), emax, F(6)):
    dr = d1c / eps - (1 - w1c)
    P(f'    eps={eps} (eps/eps_max={float(eps / emax):.3f}): Delta/eps={dr}'
      f'  => (G-4)_w0 holds iff w0 <= {w1c - dr}')

# --------------------------------------------- [8] the scheme test (exact)
P('-' * 84)
P('[8] scheme test: nu -> nu + sigma*delta_{s_UV}  (what contact subtraction / flow moves)')
P('    eps and s1 do not move; Delta/eps and w1 both move.  Exact rationals.')
P('-' * 84)
base = [(F(3, 5), F(1)), (F(2, 5), F(4))]
for sig, suv in ((F(0), F(0)), (F(1), F(10) ** 4), (F(10 ** 4), F(10) ** 4)):
    at = list(base) + ([(sig, suv)] if sig > 0 else [])
    tot = sum(w for w, s in at)
    W = [(w / tot, s) for w, s in at]
    s1 = F(1)
    w1 = sum(w for w, s in W if s == s1)
    ebar = sum(w * s for w, s in W)
    eps = F(3)
    Delta = ebar - w1 * s1 - (1 - w1) * F(4)
    P(f'    sigma={sig}: w1={w1} ({float(w1):.5f})  eps={eps}  '
      f'Ebar={float(ebar):.4e}  Delta={float(Delta):.4e}'
      f'  Delta/eps={float(Delta / eps):.4e}')

P('=' * 84)
P('ALL SECTIONS DONE (deterministic; no sampling, no new MC)')
P('=' * 84)
os.makedirs(HERE, exist_ok=True)
with open(os.path.join(HERE, 'hv4_deltaeps.log'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(OUT) + '\n')
