#!/usr/bin/env python3
"""
ns_nosub.py -- TARGET (NS):   Delta_dyn + W_geo >= lambda      ("no-subtraction" target)

Delta_dyn := E[H]-Var(G)  (= action-only  log-concavity of the conditioned marginal)
W_geo     := E[Phi'']-Var(Phi')
X         := -2Cov(G,Phi')      =>   G'' = Delta_dyn + W_geo + X   and   NOSUB := Delta_dyn+W_geo = G''-X.

PART A  (parse only, NO re-run of any archived witness):
        na_plaquette.log  [P4]        SU(2) 1 plaquette, th=pi, a0=0,+-1
        na2_su3_na.log    [N=3]       SU(3) 1 plaquette, 4 spectra x 5 a0 x 4 beta
        na2_chain.log     [L=2,3,4]   SU(2)^L chain, chain-average coarse field
        -> NOSUB table, sign census, and  min W_geo  vs  the quantitative-Prekopa floor
           lam(v) = inf_xi <v,Hess Phi v> = ||v||^2/6  (SU(2), flat chart, principal branch).
           Chain average mode: v = (1/L)(e3,...,e3)  =>  ||v||^2 = 1/L  =>  lam = 1/(6L).

PART B  (exact quadrature, model A, theta=pi only, ~40 moments() calls = seconds):
        self-test vs P4 archived rows (beta=1,2.5,4,6  x  a0=0,+-1)   [witness validation]
        beta_c^nosub(a0) and beta_c^Gpp(a0) by bisection on [0,7]
        elimination identity  Delta_dyn = dE[G_S]/da + Cov      (tested numerically)

Run: py -3.12 -I ns_nosub.py     (single process, no parallelism, no git)
"""
import json
import re
import numpy as np

OUT = {}


# ============================================================ PART A: parse
def numtok(line):
    line = re.sub(r'\[[^\]]*\]', ' ', line)          # drop [theta spectra]
    line = re.sub(r'(\d)(-\d)', r'\1 \2', line)      # split glued  "0.50-0.50"
    return [float(t) for t in re.findall(r'[-+]?\d+\.\d+(?:[eE][-+]?\d+)?', line)]


def parse_plaquette(fn):
    rows, on = [], False
    for ln in open(fn, encoding='utf-8', errors='ignore'):
        if ln.startswith('[P4]'):
            on = True
            continue
        if on and ln.startswith('[P5]'):
            break
        if on:
            t = numtok(ln)
            if len(t) == 9:
                th, beta, a0, D, W, X, Gpp, fd, dif = t
                rows.append(dict(fam='SU2plaq', th=th, beta=beta, a0=a0,
                                 D=D, W=W, X=X, Gpp=Gpp, nref='P4'))
    return rows


def parse_su3(fn):
    rows, on = [], False
    for ln in open(fn, encoding='utf-8', errors='ignore'):
        if '[N=3]' in ln:
            on = True
            continue
        if on and ln.strip().startswith('====='):
            continue
        if on:
            spec = re.search(r'\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]', ln)
            t = numtok(ln)
            if spec and len(t) >= 9:
                beta, a0, X, se, XM2, D, W, Gpp = t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]
                rows.append(dict(fam='SU3plaq', B=tuple(float(s) for s in spec.groups()),
                                 beta=beta, a0=a0, D=D, W=W, X=X, Gpp=Gpp, nref='N=3'))
    return rows


def parse_chain(fn):
    rows, on = [], False
    for ln in open(fn, encoding='utf-8', errors='ignore'):
        m = re.match(r'\[L=(\d) (open|ring)\]', ln)
        if m:
            on = True
            L, top = int(m.group(1)), m.group(2)
            coarse = None
            continue
        if on and 'coarse field = endpoint' in ln:
            coarse = 'endpoint'
            continue
        if on and 'coarse field = chain average' in ln:
            coarse = 'avg'
            continue
        if on:
            t = numtok(ln)
            if coarse and len(t) == 8:
                beta, a0, X, se, XM2, D, W, Gpp = t
                rows.append(dict(fam='chainL%d' % L, top=top, coarse=coarse,
                                 beta=beta, a0=a0, D=D, W=W, X=X, Gpp=Gpp, nref='[L=%d]' % L))
    return rows


def tbl(rows, title, keep=None):
    print('\n' + '-' * 118)
    print(title)
    print('-' * 118)
    print('%10s %6s %6s %6s %10s %10s %10s %10s   %s'
          % ('fam', 'beta', 'a0', 'extra', 'Ddyn', 'Wgeo', 'X', 'NOSUB', 'sign'))
    for r in rows:
        if keep and not keep(r):
            continue
        N = r['D'] + r['W']
        print('%10s %6.2f %6.2f %6s %10.5f %10.5f %10.5f %10.5f   %s'
              % (r['fam'], r['beta'], r['a0'], '', r['D'], r['W'], r['X'], N,
                 '+' if N > 0 else ('-' if N < 0 else '0')))
        r['NOSUB'] = N


# ============================================================ PART B: model A
class J:
    __slots__ = ('f', 'd', 'dd')

    def __init__(s, f, d=0.0, dd=0.0): s.f, s.d, s.dd = f, d, dd

    def __add__(a, b):
        b = jt(b); return J(a.f + b.f, a.d + b.d, a.dd + b.dd)
    __radd__ = __add__

    def __neg__(a): return J(-a.f, -a.d, -a.dd)

    def __sub__(a, b):
        b = jt(b); return J(a.f - b.f, a.d - b.d, a.dd - b.dd)

    def __rsub__(a, b): return jt(b) - a

    def __mul__(a, b):
        b = jt(b)
        return J(a.f * b.f, a.d * b.f + a.f * b.d, a.dd * b.f + 2 * a.d * b.d + a.f * b.dd)
    __rmul__ = __mul__

    def __truediv__(a, b):
        b = jt(b); i = 1.0 / b.f; i2 = i * i
        return J(a.f * i, a.d * i - a.f * b.d * i2,
                 a.dd * i - 2 * a.d * b.d * i2 - a.f * (b.dd * i2 - 2 * b.d * b.d * i2 * i))

    def __rtruediv__(a, b): return jt(b) / a


def jt(x): return x if isinstance(x, J) else J(x)


def jsqrt(a):
    g = np.sqrt(a.f)
    return J(g, a.d / (2 * g), a.dd / (2 * g) - a.d * a.d / (4 * g * g * g))


def jsin(a): return J(np.sin(a.f), np.cos(a.f) * a.d, -np.sin(a.f) * a.d * a.d + np.cos(a.f) * a.dd)


def jcos(a): return J(np.cos(a.f), -np.sin(a.f) * a.d, -np.cos(a.f) * a.d * a.d - np.sin(a.f) * a.dd)


def jlog(a): return J(np.log(a.f), a.d / a.f, a.dd / a.f - (a.d / a.f) ** 2)


def plaquette(xi, v, beta, theta):
    x = [J(xi[:, i], v[i], 0.0) for i in range(3)]
    rho = jsqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    h = rho * 0.5
    u = x[2] / rho
    S = beta * (1.0 - jcos(h) * np.cos(theta / 2) + jsin(h) * np.sin(theta / 2) * u)
    Phi = 2.0 * jlog(rho * 0.5) - 2.0 * jlog(jsin(h))
    return S, Phi, rho.f


def onb(v):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    a = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(v, a); e1 /= np.linalg.norm(e1)
    return v, e1, np.cross(v, e1)


def slice_nodes(a0, v, ns, nph):
    v, e1, e2 = onb(v)
    R = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-14))
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5 * R * (xs + 1.0); wsv = 0.5 * R * ws
    ph = np.pi * (xp + 1.0); wpv = np.pi * wp
    S, PH = np.meshgrid(s, ph, indexing='ij')
    WS, WPH = np.meshgrid(wsv, wpv, indexing='ij')
    xi = (a0 * v + (S * np.cos(PH))[..., None] * e1 + (S * np.sin(PH))[..., None] * e2)
    return xi.reshape(-1, 3), (WS * S * WPH).reshape(-1), v


def moments(a0, v, beta, theta, ns=4000, nph=48):
    xi, w, vu = slice_nodes(a0, v, ns, nph)
    S, Phi, rho = plaquette(xi, vu, beta, theta)
    lw = np.log(w) - S.f - Phi.f
    mx = lw.max(); ww = np.exp(lw - mx); Z = ww.sum()
    E = lambda x: (ww * x).sum() / Z
    G, H, P1, P2 = S.d, S.dd, Phi.d, Phi.dd
    EG, EP = E(G), E(P1)
    return dict(logrho=np.log(Z) + mx, EH=E(H), VarG=E(G * G) - EG * EG,
                EPhi2=E(P2), VarPhi=E(P1 * P1) - EP * EP,
                Cov=E(G * P1) - EG * EP, EG=EG, EPhi=EP)


def geom(a0, beta, theta, v=(0., 0., 1.), ns=4000, nph=48):
    m = moments(a0, v, beta, theta, ns, nph)
    D = m['EH'] - m['VarG']; W = m['EPhi2'] - m['VarPhi']; X = -2.0 * m['Cov']
    return dict(D=D, W=W, X=X, Gpp=D + W + X, EH=m['EH'], VarG=m['VarG'],
                EG=m['EG'], Cov=m['Cov'])


def bisect(f, lo, hi, it=40):
    flo = f(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi); fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < 1e-6:
            break
    return 0.5 * (lo + hi)


def main():
    W = lambda s: (print(s), OUT.setdefault('log', []).append(s))
    W('=' * 118)
    W('NS  TARGET:  Delta_dyn + W_geo >= lambda        (NOSUB := Delta_dyn + W_geo = G" - X)')
    W('=' * 118)

    # ---------------- PART A
    p = parse_plaquette('na_plaquette.log')
    s3 = parse_su3('na2_su3_na.log')
    ch = parse_chain('na2_chain.log')
    W('\n[PARSED]  SU2plaq=%d   SU3plaq=%d   chain=%d   (no witness re-run)' % (len(p), len(s3), len(ch)))

    # P4 self-table, th=pi only
    tbl([r for r in p if abs(r['th'] - np.pi) < 1e-3], 'A1. SU(2) 1-plaquette, th=pi (P4) -- the (a0,beta) grid')

    # SU(3): collapse by spectrum
    for spec in sorted({r['B'] for r in s3}):
        tbl([r for r in s3 if r['B'] == spec], 'A2. SU(3) 1-plaquette, B=diag(e^{i%s})' % (spec,))

    # chain average only
    tbl([r for r in ch if r['coarse'] == 'avg'], 'A3. SU(2)^L chain, CHAIN-AVERAGE coarse field (the only Phi\'(a=0)!=0 family)')
    tbl([r for r in ch if r['coarse'] == 'endpoint' and r['fam'] == 'chainL2'],
        'A4. control: ENDPOINT coarse field, L=2 (Prop. D locus, X~0)')

    # sign census in the physical window beta<=4.3
    allr = p + s3 + [r for r in ch if r['coarse'] == 'avg']
    for r in allr:
        r['NOSUB'] = r['D'] + r['W']
    win = [r for r in allr if r['beta'] <= 4.3 + 1e-9]
    neg = [r for r in win if r['NOSUB'] < 0]
    W('\n[AP] SIGN CENSUS  beta<=4.3 : %d points, NOSUB<0 at %d (%.0f%%), NOSUB>0 at %d'
      % (len(win), len(neg), 100.0 * len(neg) / max(len(win), 1), len(win) - len(neg)))
    for r in neg:
        W('      NOSUB<0 : %-10s beta=%4.2f a0=%+5.2f  Ddyn=%+8.5f Wgeo=%+8.5f NOSUB=%+8.5f  Gpp=%+8.5f  X=%+8.5f'
          % (r['fam'] + ('/' + ('/'.join('%.2f' % b for b in r['B'])) if 'B' in r else ''),
             r['beta'], r['a0'], r['D'], r['W'], r['NOSUB'], r['Gpp'], r['X']))
    resc = [r for r in win if r['NOSUB'] < 0 < r['Gpp']]
    W('   rescue window (NOSUB<0<Gpp, i.e. X alone keeps G" positive): %d points found in the archive' % len(resc))

    # quantitative-Prekopa floor
    W('\n[AQ] QUANTITATIVE PREKOPA CHECK:   W_geo >= lam(v) = (1/6)||v||^2  ?')
    W('     %-22s %8s %10s %10s %8s' % ('family', 'n', 'min Wgeo', 'lam(v)', 'pass'))
    fams = [('SU2plaq th=pi (||v||^2=1)', [r for r in p if abs(r['th'] - np.pi) < 1e-3], 1.0),
            ('SU2plaq th=0,pi/2 (||v||^2=1)', [r for r in p if abs(r['th'] - np.pi) > 1e-3], 1.0),
            ('SU3plaq (||v||^2=1)', s3, 1.0)]
    for L in (2, 3, 4):
        fams.append(('chainL%d avg (||v||^2=1/L)' % L,
                     [r for r in ch if r['coarse'] == 'avg' and r['fam'] == 'chainL%d' % L], 1.0 / L))
    for nm, rs, nv2 in fams:
        if not rs:
            continue
        mn = min(r['W'] for r in rs)
        W('     %-22s %8d %10.5f %10.5f %8s' % (nm, len(rs), mn, nv2 / 6.0, 'PASS' if mn >= nv2 / 6.0 else 'FAIL'))
    W('     (SU(2) exponent in the flat chart: min of Phi\'\'=1/2 csc^2(rho/2)-2/rho^2 and Phi\'/rho = 1/6)')

    # ---------------- PART B
    W('\n' + '=' * 118)
    W('PART B  exact quadrature, model A (SU(2) 1 plaquette), theta=pi, v=e3, ns=4000 nph=48')
    W('=' * 118)
    thr = np.pi
    W('\n[B0] SELF-TEST vs archived P4 rows (witness must be reproduced, not re-derived):')
    W('     %6s %6s | %10s %10s | %10s %10s' % ('beta', 'a0', 'Ddyn(arch)', 'Ddyn(now)', 'NOSUB(arch)', 'NOSUB(now)'))
    worst = 0.0
    for beta in (1.0, 2.5, 4.0, 6.0):
        for a0 in (0.0, -1.0, 1.0):
            arc = [r for r in p if abs(r['th'] - thr) < 1e-3 and abs(r['beta'] - beta) < 1e-9
                   and abs(r['a0'] - a0) < 1e-9]
            g = geom(a0, beta, thr)
            if arc:
                a = arc[0]
                worst = max(worst, abs(g['D'] - a['D']), abs(g['D'] + g['W'] - (a['D'] + a['W'])))
                W('     %6.2f %6.2f | %10.5f %10.5f | %10.5f %10.5f'
                  % (beta, a0, a['D'], g['D'], a['D'] + a['W'], g['D'] + g['W']))
    W('     max |now - archive| over 12 rows = %.1e   (P4 is printed to 5 decimals; 5e-6 is rounding)' % worst)
    OUT['selftest_max_dev'] = worst

    W('\n[B1] beta-scan, theta=pi:  NOSUB = Delta_dyn+W_geo  and  Gpp = G"+ ... at fixed a0')
    betas = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.3, 5.0, 6.0]
    a0s = [0.0, 0.5, 1.0, 1.5, 2.0]
    scan = {}
    for a0 in a0s:
        W('   a0=%+.2f' % a0)
        W('      %5s %10s %10s %10s %10s %10s' % ('beta', 'Ddyn', 'Wgeo', 'X', 'NOSUB', 'Gpp'))
        for beta in betas:
            g = geom(a0, beta, thr)
            scan['%.2f_%.2f' % (a0, beta)] = g
            W('      %5.2f %10.5f %10.5f %10.5f %10.5f %10.5f'
              % (beta, g['D'], g['W'], g['X'], g['D'] + g['W'], g['Gpp']))
    OUT['scan'] = scan

    W('\n[B2] ZEROS at theta=pi (bisection on [0,7])   --  compare beta* = 4.30661')
    W('      %6s %14s %14s %10s' % ('a0', 'beta_c^NOSUB', 'beta_c^Gpp', 'gap'))
    z = {}
    for a0 in a0s:
        def f1(b):
            g = geom(a0, b, thr); return g['D'] + g['W']

        def f2(b):
            return geom(a0, b, thr)['Gpp']
        try:
            b1 = bisect(f1, 1e-3, 7.0)
        except Exception:
            b1 = float('nan')
        try:
            b2 = bisect(f2, 1e-3, 7.0)
        except Exception:
            b2 = float('nan')
        z['%.2f' % a0] = (b1, b2)
        W('      %6.2f %14.5f %14.5f %10s' % (a0, b1, b2,
                                              ('NOSUB first' if b1 < b2 else 'Gpp first')))
    OUT['zeros'] = z
    b0 = z['0.00'][0]
    W('      [exact] a0=0: X==0 (Prop.D, rank-1+single link) => NOSUB == Gpp identically,')
    W('              zero = sqrt(0.2051152971/0.01105929) = %.5f   (beta*=4.30661)' %
      np.sqrt(0.2051152971 / 0.01105929))

    W('\n[B3] ELIMINATION TEST:  Delta_dyn = dE[G_S]/da + Cov   (na2 identity, my sign convention)')
    h = 1e-3
    W('      %6s %6s %12s %12s %12s' % ('beta', 'a0', 'dEG/da', 'Cov', 'Ddyn-check'))
    for (beta, a0) in ((2.5, 1.0), (4.0, 1.0), (2.5, 0.5), (6.0, -1.0)):
        gm = geom(a0 - h, beta, thr); gp = geom(a0 + h, beta, thr); g0 = geom(a0, beta, thr)
        dEG = (gp['EG'] - gm['EG']) / (2 * h)
        pred = dEG + g0['Cov']
        W('      %6.2f %6.2f %12.6f %12.6f %12.6f   (Ddyn=%+.6f)' % (beta, a0, dEG, g0['Cov'], pred, g0['D']))
    W('      => Delta_dyn is a TOTAL DERIVATIVE **plus Cov**: Cov cannot be eliminated.')
    W('         Equivalently  NOSUB = G" + 2Cov = dE[G]/da + 2Cov.')

    W('\n[BT] BETA-SCALING TEST of the no-go:  Ddyn(beta) = beta*E[H1] - beta^2*Var(G1),  W_geo bounded')
    W('      %6s %10s %10s %12s %12s' % ('a0', 'beta', 'Ddyn', 'Ddyn/beta^2', 'Wgeo'))
    for a0 in (0.0, 1.0):
        for beta in (2.0, 3.0, 4.0, 5.0, 6.0):
            g = geom(a0, beta, thr)
            W('      %6.2f %10.2f %10.5f %12.5f %12.5f' % (a0, beta, g['D'], g['D'] / beta ** 2, g['W']))

    json.dump(OUT, open('ns_nosub.json', 'w'), indent=1)
    txt = '\n'.join(OUT['log']) + '\n'
    open('ns_nosub.log', 'w', encoding='utf-8').write(txt)
    print('\n[wrote ns_nosub.log / ns_nosub.json]   Clay 0%. No QED.')


if __name__ == '__main__':
    main()
