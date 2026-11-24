"""
wi_embedded.py -- WI-EMBEDDED: atom embedded in a continuum.

Resolve WI-LATTICE ("isolation suffices") vs WI-ADMISSIBLE ("1% continuum moves A") on the
REAL physical case: an ATOM riding on a background CONTINUUM (not a pure single level, not a
pure cut).

Model (commensurate, omega=E1=1).  u=e^{-tau}; Z = sum_n c_n u^n (polynomial, c_n>=0):
  - atom at level 1 (E1=1), weight f^2;
  - continuum at levels n=Eth..Nmax (two-particle threshold Eth=2), total weight eps;
  - VEV w0 at level 0.
|Im tau| = |arg u| / omega ;  A = pi/|Im tau| = pi / |arg u|.

Two metrics:
  (N) NEAREST-ZERO  (WI-ADMISSIBLE's metric): the zero with smallest |Im tau| = smallest |arg|
      => largest A.  Pure atom: arg=pi => A=E1=1.
  (T) ATOM-TRACKING (charitable): the zero whose arg is CLOSEST to pi (the atom-like zero);
      its A=pi/|arg|.  Pure atom: exactly pi => A=E1.  With a continuum the atom root is
      perturbed -> A drifts.  Tracking it is the BEST-CASE test of the WI signature.

Reweighting: vary the atom weight f^2 (pure atom is provably invariant to this); also a
fuller family that moves the continuum coupling too.
"""
import numpy as np

def roots_of(levels, weights, w0):
    Nmax = int(max(levels))
    coeffs = np.zeros(Nmax + 1); coeffs[0] = w0
    for n, w in zip(levels, weights):
        coeffs[int(n)] += w
    return np.roots(coeffs[::-1])

def A_nearest(levels, weights, w0):
    r = roots_of(levels, weights, w0)
    if r.size == 0: return np.inf
    a = np.abs(np.angle(r)); a = a[a > 1e-6]
    return np.pi / a.min() if a.size else np.inf

def A_atomtrack(levels, weights, w0, E1=1.0):
    r = roots_of(levels, weights, w0)
    if r.size == 0: return np.inf
    a = np.abs(np.angle(r))
    a = a[a > 1e-6]
    if a.size == 0: return np.inf
    # atom-like zero: arg closest to pi (E1=omega=1 -> arg=pi)
    idx = np.argmin(np.abs(a - np.pi))
    return np.pi / a[idx]

def cont_w(eps, Eth, Nmax, shape='rise'):
    ns = np.arange(Eth, Nmax + 1)
    sh = (ns - Eth + 1).astype(float) if shape == 'rise' else np.ones(len(ns))
    sh = sh / sh.sum()
    return [(int(n), eps*s) for n, s in zip(ns, sh)]

def build(f2, eps, w0=1.0, Eth=2, Nmax=16, shape='rise'):
    cw = cont_w(eps, Eth, Nmax, shape)
    levels = [1] + [n for n, _ in cw]
    weights = [f2] + [w for _, w in cw]
    return levels, weights, w0

F2S = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]

def spread(eps, metric, **kw):
    As = []
    for f2 in F2S:
        lv, wt, w0 = build(f2, eps, **kw)
        As.append((A_nearest if metric=='N' else A_atomtrack)(lv, wt, w0))
    return max(As) - min(As), As

print("="*76)
print("PURE ATOM (eps=0): both metrics must give A=E1=1 for all f^2, spread 0")
print("="*76)
lN, aN = spread(0.0, 'N'); lT, aT = spread(0.0, 'T')
print("  nearest:", [f"{x:.6f}" for x in aN], " spread", f"{lN:.1e}")
print("  track  :", [f"{x:.6f}" for x in aT], " spread", f"{lT:.1e}")

print()
print("="*76)
print("Atom + continuum (Eth=2, rise, Nmax=16): A and spread vs eps")
print("="*76)
print(f"{'eps':>8} | {'A_N(f2=1)':>9} {'spN':>7} | {'A_T(f2=1)':>9} {'spT':>7}")
for eps in [0.0,1e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.3,1.0]:
    sN,_ = spread(eps,'N'); sT,aT = spread(eps,'T')
    _,aN = spread(eps,'N')
    print(f"{eps:8.0e} | {aN[2]:9.4f} {sN:7.4f} | {aT[2]:9.4f} {sT:7.4f}")

print()
print("="*76)
print("Charitable threshold: eps* where ATOM-TRACKING drift crosses 10% of E1 (=0.1)")
print("="*76)
def eps_star(metric, tol, tol_nearest=False, shape='rise', Eth=2):
    lo,hi=0.0,1.0
    for _ in range(60):
        m=0.5*(lo+hi)
        s,_=spread(m, metric, shape=shape, Eth=Eth)
        if s<tol: lo=m
        else: hi=m
    return 0.5*(lo+hi)
for shape in ['rise','flat']:
    for Eth in [2,3]:
        eN=eps_star('N',0.1,shape=shape,Eth=Eth); eT=eps_star('T',0.1,shape=shape,Eth=Eth)
        print(f"  shape={shape:4s} Eth={Eth}:  eps*_nearest={eN:.3f}   eps*_atomtrack={eT:.3f}")

print()
print("="*76)
print("Fuller family: f2(c)=1+3c, eps(c)=eps0*(1+0.5c), c in [0,1]  (atom+continuum move together)")
print("="*76)
for eps0 in [1e-2,3e-2,0.1,0.3]:
    As=[]
    for c in np.linspace(0,1,11):
        lv,wt,w0=build(1+3*c, eps0*(1+0.5*c))
        As.append(A_atomtrack(lv,wt,w0))
    print(f"  eps0={eps0:5.2f}: A_track in [{min(As):.3f},{max(As):.3f}] spread={max(As)-min(As):.3f}")

print()
print("="*76)
print("PURE CUT control: no atom; vary total continuum weight eps, track arg-closest-to-pi.")
print("A_track must MOVE with eps (no conserved zero) -- contrast the atom case (A_T ~ E1 fixed).")
print("="*76)
for eps in [0.1,0.3,1.0,3.0,10.0]:
    lv=[n for n,_ in cont_w(eps,2,16)]; wt=[w for _,w in cont_w(eps,2,16)]
    print(f"  eps={eps:5.2f}: A_track={A_atomtrack(lv,wt,1.0):.4f}   A_nearest={A_nearest(lv,wt,1.0):.4f}")

print()
print("="*76)
print("DISCRIMINATOR = conserved zero across reweighting (track), NOT nearest zero.")
print("="*76)
