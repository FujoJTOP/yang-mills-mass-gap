#!/usr/bin/env python3
"""
na_plaquette.py -- MODEL A: SU(2) one plaquette + fixed boundary holonomy B=exp(i theta sig3/2).

ADJUDICATES (NA):   X := -2 Cov_mu(G, Phi' | a0) >= 0   <=>  Cov <= 0.

Every number printed MUST be read as the 4-tuple
      (model, conditioning, v, a0)
because (NA) is only a statement once those are fixed.

  model        : one plaquette P=UB, free link U=exp(xi), xi in R^3, |xi|<2 pi,
                 S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  u = xi3/rho,
                 flat-chart density  e^{-S} J d^3 xi,  J ~ sin^2(rho/2)/rho^2,  Phi = -log J.
  conditioning : <v,xi> = a0   (a0=0 is the "coarse field = 0" point; since coarse=fine here,
                 conditioning at the coarse field zero forces <v,xi>=0 pointwise)
  v            : unit direction in su(2) (e3 longitudinal / e1 transverse)
  a0           : value of the conditioned coarse coordinate

Object (exact, Lemma G1):
   G''(v;a0) := -d^2/da^2 log rho_v(a)|_{a0} = E[V''|a0] - Var(V'|a0)
              = Delta_dyn + W_geo + X,
   Delta_dyn = E[H]-Var(G),  W_geo = E[Phi'']-Var(Phi'),  X = -2Cov(G,Phi'),
   G = d_v S, H = d_v^2 S, Phi' = d_v Phi, Phi'' = d_v^2 Phi.

DERIVATIVES: computed by an exact second-order forward jet (no finite differences in xi).
Three independent routes to X:
   (1) moment/quadrature   X = -2(E[G Phi'] - E[G]E[Phi'])
   (2) exact identity      X = 2 d/da0 E[Phi'|a0] - 2 W_geo        (5-point FD in a0)
   (3) [-log rho_v]''      5-point FD of the marginal (validates the whole decomposition)

Run: py -3.12 na_plaquette.py
"""
import numpy as np

# ---------------------------------------------------------------- 2-jets
class J:
    __slots__ = ('f', 'd', 'dd')
    def __init__(s, f, d=0.0, dd=0.0): s.f, s.d, s.dd = f, d, dd
    def __add__(a, b):
        b = jt(b); return J(a.f+b.f, a.d+b.d, a.dd+b.dd)
    __radd__ = __add__
    def __neg__(a): return J(-a.f, -a.d, -a.dd)
    def __sub__(a, b):
        b = jt(b); return J(a.f-b.f, a.d-b.d, a.dd-b.dd)
    def __rsub__(a, b): return jt(b)-a
    def __mul__(a, b):
        b = jt(b)
        return J(a.f*b.f, a.d*b.f+a.f*b.d, a.dd*b.f+2*a.d*b.d+a.f*b.dd)
    __rmul__ = __mul__
    def __truediv__(a, b):
        b = jt(b); i = 1.0/b.f; i2 = i*i
        return J(a.f*i, a.d*i-a.f*b.d*i2,
                 a.dd*i-2*a.d*b.d*i2-a.f*(b.dd*i2-2*b.d*b.d*i2*i))
    def __rtruediv__(a, b): return jt(b)/a

def jt(x): return x if isinstance(x, J) else J(x)

def jsqrt(a):
    g = np.sqrt(a.f)
    return J(g, a.d/(2*g), a.dd/(2*g) - a.d*a.d/(4*g*g*g))

def jsin(a): return J(np.sin(a.f), np.cos(a.f)*a.d, -np.sin(a.f)*a.d*a.d+np.cos(a.f)*a.dd)
def jcos(a): return J(np.cos(a.f), -np.sin(a.f)*a.d, -np.cos(a.f)*a.d*a.d-np.sin(a.f)*a.dd)
def jlog(a): return J(np.log(a.f), a.d/a.f, a.dd/a.f-(a.d/a.f)**2)

# ---------------------------------------------------------------- model A
def plaquette(xi, v, beta, theta):
    """jets of (S, Phi) along unit v at base points xi (N,3). Returns (S, Phi, rho)."""
    x = [J(xi[:, i], v[i], 0.0) for i in range(3)]
    rho = jsqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    h = rho*0.5
    u = x[2]/rho
    S = beta*(1.0 - jcos(h)*np.cos(theta/2) + jsin(h)*np.sin(theta/2)*u)
    Phi = 2.0*jlog(rho*0.5) - 2.0*jlog(jsin(h))
    return S, Phi, rho.f

def onb(v):
    v = np.asarray(v, float); v = v/np.linalg.norm(v)
    a = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(v, a); e1 /= np.linalg.norm(e1)
    return v, e1, np.cross(v, e1)

def slice_nodes(a0, v, ns, nph):
    """nodes on {<v,xi>=a0} inside the ball |xi|<2pi, with Lebesgue weight s ds dphi."""
    v, e1, e2 = onb(v)
    R = np.sqrt(max(4*np.pi**2 - a0*a0, 1e-14))
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5*R*(xs+1.0); wsv = 0.5*R*ws
    ph = np.pi*(xp+1.0); wpv = np.pi*wp
    S, PH = np.meshgrid(s, ph, indexing='ij')
    WS, WPH = np.meshgrid(wsv, wpv, indexing='ij')
    xi = (a0*v + (S*np.cos(PH))[..., None]*e1 + (S*np.sin(PH))[..., None]*e2)
    return xi.reshape(-1, 3), (WS*S*WPH).reshape(-1), v

def moments(a0, v, beta, theta, ns=4000, nph=48):
    """conditioned moments at <v,xi>=a0."""
    xi, w, vu = slice_nodes(a0, v, ns, nph)
    S, Phi, rho = plaquette(xi, vu, beta, theta)
    lw = np.log(w) - S.f - Phi.f
    mx = lw.max()
    ww = np.exp(lw-mx); Z = ww.sum()
    E = lambda x: (ww*x).sum()/Z
    G, H, P1, P2 = S.d, S.dd, Phi.d, Phi.dd
    EG, EP = E(G), E(P1)
    return dict(logrho=np.log(Z)+mx, EH=E(H), VarG=E(G*G)-EG*EG,
                EPhi2=E(P2), VarPhi=E(P1*P1)-EP*EP,
                Cov=E(G*P1)-EG*EP, EG=EG, EPhi=EP,
                absPhi1=E(np.abs(P1)), absPhi1max=np.abs(P1).max(), rho=rho)

def decompose(m):
    X = -2.0*m['Cov']
    D = m['EH']-m['VarG']; W = m['EPhi2']-m['VarPhi']
    return dict(Ddyn=D, Wgeo=W, X=X, Gpp=D+W+X)

def fd_curv(a0, v, beta, theta, da=3e-4, ns=4000, nph=48):
    """G'' = -d^2/da^2 log rho_v(a) at a0, 5-point."""
    l = [moments(a0+k*da, v, beta, theta, ns, nph)['logrho'] for k in (-2, -1, 0, 1, 2)]
    return -(-l[0]+16*l[1]-30*l[2]+16*l[3]-l[4])/(12*da*da)

def x_fd(a0, v, beta, theta, da=3e-4, ns=4000, nph=48):
    """X from the exact identity  X = 2 d/da0 E[Phi'|a0] - 2 W_geo."""
    e = [moments(a0+k*da, v, beta, theta, ns, nph)['EPhi'] for k in (-2, -1, 1, 2)]
    dE = (e[0]-8*e[1]+8*e[2]-e[3])/(12*da)
    m = moments(a0, v, beta, theta, ns, nph)
    W = m['EPhi2']-m['VarPhi']
    return 2.0*dE-2.0*W

E3 = np.array([0.0, 0, 1.0]); E1 = np.array([1.0, 0, 0])

def hdr(*a): print(*a, flush=True)

def main():
    hdr("="*96)
    hdr("MODEL A : SU(2) one plaquette, boundary B=exp(i th sig3/2), flat chart, Phi=-log J")
    hdr("(NA): Cov(G,Phi'|a0)<=0  <=>  X=-2Cov >= 0     [4-tuple: model, cond, v, a0]")
    hdr("="*96)

    # ---- [A0] validate the jet machinery against the archived exact constants
    ns = 4000
    m = moments(0.0, E3, 1e-12, np.pi, ns, 96)
    A = m['EPhi2']; B = m['VarG']
    hdr("\n[A0] validation at (th=pi, v=e3, a0=0), beta->0 :  A=E[Phi''] ,  B=Var(G)/beta^2")
    hdr(f"     A = {A:.10f}   (archived 0.20511530)")
    hdr(f"     B = {B:.10f}   (archived 0.01105929)   beta* = {np.sqrt(A/B):.5f} (archived 4.30661)")

    # ---- [P1] the archived table, three methods
    hdr("\n[P1] (th=pi, v=e3, a0=0):  G'' = A - B beta^2   vs FD(-d^2 log rho) vs decomposition")
    hdr(f"     {'beta':>6}{'A-B b^2':>12}{'Gpp decomp':>12}{'fd(-d2logrho)':>15}"
        f"{'Ddyn':>10}{'Wgeo':>10}{'X':>11}")
    for beta in (0.0, 0.5, 1.0, 2.5, 4.0, 6.0):
        be = max(beta, 1e-12)
        ex = A-B*be*be
        m = moments(0.0, E3, be, np.pi, ns, 96); d = decompose(m)
        F = fd_curv(0.0, E3, be, np.pi, ns=ns, nph=96)
        hdr(f"     {beta:6.2f}{ex:12.6f}{d['Gpp']:12.6f}{F:15.6f}"
            f"{d['Ddyn']:10.5f}{d['Wgeo']:10.5f}{d['X']:11.2e}")

    # ---- [P2] THE DEGENERACY: a0=0 kills Phi' pointwise for every v, th, beta
    hdr("\n[P2] DEGENERACY TEST at a0=0 (the coarse-field-zero point of this model):")
    hdr("     Phi is radial => d_v Phi = Phi'(rho) <v,xi>/rho == 0 identically on {<v,xi>=0}.")
    hdr(f"     {'v':>5}{'th':>7}{'beta':>6}{'max|Phi1|':>12}{'E|Phi1|':>11}"
        f"{'Var(Phi1)':>12}{'Cov(G,Phi1)':>13}{'X':>12}")
    worst = 0.0
    for vname, vv in (('e3', E3), ('e1', E1)):
        for th in (0.0, np.pi/2, np.pi):
            for beta in (0.0, 1.0, 2.5, 6.0):
                m = moments(0.0, vv, max(beta,1e-12), th, 3000, 64)
                worst = max(worst, abs(m['Cov']), m['absPhi1max'])
                hdr(f"     {vname:>5}{th:7.3f}{beta:6.2f}{m['absPhi1max']:12.2e}"
                    f"{m['absPhi1']:11.2e}{m['VarPhi']:12.2e}{m['Cov']:13.2e}"
                    f"{-2*m['Cov']:12.2e}")
    hdr(f"     ==> max|Cov| over the whole a0=0 table = {worst:.3e}  (exact zero)")

    # ---- [P3] off-origin sign scan: (NA) has content only here
    hdr("\n[P3] OFF-ORIGIN SIGN SCAN (v=e3).  X=-2Cov(G,Phi'|a0);  X>=0 <=> (NA) holds.")
    for th in (0.0, np.pi/2, np.pi):
        hdr(f"\n  theta={th:.3f}")
        hdr(f"     {'beta':>5}" + "".join(f"{a:>10.2f}" for a in (-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3)))
        for beta in (0.0, 0.5, 1.0, 2.5, 4.0, 6.0):
            be = max(beta, 1e-12)
            row = ""
            for a0 in (-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3):
                m = moments(a0, E3, be, th, 2000, 32)
                row += f"{-2*m['Cov']:>10.4f}"
            hdr(f"     {beta:5.2f}"+row)

    # ---- [P4] full four-term decomposition on the requested grid (a0=0 and a0!=0)
    hdr("\n[P4] DECOMPOSITION GRID  (Ddyn, Wgeo, X, G'') ; v=e3")
    hdr(f"     {'th':>6}{'beta':>6}{'a0':>6}{'Ddyn':>10}{'Wgeo':>10}{'X':>10}"
        f"{'Gpp':>10}{'fd':>10}{'dyn+geo+X-Gpp':>15}")
    for th in (0.0, np.pi/2, np.pi):
        for beta in (0.0, 0.5, 1.0, 2.5, 4.0, 6.0):
            for a0 in (0.0, -1.0, 1.0):
                be = max(beta, 1e-12)
                m = moments(a0, E3, be, th, 3000, 48); d = decompose(m)
                F = fd_curv(a0, E3, be, th, ns=3000, nph=48)
                hdr(f"     {th:6.3f}{beta:6.2f}{a0:6.2f}{d['Ddyn']:10.5f}{d['Wgeo']:10.5f}"
                    f"{d['X']:10.5f}{d['Gpp']:10.5f}{F:10.5f}{-F+d['Ddyn']+d['Wgeo']+d['X']:15.2e}")

    # ---- [P5] task 2: the conditional variance Var(<v,grad S>|0) -- the (M)/H1 direction
    hdr("\n[P5] (M)/H1 DIRECTION: Var(G|a0) and E[H|a0] on the same grid (v=e3)")
    hdr(f"     {'th':>6}{'beta':>6}{'Var(G|0)':>12}{'E[H|0]':>11}{'Ddyn|0':>11}"
        f"{'Gpp|0':>10}{'-Var exposed':>14}")
    for th in (0.0, np.pi/2, np.pi):
        for beta in (0.0, 0.5, 1.0, 2.5, 4.0, 6.0):
            be = max(beta, 1e-12)
            m = moments(0.0, E3, be, th, 3000, 48); d = decompose(m)
            hdr(f"     {th:6.3f}{beta:6.2f}{m['VarG']:12.5f}{m['EH']:11.5f}{d['Ddyn']:11.5f}"
                f"{d['Gpp']:10.5f}{-m['VarG']+m['EH']:14.5f}")

    # ---- [P6] three-route cross-check of X at a few off-origin points
    hdr("\n[P6] THREE ROUTES TO X  (X_mom / X_fd = 2 dE[Phi']/da0 - 2Wgeo / and the FD of Gpp)")
    hdr(f"     {'th':>6}{'beta':>6}{'a0':>6}{'X_mom':>11}{'X_fd':>11}{'diff':>10}"
        f"{'Gpp':>10}{'fd(Gpp)':>10}{'|dGpp|':>10}")
    for th, beta, a0 in ((np.pi, 2.5, 1.0), (np.pi, 6.0, 1.0), (np.pi/2, 2.5, -1.0),
                         (np.pi, 6.0, -2.0), (0.0, 2.5, 1.0), (np.pi, 6.0, 2.5)):
        Xm = -2*moments(a0, E3, beta, th, 4000, 64)['Cov']
        Xf = x_fd(a0, E3, beta, th, ns=4000, nph=64)
        Gp = decompose(moments(a0, E3, beta, th, 4000, 64))
        F = fd_curv(a0, E3, beta, th, ns=4000, nph=64)
        hdr(f"     {th:6.3f}{beta:6.2f}{a0:6.2f}{Xm:11.5f}{Xf:11.5f}{Xm-Xf:10.1e}"
            f"{Gp['Gpp']:10.5f}{F:10.5f}{Gp['Gpp']-F:10.1e}")

    # ---- [P7] quadrature convergence (honest error bars)
    hdr("\n[P7] QUADRATURE CONVERGENCE of the key constants (th=pi,v=e3,a0=1,beta=2.5)")
    for ns2 in (1000, 2000, 4000, 8000):
        m = moments(1.0, E3, 2.5, np.pi, ns2, 48); d = decompose(m)
        hdr(f"     ns={ns2:5d}  Cov={m['Cov']:+.10f}  X={d['X']:+.10f}  Gpp={d['Gpp']:+.10f}"
            f"  EPhi2={m['EPhi2']:.8f}")

    hdr("\nCLay 0%.  No QED.")

if __name__ == '__main__':
    main()
