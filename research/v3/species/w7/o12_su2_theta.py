#!/usr/bin/env python3
r"""
o12_su2_theta.py -- O12: FALSIFY H1.  Candidate 1 (N=2 large beta, is the failure only the
E[H]=0 accident?) and Candidate 2 (theta=pi) and Candidate 5 (boundary chart) on the exact
SU(2) one-plaquette model A of na_plaquette.py.

  Delta := E[curv | 0] - Var[force | 0] = G''(0) = Ddyn + Wgeo + X ,
  Ddyn = E[H]-Var(G), Wgeo = E[Phi'']-Var(Phi'), X = -2Cov(G,Phi').

Model A: S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u), u = xi3/rho,
         flat chart, J = sin^2(rho/2)/rho^2, Phi = -log J, v = e3, a0 = 0.
On a0=0: Phi is radial => Phi' = 0 pointwise => X = 0 exactly, all theta, all beta.
         G = beta sin(th/2) sin(rho/2)/rho ,  H = beta cos(th/2) sin(rho/2)/(2 rho).

Question: after removing the theta=pi accident E[H]==0, is there a beta with Delta<0 and E[H]>0?
Run: py -3.12 -I o12_su2_theta.py
"""
import numpy as np

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
        return J(a.f*i, a.d*i-a.f*b.d*i2, a.dd*i-2*a.d*b.d*i2-a.f*(b.dd*i2-2*b.d*b.d*i2*i))
    def __rtruediv__(a, b): return jt(b)/a

def jt(x): return x if isinstance(x, J) else J(x)
def jsqrt(a):
    g = np.sqrt(a.f); return J(g, a.d/(2*g), a.dd/(2*g)-a.d*a.d/(4*g*g*g))
def jsin(a): return J(np.sin(a.f), np.cos(a.f)*a.d, -np.sin(a.f)*a.d*a.d+np.cos(a.f)*a.dd)
def jcos(a): return J(np.cos(a.f), -np.sin(a.f)*a.d, -np.cos(a.f)*a.d*a.d-np.sin(a.f)*a.dd)
def jlog(a): return J(np.log(a.f), a.d/a.f, a.dd/a.f-(a.d/a.f)**2)

E3 = np.array([0.0, 0.0, 1.0])

def plaquette(xi, v, beta, theta):
    x = [J(xi[:, i], v[i], 0.0) for i in range(3)]
    rho = jsqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    h = rho*0.5; u = x[2]/rho
    S = beta*(1.0 - jcos(h)*np.cos(theta/2) + jsin(h)*np.sin(theta/2)*u)
    Phi = 2.0*jlog(rho*0.5) - 2.0*jlog(jsin(h))
    return S, Phi

def slice_nodes(a0, v, ns, nph):
    R = np.sqrt(max(4*np.pi**2 - a0*a0, 1e-14))
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5*R*(xs+1.0); wsv = 0.5*R*ws
    ph = np.pi*(xp+1.0); wpv = np.pi*wp
    S_, PH = np.meshgrid(s, ph, indexing='ij')
    WS, WPH = np.meshgrid(wsv, wpv, indexing='ij')
    e1 = np.array([1.0, 0.0, 0.0])
    xi = (a0*np.asarray(v, float) + (S_*np.cos(PH))[..., None]*e1
          + (S_*np.sin(PH))[..., None]*np.array([0.0, 1.0, 0.0]))
    return xi.reshape(-1, 3), (WS*S_*WPH).reshape(-1)

def stats(beta, theta, a0=0.0, ns=2000, nph=32):
    xi, w = slice_nodes(a0, E3, ns, nph)
    S, Phi = plaquette(xi, E3, beta, theta)
    lw = np.log(w) - S.f - Phi.f
    mx = lw.max(); ww = np.exp(lw-mx); Z = ww.sum()
    E = lambda x: (ww*x).sum()/Z
    G, H, P1, P2 = S.d, S.dd, Phi.d, Phi.dd
    EG, EP = E(G), E(P1)
    VarG = E(G*G)-EG*EG
    Cov = E(G*P1)-EG*EP
    return dict(logrho=np.log(Z)+mx, EH=E(H), VarG=VarG,
                Wgeo=E(P2)-(E(P1*P1)-EP*EP), X=-2*Cov, absPhi1=np.abs(P1).max())

def delta(beta, theta, **kw):
    m = stats(beta, theta, **kw)
    return m['EH'] - m['VarG'] + m['Wgeo'] + m['X'], m

def fd_curv(beta, theta, a0=0.0, da=3e-4, ns=2000, nph=32):
    l = [stats(beta, theta, a0+k*da, ns, nph)['logrho'] for k in (-2, -1, 0, 1, 2)]
    return -(-l[0]+16*l[1]-30*l[2]+16*l[3]-l[4])/(12*da*da)

def bstar(theta, lo=0.0, hi=40.0, ns=700, nph=16):
    """first beta>0 where Delta changes sign from + to - (scan + bisection)."""
    f = lambda b: delta(b, theta, ns=ns, nph=nph)[0]
    if theta == 0.0:
        return float('inf')
    prev = f(lo+1e-6); a = lo+1e-6
    b = a
    while b < hi:
        b = min(b+1.0, hi)
        v = f(b)
        if v < 0 and prev > 0:
            for _ in range(40):
                m = 0.5*(a+b); vm = f(m)
                if vm > 0: a, prev = m, vm
                else: b = m
            return 0.5*(a+b)
        a, prev = b, v
    return float('nan')

def main():
    print("="*104)
    print("O12 / CAND 1,2 : SU(2) one plaquette, a0=0, v=e3.  Delta = E[H]-Var(G)+Wgeo+X = G''(0)")
    print("="*104)
    # --- validation against the archive (na_plaquette.log [A0])
    A = stats(1e-12, np.pi, ns=4000, nph=96)
    print(f"\n[V] beta->0, theta=pi :  Wgeo = E[Phi''] = {A['Wgeo']:.8f}  (archived 0.20511530)")
    print(f"                              Var(G)/b^2 = {A['VarG']:.8f}  (archived 0.01105929)"
          f"  => beta* = {np.sqrt(A['Wgeo']/A['VarG']):.5f} (archived 4.30661)")

    print("\n[F1] theta x beta grid.  EH>0 everywhere except theta=pi (exact zero).")
    print(f"     {'theta':>7}{'beta':>6}{'E[H]':>11}{'Var(G)':>11}{'Wgeo':>10}{'X':>9}"
          f"{'Delta':>11}{'fd(-d2logrho)':>15}{'sign':>6}")
    thetas = [0.0, np.pi/4, np.pi/2, 3*np.pi/4, 0.90*np.pi, 0.95*np.pi, np.pi]
    betas = [1.0, 2.5, 4.0, 5.0, 6.0, 10.0]
    neg = []
    for th in thetas:
        for be in betas:
            D, m = delta(be, th, ns=1500, nph=24)
            F = fd_curv(be, th, ns=900, nph=16) if be in (2.5, 5.0) else float('nan')
            sgn = 'NEG' if D < 0 else 'pos'
            if D < 0:
                neg.append((th, be, D, m['EH']))
            print(f"     {th:7.3f}{be:6.2f}{m['EH']:11.5f}{m['VarG']:11.5f}{m['Wgeo']:10.5f}"
                  f"{m['X']:9.2e}{D:11.5f}{F:15.5f}{sgn:>6}")

    print("\n[F2] beta*(theta): the Delta=0 crossing (only theta=pi has E[H] identically 0)")
    print(f"     {'theta':>8}{'beta*':>9}{'1-tanh?':>10}")
    for th in thetas:
        bs = bstar(th)
        print(f"     {th:8.3f}{bs:9.4f}")

    print("\n[F3] does Delta<0 survive after E[H]>0 is enforced?  (theta=pi excluded)")
    phys = [(t, b, D, eh) for (t, b, D, eh) in neg if b <= 3.0]
    print(f"     number of (theta,beta) cells with Delta<0 over the grid: {len(neg)}")
    print(f"     of those with E[H]>1e-9 : {sum(1 for r in neg if r[3] > 1e-9)}")
    print(f"     of those with beta<=3.0 : {len(phys)}   {phys}")
    print("\n     reading: at theta<pi E[H] = beta cos(th/2) * 0.2212... > 0 strictly, so any")
    print("     Delta<0 there is NOT the E[H]==0 accident.  The question is only WHERE (beta*).")
    print("\nClay 0%. No QED.")

if __name__ == '__main__':
    main()
