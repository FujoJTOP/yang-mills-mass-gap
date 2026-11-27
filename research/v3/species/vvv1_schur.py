#!/usr/bin/env python3
"""vvv1_schur.py -- VVV1-SCHUR witness. Deterministic, zero MC, no sampling.

Toy (a) : nu(ds) = s^alpha ds on [0,1], alpha known.
Toy (b) : SSS5's own witness nu_eps = (1-eps) d_{s=1} + eps d_{s=eps^2}.
Object  : Schur test  limsup|a_n| < 1  (SSS5 sect 1.2)  vs  alpha.
  a_n := Schur parameters (Schur algorithm) of the Caratheodory function
         F(z) = 1 + 2 sum_{k>=1} c_k z^k  of the pushed-forward circle measure,
         phi = (F-1)/(F+1),  a_n = gamma_n  (|a_n| <= 1).
Run: py -3.12 -I research/v3/species/vvv1_schur.py
"""
import mpmath as mp
import sympy as sp
from fractions import Fraction as F

mp.mp.dps = 40


# ---------- exact arithmetic: nu = s^alpha ds on [0,1] ----------------------
def Delta(n, a):
    """det(1/(i+j+1+a))_{0<=i,j<=n} (Cauchy) -- exact Fraction."""
    num, fk = F(1), F(1)
    for k in range(1, n + 1):
        fk *= k
        num *= fk * fk
    den = F(1)
    for i in range(n + 1):
        for j in range(n + 1):
            den *= (i + j + 1 + a)
    return num / den


def X_exact(n, a):                       # h_n/h_{n-1}, h_n = ||P_n||^2 monic
    D = [Delta(k, a) for k in range(n + 1)]
    h = [D[k] / D[k - 1] for k in range(1, n + 1)]
    return h[-1] / h[-2]


def X_formula(n, a):
    return (F(n * n) * (n + a) ** 2
            / (F(2 * n + 1 + a) * F(2 * n + a) ** 2 * F(2 * n - 1 + a)))


print("[A] exact  X_n = h_n/h_{n-1}  vs closed form (Cauchy determinant)")
for a in (F(0), F(1, 2), F(1), F(2)):
    ok = all(X_exact(n, a) == X_formula(n, a) for n in range(2, 7))
    print(f"    alpha={str(a):>4}: exact==formula n=2..6 : {ok}"
          f"   X_10={float(X_formula(10,a)):.12f}")

n, al, t = sp.symbols('n alpha t', positive=True)
Xs = n**2 * (n + al)**2 / ((2*n + 1 + al) * (2*n + al)**2 * (2*n - 1 + al))
SER = sp.series(Xs.subs(n, 1/t), t, 0, 4).removeO().expand()
print("[A] X_n =", sp.simplify(SER), "  => alpha enters only at O(1/n^2), even in alpha")

# ---------- finite part of the edge transform -------------------------------
print("[B] G(z)+log(1-z) -> -(psi(alpha+1)+gamma) , z=1-eps  (sign of alpha)")
for a in (mp.mpf(0), mp.mpf('0.5'), mp.mpf(1), mp.mpf(2)):
    tgt = -(mp.digamma(a + 1) + mp.euler)
    vals = []
    for e in (mp.mpf('1e-8'), mp.mpf('1e-16'), mp.mpf('1e-30')):
        z = 1 - e
        G = (-mp.log(1-z)/z if a == 0
             else mp.quad(lambda s: s**a/(1 - z*s), [0, 1]))
        vals.append(G + mp.log(1-z))
    print(f"    alpha={float(a):<4} target={mp.nstr(tgt,8):>11}"
          f"  vals={[mp.nstr(v,8) for v in vals]}")


# ---------- series helpers + Schur algorithm --------------------------------
def smul(u, v, N):
    w = [mp.mpc(0)]*(N+1)
    for i, ui in enumerate(u):
        if ui == 0:
            continue
        for j, vj in enumerate(v[:N+1-i]):
            w[i+j] += ui*vj
    return w


def sinv(u, N):                          # 1/u, u[0] != 0
    w = [mp.mpc(0)]*(N+1)
    w[0] = 1/u[0]
    for k in range(1, N+1):
        w[k] = -sum(u[j]*w[k-j] for j in range(1, k+1))/u[0]
    return w


def schur_params(c, N):
    """c[k]=c_k/c_0 (k=0..N).  Returns gamma_0..gamma_{N-1} of phi=(F-1)/(F+1)."""
    Fs = [mp.mpf(2)*c[k] for k in range(N+1)]
    Fs[0] = 1
    num = [Fs[k] - (1 if k == 0 else 0) for k in range(N+1)]
    den = [Fs[k] + (1 if k == 0 else 0) for k in range(N+1)]
    phi = smul(num, sinv(den, N), N)
    g, out = phi, []
    for _ in range(N):
        gam = g[0]
        out.append(gam)
        if abs(gam) > 1 - mp.mpf('1e-30'):
            break
        psi = [g[k+1] for k in range(N)] + [mp.mpc(0)]          # (g-gam)/z
        d = [1 - gam.conjugate()*g[k] for k in range(N+1)]      # 1-gam*g
        g = smul(psi, sinv(d, N), N)
    return out


def fourier(w, N):
    c = [mp.quad(lambda tt: w(tt)*mp.cos(k*tt), [0, mp.pi/2, mp.pi])/mp.pi
         for k in range(N+1)]
    return [c[k]/c[0] for k in range(N+1)]


def ctrl():
    N = 12
    flat = schur_params([mp.mpc(0)]*(N+1), N)
    th = mp.mpf('0.7')
    atom = schur_params([mp.mpc(1)] + [mp.exp(-mp.mpc(0, 1)*th*k)
                                       for k in range(1, N+1)], N)
    print(f"[C] controls: flat max|a_n|={mp.nstr(max(abs(x) for x in flat),3)} (=0);"
          f" 1-atom |a_0|={mp.nstr(abs(atom[0]),3)} (=0) |a_1|={mp.nstr(abs(atom[1]),8)}"
          f" (=1) then recursion stops (extreme point), #{len(atom)} params")


ctrl()

# ---------- toy (a): arc map s=(1-cos th)/2, w = s^alpha sin(th)/2 ----------
print("[D] toy (a) ARC map  nu(ds)=s^alpha ds -> beta=alpha+1/2 (th=0), 1/2 (th=pi)")
print("    exact mechanism: w=0 off the arc => Re F=0 there => |phi|=1 there => |a_n|->1")
mp.mp.dps = 60
N = 40
for a in (mp.mpf(0), mp.mpf('0.5'), mp.mpf(1), mp.mpf(2)):
    w = (lambda aa: (lambda th: mp.sin(th/2)**(2*aa)*mp.sin(th)/2))(a)
    g = schur_params(fourier(w, N), N)
    seq = []
    for k in range(0, min(6, len(g))):
        seq.append(f"{k}:{mp.nstr(abs(g[k]),6)}")
    print(f"    alpha={float(a):<4} |a_n| = {' '.join(seq)}"
          f"   #{len(g)} params before 1-|a|^2<1e-30")
    if len(g) >= 5:
        print(f"              1-|a_n| = "
              f"{[mp.nstr(1-abs(g[k]),3) for k in range(1,5)]}")

# ---------- toy (b): SSS5's nu_eps  (truth alpha = -1) ----------------------
print("[E] toy (b) nu_eps=(1-eps)d_{s=1}+eps d_{s=eps^2}, arc map; truth alpha=-1")
for e in (mp.mpf('1e-2'), mp.mpf('1e-4'), mp.mpf('1e-6')):
    th2 = mp.acos(1 - 2*e**2)
    p1, p2 = 1 - e, e
    c = [mp.mpc(1)] + [(p1*mp.exp(-mp.mpc(0, 1)*mp.pi*k)
                        + p2*mp.exp(-mp.mpc(0, 1)*th2*k))/(p1 + p2)
                       for k in range(1, 17)]
    g = schur_params(c, 16)
    seq = [f"{k}:{mp.nstr(abs(g[k]),8)}" for k in range(min(4, len(g)))]
    print(f"    eps={mp.nstr(e,2):>7} |c_1|={mp.nstr(abs(c[1]),8)} (=1-2eps)"
          f"  |a_n| {seq}  #{len(g)} params")
print("[E] exact: 2 atoms => F rational deg 2 => phi inner rational => |a_2|=1 exactly ;")
print("     hence limsup=1 at EVERY eps>0 with alpha=-1 (test cannot see alpha here)")

