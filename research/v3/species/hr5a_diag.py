#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""hr5a_diag.py -- HR5A: chi = Var(a) 的 UV/IR 诊断（精确求积，零新 MC）。

约定（与 GY3/M1 逐字一致）：U = exp(i xi.sigma/2)，xi = phi n, phi=|xi| in [0,2pi)，
a0 = cos(phi/2) = tr U/2，Haar 类密度 ∝ sin^2(phi/2) dphi，
xi_3 = phi n_3，E[n_3^2]=1/3  =>  Var(xi_3) = (1/3) E[phi^2]。

1D 单 link Wilson 类边缘：w_kappa(phi) ∝ sin^2(phi/2) exp(kappa cos(phi/2))；
  Elitzur/中心对称化（单 link 边缘在 U -> -U 下必须偶）：w ∝ sin^2(phi/2) cosh(kappa cos(phi/2))。
i.i.d. 边缘下 Var(a) = Var(xi_3) 逐分量（a = N^{-1/2} sum xi）。
one-link plaquette: <a0> = I_2(kappa)/I_1(kappa)  (exact).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import iv

NQ = 3000
T, WT = leggauss(NQ)                     # t in [-1,1]
PHI = 0.5 * (T + 1.0) * (2.0 * np.pi)    # phi in [0,2pi]
WPHI = WT * np.pi                        # int_0^{2pi} . dphi
SIN2 = np.sin(0.5 * PHI) ** 2


def stat(kappa=0.0, sym=False):
    """returns (a1=<a0>, phi2=<phi^2>, var=Var(xi_3), Z)."""
    c = np.cos(0.5 * PHI)
    w = SIN2 * (np.cosh(kappa * c) if sym else np.exp(kappa * c))
    Z = float(np.dot(WPHI, w))
    a1 = float(np.dot(WPHI, w * c)) / Z
    phi2 = float(np.dot(WPHI, w * PHI ** 2)) / Z
    return a1, phi2, phi2 / 3.0, Z


def kappa_from_plaq(up):
    lo, hi = 0.0, 80.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if stat(mid)[0] < up:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


HAAR = (4 * np.pi ** 2 - 6) / 9.0
print("== [A] Haar exact ==")
a1, phi2, var, Z = stat(0.0)
print(f"  <a0>={a1:.2e}  <phi^2>={phi2:.9f}   closed 4pi^2/3-2={4*np.pi**2/3-2:.9f}")
print(f"  Var(xi_3)={var:.9f}   (4pi^2-6)/9={HAAR:.9f}   1/Var={1/var:.6f}")
print(f"  Var(a) per comp (iid) = Var(xi_3) = {var:.9f}  (same number)")
print(f"  quadrature err = {abs(var-HAAR):.2e}")

BET = [0.5, 1.5, 2.5]
UPL = [0.1234, 0.3619, 0.6526]   # 任务给定 "plaq"（数值上是 <1/2 tr U_p> 口径，见 §2 核对）
KSTAR = {}

print("\n== [B] one-link plaquette g(beta)=I2/I1 vs given plaq ==")
for b, u in zip(BET, UPL):
    g = iv(2, b) / iv(1, b)
    print(f"  beta={b}: g(beta)={g:.6f}  given={u:.4f}  diff={g-u:+.2e}  rel={(g-u)/u:+.2e}")

print("\n== [C] 1D marginal: Var(xi_3), symmetrized version ==")
print(f"  {'kappa':>7} {'g(k)':>8} {'Var':>9} {'Var_sym':>9} {'<phi2>':>9}")
for k in [0.0, 0.5, 1.0, 1.5, 2.5, 3.0, 4.0, 6.0, 10.0]:
    s = stat(k)
    ss = stat(k, sym=True)
    print(f"  {k:7.2f} {s[0]:8.4f} {s[2]:9.4f} {ss[2]:9.4f} {s[1]:9.4f}")

print("\n== [D] beta -> kappa (plaq-matched) and predicted Var(a) ==")
print(f"  {'beta':>5} {'k*':>7} {'g(k*)':>7} {'Var(k*)':>9} {'Varsym(k*)':>10} {'kMFT':>7} {'VarMFT':>8}")
for b, u in zip(BET, UPL):
    ks = kappa_from_plaq(u)
    KSTAR[b] = ks
    A = stat(ks)
    Bv = stat(ks, sym=True)
    um = u ** 0.25
    kd = 6.0 * b * um ** 3
    Dv = stat(kd)
    print(f"  {b:5.2f} {ks:7.4f} {A[0]:7.4f} {A[2]:9.4f} {Bv[2]:10.4f} {kd:7.3f} {Dv[2]:8.4f}")

print("\n== [E] beta=kappa (literal one-link Wilson marginal) ==")
for b in BET:
    s = stat(b)
    ss = stat(b, sym=True)
    print(f"  beta={b}: Var={s[2]:.4f}  Var_sym={ss[2]:.4f}  g(beta)={s[0]:.4f}")

print("\n== [F] measured chi_flat (GY3) vs Haar floor ==")
M = {0.5: [3.770, 3.677, 3.526], 1.5: [3.636, 3.607, 3.863], 2.5: [3.731, 3.753, 3.816]}
SE = {0.5: [0.117, 0.119, 0.106], 1.5: [0.101, 0.109, 0.126], 2.5: [0.067, 0.127, 0.157]}
WMU = {}
for b in BET:
    v = np.array(M[b]); se = np.array(SE[b]); w = 1 / se ** 2
    mu = float(np.sum(w * v) / np.sum(w)); smu = float(1 / np.sqrt(np.sum(w)))
    WMU[b] = (mu, smu)
    print(f"  beta={b}: wmean={mu:.4f}+-{smu:.4f}  spread={v.max()-v.min():.3f} "
          f"dev/Haar={(mu/HAAR-1):+.4f}")

print("\n== [G] predicted ratio Var(2.5)/Var(0.5) vs measured ==")
m0 = WMU[0.5][0]; m2 = WMU[2.5][0]
print(f"  measured: {m0:.4f} -> {m2:.4f}  ratio={m2/m0:.4f}  rel={m2/m0-1:+.4f}  "
      f"se(ratio)~{np.hypot(WMU[0.5][1]/m0, WMU[2.5][1]/m2)*m2/m0:.4f}")
for tag, f in [("A  k=beta", lambda b: stat(b)[2]),
               ("A' k=plaq", lambda b: stat(KSTAR[b])[2]),
               ("B  sym k=beta", lambda b: stat(b, sym=True)[2]),
               ("B' sym k=plaq", lambda b: stat(KSTAR[b], sym=True)[2]),
               ("C  one-plaquette", lambda b: HAAR),
               ("D  MFT k=6b u^3", lambda b: stat(6.0*b*(dict(zip(BET,UPL))[b]**0.25)**3)[2])]:
    v0, v2 = f(0.5), f(2.5)
    print(f"  {tag:16s}: {v0:.4f} -> {v2:.4f}  ratio={v2/v0:.4f}  rel={v2/v0-1:+.4f}")
print("\n  ranges at beta=2.5: A..B ->", f"{stat(2.5)[2]:.3f} .. {stat(2.5,sym=True)[2]:.3f}",
      " (sym big-k limit 2*pi^2/3 =", f"{2*np.pi**2/3:.4f})")
