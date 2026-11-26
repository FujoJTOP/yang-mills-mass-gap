# -*- coding: utf-8 -*-
r"""
v06_mp.py -- V06 高精度复核 (mpmath 40 dps, 独立于 numpy 版本)
β=2.5, θ=π/2, 单 plaquette SU(2); 全测度 dHaar=4π sin²ψ dψ du, a=0 切片 dμ=2π sin²ψ/ψ dψ.
"""
import math
from mpmath import mp, mpf, sin, cos, exp, sqrt, pi as MPI, legendre, quad

mp.dps = 40
B = mpf("2.5")
CT = sqrt(2) / 2
ST = sqrt(2) / 2


def kern(psi, u):
    rho = 2 * psi
    sp, cp = sin(psi), cos(psi)
    S = B * (1 - CT * cp + ST * sp * u)
    Sr = B * (CT * sp / 2 + ST * cp * u / 2)
    Srr = B * (CT * cp / 4 - ST * sp * u / 4)
    Sru = B * ST * cp / 2
    Su = B * ST * sp
    G = u * Sr + (1 - u * u) / rho * Su
    H = (u * u * Srr + 2 * u * (1 - u * u) / rho * Sru
         + (1 - u * u) / rho * Sr - 3 * u * (1 - u * u) / rho ** 2 * Su)
    Pr = 2 / rho - cp / sp
    Prr = -2 / rho ** 2 + 1 / (2 * sp * sp)
    P1 = u * Pr
    P2 = u * u * Prr + (1 - u * u) / rho * Pr
    return S, G, H, P1, P2


def run(pts_psi, pts_u, slice_only=False):
    # Gauss-Legendre 节点/权重 (mpmath)
    def gl(n):
        x = [legendre(n, mpf(i) / mpf(n)) for i in range(n)]
        # 用 mpmath 的 gauss-legendre 积分求权重: w_i = 2/((1-x^2) P'^2)  -- 直接用 quad 反解太慢,
        # 改用 numpy 节点再转 mpf (节点误差 1e-16, 与 1e-13 的目标相容; 关键是把*被积函数*算到 40 dps)
        import numpy as np
        xn, wn = np.polynomial.legendre.leggauss(n)
        return [mpf(float(t)) for t in xn], [mpf(float(t)) for t in wn]

    xs, ws = gl(pts_psi)
    if slice_only:
        ys, wy = [mpf(0)], [mpf(1)]
    else:
        ys, wy = gl(pts_u)
    psi = [(x + 1) * MPI / 2 for x in xs]
    wpsi = [(w * MPI / 2) for w in ws]
    acc = {}
    keys = ("EH", "EP2", "EG", "EG2", "EP1", "EP1s", "EGP1")
    for k in keys:
        acc[k] = mpf(0)
    Z = mpf(0)
    for i, p in enumerate(psi):
        for j, u in enumerate(ys):
            S, G, H, P1, P2 = kern(p, u)
            if slice_only:
                w = wpsi[i] * wy[j] * sin(p) ** 2 / p * exp(-S)
            else:
                w = wpsi[i] * wy[j] * sin(p) ** 2 * exp(-S)
            Z += w
            acc["EH"] += w * H
            acc["EP2"] += w * P2
            acc["EG"] += w * G
            acc["EG2"] += w * G * G
            acc["EP1"] += w * P1
            acc["EP1s"] += w * P1 * P1
            acc["EGP1"] += w * G * P1
    d = {k: acc[k] / Z for k in keys}
    d["VarG"] = d["EG2"] - d["EG"] ** 2
    d["VarP1"] = d["EP1s"] - d["EP1"] ** 2
    d["Cov"] = d["EGP1"] - d["EG"] * d["EP1"]
    d["EF"] = d["EG"] + d["EP1"]
    d["S"] = d["EH"] + d["EP2"]
    d["D"] = d["VarG"] + d["VarP1"] + 2 * d["Cov"]
    return d


def show(lab, d):
    print("[" + lab + "]  (mpmath, 40 dps)")
    for k in ("EH", "EP2", "VarG", "VarP1", "Cov", "EF"):
        print("   %-8s %s" % (k, mp.nstr(d[k], 22)))
    print("   S        %s" % mp.nstr(d["S"], 22))
    print("   D        %s" % mp.nstr(d["D"], 22))
    print("   S-D      %s" % mp.nstr(d["S"] - d["D"], 6))
    print("")


show("F full", run(200, 100))
show("S slice a=0", run(400, 1, slice_only=True))
