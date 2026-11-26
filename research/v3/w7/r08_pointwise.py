# -*- coding: utf-8 -*-
"""r08_diag.py -- 定位 |Δ|^2Φ'' 或 |Δ|^2Φ'^2 的多项式编码错误 (点态核对)"""
import mpmath as mp
mp.mp.dps = 50

LAM = (1, 0, -1)


def J(t):
    s = 1
    for i, j in ((0, 1), (0, 2), (1, 2)):
        s *= 4 * mp.sin((t[i] - t[j]) / 2) ** 2
    return s


def phi_p(t):
    return sum(-mp.mpf(LAM[i] - LAM[j]) / mp.tan((t[i] - t[j]) / 2) for i, j in ((0, 1), (0, 2), (1, 2)))


def phi_pp(t):
    return sum(mp.mpf(1) / 2 * (LAM[i] - LAM[j]) ** 2 / mp.sin((t[i] - t[j]) / 2) ** 2
               for i, j in ((0, 1), (0, 2), (1, 2)))


def Jtau(t, tau):
    return J([t[i] + tau * LAM[i] for i in range(3)])


t = [mp.mpf('0.7'), mp.mpf('-1.1')]
t3 = -(t[0] + t[1])
th = [t[0], t[1], t3]
JJ = J(th)
Pp = phi_p(th)
Ppp = phi_pp(th)
d2J = mp.diff(lambda tau: Jtau(th, tau), 0, 2)
print("J            =", mp.nstr(JJ, 20))
print("Phi'         =", mp.nstr(Pp, 20))
print("Phi''        =", mp.nstr(Ppp, 20))
print("J(Phi'^2-Phi'')   =", mp.nstr(JJ * (Pp ** 2 - Ppp), 20))
print("d^2_tau J (num)   =", mp.nstr(d2J, 20))
print("residual          =", mp.nstr(JJ * (Pp ** 2 - Ppp) - d2J, 5))
