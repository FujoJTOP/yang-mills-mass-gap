# -*- coding: utf-8 -*-
r"""
t06_su3_gap.py -- T06: SU(3), theta=pi, beta=2.5, a=0 的 gap 复核 + 加密 2x/4x + 区间.

框架 (r12, 采用):
  <f> = int |Delta|^2 e^{-S_W} f (da db/3)/Z ,  |Delta|^2 = prod_p 4 sin^2(d_p/2), d=(a,b,a+b)
  Phi = sum_p g(d_p), g(d) = -2 log[sin(d/2)/(d/2)]
  S = <H + Phi''>,  D = Var(G + Phi'),  gap = S - D
  G = (beta/3) sum_j sin(phi_j - x_j) u_j,  H = (beta/3) sum_j cos(phi_j - x_j) u_j^2
  u = torus 方向 (sum u = 0);  lambda=(1,0,-1) <=> u=(1,0,-1) <=> (a,b) 平面 v=(1,1)
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 40
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


BETA = 2.5
PHI = (math.pi, math.pi, 0.0)


def g(d):
    return -2.0 * mp.log(mp.sin(d / 2) / (d / 2))


def g1(d):
    return 2 / d - mp.cos(d / 2) / mp.sin(d / 2)


def g2(d):
    s = mp.sin(d / 2)
    return -2 / d ** 2 + 1 / (2 * s * s)


def part1():
    P("=" * 104)
    P("t06_su3_gap.log -- T06: SU(3) theta=pi, beta=2.5, a=0 的 gap 复核 / 加密 / 区间")
    P("=" * 104)
    P()
    P("[1] Phi 的核 g(d) = -2 log[sin(d/2)/(d/2)] 在 d->0 (task 点 2):")
    P("    %-8s %-22s %-22s %-22s" % ("d", "g(d)/(d^2/12)", "g'' (应 1/6)", "g(d)-(d^2/12+d^4/1440)"))
    for dd in ("1e-1", "1e-2", "1e-3", "1e-5"):
        D = mp.mpf(dd)
        P("    %-8s %-22s %-22s %-22s" % (dd, mp.nstr(g(D) / (D ** 2 / 12), 14),
                                          mp.nstr(g2(D), 12), mp.nstr(g(D) - (D ** 2 / 12 + D ** 4 / 1440), 6)))
    P("    => g(d) = d^2/12 + d^4/1440 + O(d^6),  g'(0)=0,  g''(0)=1/6:  d^2 起步的偶解析函数.")
    P("    [T] g 在 d=0 无奇性 => Phi=sum_p g(d_p) 在 d_p=0 处解析; 无 'd=0 奇性障碍'.")
    P("        边界上 |Delta|^2 ~ d_p^2 使 T*Phi' (O(d^2)), T*Phi'' (O(d^2)) 均光滑 (见 [3] 逐点表).")
    P()


def integ(N, u=(1, 0, -1), weight='SW', tau=0.0, z3=False, BLK=512, fold=False):
    """中点网格 N x N on (a,b) in (0,2pi)^2, a=x1-x2, b=x2-x3.  分块防内存.
    z3=True: 再对 Z3 中心像 x -> x + k*(2pi/3)(1,1,1) 取平均 (补回 torus 的 2/3).
    fold=True: d_p 折到 (-pi,pi] (g 含 log(d/2) => Phi 不是 torus 函数, lift 约定敏感)."""
    h = 2 * math.pi / N
    ax = (np.arange(N) + 0.5) * h
    u1, u2, u3 = u
    dd = {'12': u1 - u2, '23': u2 - u3, '13': u1 - u3}
    va, vb = dd['12'], dd['23']
    ph1, ph2, ph3 = PHI
    T3 = 0.0
    ks = (0.0, 2 * math.pi / 3, 4 * math.pi / 3) if z3 else (0.0,)
    norm_k = 3.0 if z3 else 1.0
    acc = {k: 0.0 for k in ('Z', 'TPhi', 'TPhi2', 'TPhipp', 'TG', 'TH', 'TG2', 'TInvG', 'TInv')}
    for i0 in range(0, N, BLK):
        A = np.broadcast_to(ax[i0:i0 + BLK, None] + tau, (min(BLK, N - i0), N))
        B = np.broadcast_to(ax[None, :] + tau, A.shape)
        d12, d23, d13 = A, B, A + B
        s12 = 4 * np.sin(d12 / 2) ** 2
        s23 = 4 * np.sin(d23 / 2) ** 2
        s13 = 4 * np.sin(d13 / 2) ** 2
        T = s12 * s23 * s13
        x1 = (2 * d12 + d23) / 3
        x2 = (d23 - d12) / 3
        x3 = -(d12 + 2 * d23) / 3
        q1, q2, q0 = {}, {}, {}
        for nm, v in (('12', d12), ('23', d23), ('13', d13)):
            vv = ((v + np.pi) % (2 * np.pi)) - np.pi if fold else v
            q1[nm] = 8 * np.sin(v / 2) ** 2 / vv - 2 * np.sin(v)
            q2[nm] = (4 * np.sin(v / 2) / vv - 2 * np.cos(v / 2)) ** 2
            q0[nm] = (2 * np.sin(v / 2) / vv) ** 2
        s = {'12': s12, '23': s23, '13': s13}
        others = {'12': ('23', '13'), '23': ('12', '13'), '13': ('12', '23')}
        TPhi = np.zeros_like(T)
        TPhi2 = np.zeros_like(T)
        TPhipp = np.zeros_like(T)
        for p in ('12', '23', '13'):
            o1, o2 = others[p]
            TPhi += dd[p] * q1[p] * s[o1] * s[o2]
            TPhi2 += dd[p] ** 2 * q2[p] * s[o1] * s[o2]
            TPhipp += dd[p] ** 2 * (2 - 2 * q0[p]) * s[o1] * s[o2]
        for pp, qq in (('12', '23'), ('12', '13'), ('23', '13')):
            rr = ({'12', '23', '13'} - {pp, qq}).pop()
            TPhi2 += 2 * dd[pp] * dd[qq] * q1[pp] * q1[qq] * s[rr]
        invd = 1 / d12 ** 2 + 1 / d23 ** 2 + 1 / d13 ** 2
        for k in ks:
            c1 = np.cos(ph1 - x1 - k)
            c2 = np.cos(ph2 - x2 - k)
            c3 = np.cos(ph3 - x3 - k)
            s1 = np.sin(ph1 - x1 - k)
            s2 = np.sin(ph2 - x2 - k)
            s3 = np.sin(ph3 - x3 - k)
            G = (BETA / 3) * (s1 * u1 + s2 * u2 + s3 * u3)
            H = (BETA / 3) * (c1 * u1 ** 2 + c2 * u2 ** 2 + c3 * u3 ** 2)
            ew = np.exp(-BETA * (1 - (c1 + c2 + c3) / 3)) if weight == 'SW' else np.ones_like(T)
            w = T * ew
            acc['Z'] += float(np.sum(w)) / norm_k
            acc['TPhi'] += float(np.sum(ew * TPhi)) / norm_k
            acc['TPhi2'] += float(np.sum(ew * TPhi2)) / norm_k
            acc['TPhipp'] += float(np.sum(ew * TPhipp)) / norm_k
            acc['TG'] += float(np.sum(w * G)) / norm_k
            acc['TH'] += float(np.sum(w * H)) / norm_k
            acc['TG2'] += float(np.sum(w * G * G)) / norm_k
            acc['TInvG'] += float(np.sum(ew * G * TPhi)) / norm_k
            acc['TInv'] += float(np.sum(w * invd)) / norm_k
    Z = acc['Z']
    o = {'Z': Z}
    o['EPhi'] = acc['TPhi'] / Z
    o['EPhi2'] = acc['TPhi2'] / Z
    o['EPhipp'] = acc['TPhipp'] / Z
    o['EG'] = acc['TG'] / Z
    o['EH'] = acc['TH'] / Z
    o['EG2'] = acc['TG2'] / Z
    o['EGPhi'] = acc['TInvG'] / Z
    o['Einv'] = acc['TInv'] / Z
    o['S'] = o['EPhipp'] + o['EH']
    o['D'] = o['EG2'] + o['EPhi2'] + 2 * o['EGPhi'] - (o['EG'] + o['EPhi']) ** 2
    o['gap'] = o['S'] - o['D']
    o['Eext'] = 2 * (va ** 2 * o['Einv'] + vb ** 2 * o['Einv'] + (va + vb) ** 2 * o['Einv'])  # 占位, 见 Eext2
    return o


def Eext(N, u=(1, 0, -1), weight='SW', z3=False, BLK=512):
    """E(v) = 2 < v_a^2/a^2 + v_b^2/b^2 + (v_a+v_b)^2/(a+b)^2>"""
    h = 2 * math.pi / N
    ax = (np.arange(N) + 0.5) * h
    u1, u2, u3 = u
    va, vb = u1 - u2, u2 - u3
    ph1, ph2, ph3 = PHI
    ks = (0.0, 2 * math.pi / 3, 4 * math.pi / 3) if z3 else (0.0,)
    nk = 3.0 if z3 else 1.0
    Z = 0.0
    Ex = 0.0
    for i0 in range(0, N, BLK):
        A = np.broadcast_to(ax[i0:i0 + BLK, None], (min(BLK, N - i0), N))
        B = np.broadcast_to(ax[None, :], A.shape)
        d12, d23, d13 = A, B, A + B
        T = 4 * np.sin(d12 / 2) ** 2 * 4 * np.sin(d23 / 2) ** 2 * 4 * np.sin(d13 / 2) ** 2
        x1 = (2 * d12 + d23) / 3
        x2 = (d23 - d12) / 3
        x3 = -(d12 + 2 * d23) / 3
        for k in ks:
            c1 = np.cos(ph1 - x1 - k)
            c2 = np.cos(ph2 - x2 - k)
            c3 = np.cos(ph3 - x3 - k)
            ew = np.exp(-BETA * (1 - (c1 + c2 + c3) / 3)) if weight == 'SW' else np.ones_like(T)
            w = T * ew
            Z += float(np.sum(w)) / nk
            Ex += float(np.sum(w * (va ** 2 / d12 ** 2 + vb ** 2 / d23 ** 2
                                    + (va + vb) ** 2 / d13 ** 2))) / nk
    return 2 * Ex / Z


def part2():
    P("[2] r12 复现 (中点 256..2048) 与加密 2x/4x (task 点 1).  v=lambda=(1,0,-1), a=0, beta=2.5:")
    P("    %-8s %-15s %-15s %-15s %-11s %-11s" % ("N", "S", "D", "gap", "<Phi'>", "<G>"))
    st = {}
    for N in (256, 512, 1024, 2048, 4096, 8192):
        o = integ(N)
        st[N] = o
        P("    %-8d %-15.8f %-15.8f %+-15.8f %+-11.3e %+-11.3e"
          % (N, o['S'], o['D'], o['gap'], o['EPhi'], o['EG']))
    P("    r12 报: S=6.80640  D=6.70468  gap=+0.10172   (N=2048, 同网格 => 应精确复现)")
    P()
    P("    分解 (N=8192): <Phi''>=%+.8f <H>=%+.8f | <G^2>=%+.8f <Phi'^2>=%+.8f 2<GPhi'>=%+.8f"
      % (st[8192]['EPhipp'], st[8192]['EH'], st[8192]['EG2'], st[8192]['EPhi2'], 2 * st[8192]['EGPhi']))
    P()
    P("    纯 Haar 权重 |Delta|^2 (r12 分支二, r12 报 S=6.44071 D=6.67651 gap=-0.23580):")
    P("    %-8s %-15s %-15s %-15s" % ("N", "S0", "D0", "gap0"))
    st2 = {}
    for N in (1024, 2048, 4096, 8192):
        o = integ(N, weight='Haar')
        st2[N] = o
        P("    %-8d %-15.8f %-15.8f %+-15.8f" % (N, o['S'], o['D'], o['gap']))
    P()
    return st, st2


def part3(st, st2):
    P("[3] 恒等式结构 + 平移不变性诊断 (task 点 3):")
    P("    严格恒等式 (对 -log T = Phi_r08 成立): <H + A''> = Var(G + A'), A=-log|Delta|^2.")
    P("    对 Phi_O09 = A + psi + const, psi = 2 sum_p log(d_p/2):")
    P("       gap = S - D = <psi''> - Var(psi') - 2 Cov(G+A', psi')")
    P("           = -2<sum_p v_p^2/d_p^2> - Var(psi') - 2 Cov(...)")
    P("    %-8s %-14s %-16s %-16s" % ("N", "D-S", "-<psi''>", "残差(D-S+<psi''>)"))
    for N in (1024, 2048, 4096, 8192):
        o = st[N]
        Em = 2 * (Eext(N) / 2)          # = 2<sum v^2/d^2> = -<psi''>
        P("    %-8d %+-14.6f %-16.6f %-16.6f" % (N, o['D'] - o['S'], -Em, (o['D'] - o['S']) + Em))
    P()
    P("    平移不变性 Z(tau)=int rho(a+tau,b+tau) 必须与 tau 无关 (真 Weyl 积分的必要条件):")
    P("    %-8s %-12s %-22s %-22s" % ("N", "tau", "Z(tau)", "Z(tau)/Z(0)-1"))
    for N in (1024, 2048):
        z0 = integ(N)['Z']
        for tau in (0.0, 0.5, 1.0, math.pi):
            z = integ(N, tau=tau)['Z']
            P("    %-8d %-12.4f %-22.10f %+.3e" % (N, tau, z, z / z0 - 1))
    P("    => 非零随 tau 漂移 => (a,b) in (0,2pi)^2 不是 torus 的基本域: 它只覆盖 x-torus 的 1/3.")
    P("       (x -> (a,b) 是 3 对 1; 缺 Z3 中心像 x -> x + k(2pi/3)(1,1,1), k=1,2.)")
    P()
    P("    Z3 平均后 (补回 2/3 覆盖):")
    P("    %-8s %-15s %-15s %-15s %-11s %-11s" % ("N", "S", "D", "gap", "<Phi'>", "<G>"))
    for N in (1024, 2048, 4096):
        o = integ(N, z3=True)
        P("    %-8d %-15.8f %-15.8f %+-15.8f %+-11.3e %+-11.3e"
          % (N, o['S'], o['D'], o['gap'], o['EPhi'], o['EG']))
    for N in (2048,):
        z0 = integ(N, z3=True)['Z']
        for tau in (0.5, math.pi):
            z = integ(N, z3=True, tau=tau)['Z']
            P("    Z3 平移核验 N=%d tau=%.3f: Z/Z0-1 = %+.3e" % (N, tau, z / z0 - 1))
    P()


def part4():
    P("[4] G 的 v 依赖 (task 点 3): v = (a,b) 平面方向, u = ((2va+vb)/3,(vb-va)/3,-(va+2vb)/3)")
    P("    %-8s %-10s %-10s %-15s %-15s %-12s" % ("psi", "va", "vb", "S", "gap", "gap/S"))
    rows = []
    for psi in (0.0, 0.4, math.pi / 4, 1.2, math.pi / 2):
        va, vb = math.cos(psi), math.sin(psi)
        u = ((2 * va + vb) / 3, (vb - va) / 3, -(va + 2 * vb) / 3)
        o = integ(2048, u=u)
        rows.append((psi, o))
        P("    %-8.4f %-10.6f %-10.6f %-15.8f %+-15.8f %-12.5f"
          % (psi, va, vb, o['S'], o['gap'], o['gap'] / o['S']))
    P("    gap 是 v 的 2 次齐次:  gap(c v) = c^2 gap(v)  (N=1024, psi=0, v=(c,0)):")
    for c in (0.5, 1.0, 2.0):
        va, vb = c, 0.0
        u = ((2 * va + vb) / 3, (vb - va) / 3, -(va + 2 * vb) / 3)
        o = integ(1024, u=u)
        P("      c=%-5.2f  gap=%+-15.8f  gap/c^2=%+-15.8f  S/c^2=%.8f" % (c, o['gap'], o['gap'] / c ** 2, o['S'] / c ** 2))
    P()
    return rows


if __name__ == "__main__":
    part1()
    st, st2 = part2()
    part3(st, st2)
    part4()
    P("[5] lift 约定 x 定义域 的敏感性 (N=2048, v=lambda):")
    P("    %-28s %-15s %-15s %-15s" % ("约定", "S", "D", "gap"))
    for lab, kw in (("d=(a,b,a+b), (0,2pi)^2 [r12]", dict()),
                    ("d 折到 (-pi,pi]", dict(fold=True)),
                    ("Z3 平均 (真 torus), lift 原样", dict(z3=True)),
                    ("Z3 平均 + d 折到 (-pi,pi]", dict(z3=True, fold=True))):
        o = integ(2048, **kw)
        P("    %-28s %-15.8f %-15.8f %+-15.8f" % (lab, o['S'], o['D'], o['gap']))
    with open("t06_su3_gap.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
