# -*- coding: utf-8 -*-
r"""
y03_vstatus.py -- Y03: `v` 的地位：把 SU(3) 的 gap(v) / E[H_vv](v) 当 (v_a,v_b) 平面上的
二次型重建，取单位圆上的本征值 ==> inf_v / sup_v。

复用 t06 的积分器 integ()（不重写；中点网格，|Delta|^2 e^{-S_W} 权重）：
  S = <H + Phi''>, D = Var(G + Phi'), gap = S - D,  E[H_vv] 由 o['EH'] 给出.
  v=(v_a,v_b) <-> u = ((2va+vb)/3,(vb-va)/3,-(va+2vb)/3) ; v=(1,1) <-> u=lambda=(1,0,-1).

两套域: z3=False = (0,2pi)^2（r12/t06 的 1/3 域）; z3=True = Z3 补全（真 Weyl 积分）.
判据: H1 要 "every coarse transverse v" => 需要 min_{||v||=1} gap(v) > 0.
本脚本只产出数; 判决与逐字引文在 research/v3/species/Y03-DIRECT.md.
"""
import sys
import math

sys.path.insert(0, r"D:\Dev\math\research\v3\w7")
import t06_su3_gap as t6  # noqa: E402

BETA = 2.5
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


def u_of(va, vb):
    return ((2 * va + vb) / 3, (vb - va) / 3, -(va + 2 * vb) / 3)


def three(N, z3):
    """v=(1,0),(0,1),(1,1) 三点：gap 与 E[H] 各三个数（二次型需 3 点）。"""
    out = {}
    for tag, (va, vb) in (("aa", (1.0, 0.0)), ("bb", (0.0, 1.0)), ("ab", (1.0, 1.0))):
        o = t6.integ(N, u=u_of(va, vb), z3=z3)
        out[tag] = (o["gap"], o["EH"], o["S"], o["D"], o["EPhipp"])
    return out


def mat(f):
    """Q(v) = A va^2 + 2B va vb + C vb^2 ; 由 Q(1,0),Q(0,1),Q(1,1) 定 A,B,C。"""
    A, B_, C = f["aa"], f["bb"], f["ab"]
    B = (C - A - B_) / 2.0
    tr = A + B_
    disc = math.sqrt(((A - B_) / 2.0) ** 2 + B * B)
    return (A, B, B_), (tr / 2.0 - disc, tr / 2.0 + disc)


def main():
    P("=" * 100)
    P("y03_vstatus.log -- Y03: v 的地位 = SU(3) 上 gap(v)/E[H](v) 的二次型本征值（theta=pi, beta=2.5, a0=0）")
    P("=" * 100)
    P("定义: Q_gap(v) = gap(v) = S(v) - D(v); Q_H(v) = E[H_vv](v).  v 为单位向量时 Q 即该方向的数.")
    P("H1 逐字量词 = 'for every coarse transverse v at k->0' => 需要 min_{||v||=1} Q_gap(v) > 0.")
    P()
    res = {}
    for z3 in (False, True):
        lab = "Z3 补全 (真 Weyl)" if z3 else "(0,2pi)^2 (1/3 域)"
        for N in ((1024, 2048) if z3 else (1024, 2048, 4096)):
            f = three(N, z3)
            (Aq, Bq, Cq), (lmin, lmax) = mat({k: v[0] for k, v in f.items()})
            (Ah, Bh, Ch), (hmin, hmax) = mat({k: v[1] for k, v in f.items()})
            res[(z3, N)] = dict(f=f, gap=(Aq, Bq, Cq, lmin, lmax), H=(Ah, Bh, Ch, hmin, hmax))
            P("[%s] N=%d" % (lab, N))
            P("   gap 型:  A=%+.8f  B=%+.8f  C=%+.8f   => lambda_min=%+.8f  lambda_max=%+.8f"
              % (Aq, Bq, Cq, lmin, lmax))
            P("   E[H]型:  A=%+.8f  B=%+.8f  C=%+.8f   => lambda_min=%+.8f  lambda_max=%+.8f"
              % (Ah, Bh, Ch, hmin, hmax))
            P("   逐点: v=(1,0): gap=%+.8f E[H]=%+.8f | v=(0,1): gap=%+.8f E[H]=%+.8f | v=(1,1): gap=%+.8f E[H]=%+.8f"
              % (f["aa"][0], f["aa"][1], f["bb"][0], f["bb"][1], f["ab"][0], f["ab"][1]))
        P()
    # 二次型校验：第 4/5 方向（不在构造里）
    z3 = True
    N = 2048
    d = res[(z3, N)]
    Aq, Bq, Cq = d["gap"][:3]
    Ah, Bh, Ch = d["H"][:3]
    P("[校验] Z3 补全 N=2048: 用 A,B,C 预测 vs 直接积分（二次型若成立，两者应一致）")
    P("   %-8s %-12s %-14s %-14s %-14s %-14s" % ("psi", "v", "gap_pred", "gap_direct", "EH_pred", "EH_direct"))
    for psi in (0.4, 1.2):
        va, vb = math.cos(psi), math.sin(psi)
        o = t6.integ(N, u=u_of(va, vb), z3=True)
        gp = Aq * va ** 2 + 2 * Bq * va * vb + Cq * vb ** 2
        hp = Ah * va ** 2 + 2 * Bh * va * vb + Ch * vb ** 2
        P("   %-8.4f (%.4f,%.4f) %+-14.8f %+-14.8f %+-14.8f %+-14.8f"
          % (psi, va, vb, gp, o["gap"], hp, o["EH"]))
    P()
    # 单位圆上的极小/极大所在方向
    P("[极值方向] Z3 补全 N=2048（theta = 0.5*atan2(2B, A-C) 给出极值方向）")
    for nm, (A, B, C) in (("gap", (Aq, Bq, Cq)), ("E[H]", (Ah, Bh, Ch))):
        th = 0.5 * math.atan2(2 * B, A - C)
        for t, w in ((th, "lambda_min/或 max@"), (th + math.pi / 2, "正交@")):
            va, vb = math.cos(t), math.sin(t)
            P("   %-5s %s psi=%.6f  v=(%+.6f,%+.6f)  Q=%+.8f"
              % (nm, w, t, va, vb, A * va ** 2 + 2 * B * va * vb + C * vb ** 2))
    P()
    # c_j = <cos x_j> 反解（E[H] 对 u_j^2 线性）—— 用 1<->2 交换对称性独立验测度
    import numpy as np
    P("[c_j 反解 + 1<->2 对称性检验]  phi=(pi,pi,0) 在 x1<->x2 下不变 => 真测度必须 c1=c2.")
    for z3 in (False, True):
        d = res[(z3, 2048)]
        M = []
        rhs = []
        for tag, (va, vb) in (("aa", (1.0, 0.0)), ("bb", (0.0, 1.0)), ("ab", (1.0, 1.0))):
            u1, u2, u3 = u_of(va, vb)
            M.append([-u1 ** 2, -u2 ** 2, u3 ** 2])
            rhs.append(d["f"][tag][1] / (BETA / 3))
        c = np.linalg.solve(np.array(M), np.array(rhs))
        P("   z3=%-5s c1=%+.8f c2=%+.8f c3=%+.8f   |c1-c2|=%.2e  => %s"
          % (str(z3), c[0], c[1], c[2], abs(c[0] - c[1]),
             "对称, 可用" if abs(c[0] - c[1]) < 1e-5 else "破 1<->2 对称 => 域不是真测度"))
    P()
    P("[判] inf_{||v||=1} gap = lambda_min; sup = lambda_max.  若 lambda_min<0 => 'every v' 的 H1 在该测度上被单方向打穿.")
    for key in sorted(res, key=lambda k: (k[0], k[1])):
        d = res[key]
        P("   z3=%-5s N=%-5d  gap: min=%+.8f max=%+.8f | E[H]: min=%+.8f max=%+.8f"
          % (str(key[0]), key[1], d["gap"][3], d["gap"][4], d["H"][3], d["H"][4]))
    with open(r"D:\Dev\math\w7\y03_vstatus.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
