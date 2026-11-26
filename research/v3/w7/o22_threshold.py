"""O22 — 阈值 beta_0 = V_inf /(c ||v||^2): 算出来并定位它在哪。

策略 (LOOP-RUNBOOK §11.64):
    Var(G_v+Phi'|0) <= V_inf(theta)      [o05 实测: beta-无关上确界]
    E[H]           >= c*beta*||v||^2      [o21 若在; 否则 c = N/12 = Phi'' 底 (乐观)]
    => H1 对一切 beta >= beta_0 := V_inf/(c||v||^2) 成立

模型: 完整复用 o05 (单 plaquette SU(2), v=e3, 条件化 xi_3=0, 精确一维求积):
    S = beta(1 - cos(r/2)cos(th/2)),  w(r)dr = 2pi sin^2(r/2)/r * e^{-S} dr
    G = beta sin(th/2) sin(r/2)/r,   H = beta cos(th/2) sin(r/2)/(2r)
    Phi''_vv = 2/r^2 - cot(r/2)/r    (r->0 极限 = 1/6 = N/12, N=2)
    Var(G_v+Phi'|0) = Var(G|0),  E[H+Phi''|0] = E[H|0] + E[Phi''|0]

只写 o22_*.py|.log（禁 git commit）。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o05_var_scaling import slice_stats  # noqa: E402  (复用 o05 的精确切片)

NPTS = 400001
BIGS = [200.0, 400.0]                    # Var 随 beta 单调增 => 这两个值给 V_inf 的数值下界
C_SU2 = 2.0 / 12.0                       # N/12, N=2  (Lemma B 的 Phi'' 底 = 0.1666667)
C_SU3 = 3.0 / 12.0                       # N/12, N=3  = 0.25 (物理组, 见 §窗口)
CHECK_THETAS = [np.pi / 2, 3 * np.pi / 4]
THETAS = [0.0, np.pi / 8, np.pi / 4, 3 * np.pi / 8, np.pi / 2, 5 * np.pi / 8,
          3 * np.pi / 4, 7 * np.pi / 8, 15 * np.pi / 16, 0.97 * np.pi,
          0.98 * np.pi, 0.99 * np.pi, 0.995 * np.pi]


def vinf(th, n=NPTS):
    """V_inf(th) = sup_beta Var(beta,th) 的数值近似 (Var 单调增 => 取最大 beta) + cH。"""
    out = [slice_stats(b, th, n=n) for b in BIGS]
    var = max(o["var"] for o in out)
    eh = out[-1]["EH"]
    return var, eh, eh / BIGS[-1], out[-1]["var"]


def beta0(th, c):
    v, _, _, _ = vinf(th)
    return v / c, v


def crossing(target, lo, hi, c, it=40):
    """在 [lo,hi] 上找 beta_0(th)=target 的根 (beta_0 在 th 上单调增, 见 [P1])。"""
    flo = beta0(lo, c)[0] - target
    fhi = beta0(hi, c)[0] - target
    if flo > 0 or fhi < 0:
        return float("nan")
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if beta0(mid, c)[0] - target < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    print("[O22] beta_0 = V_inf/(c*||v||^2); 单 plaquette SU(2) 精确求积 (复用 o05)")
    print(f"[O22] 网格 n={NPTS}; V_inf 取 beta in {BIGS} 的较大者 (Var 单调增)")
    print(f"[O22] c = N/12: SU(2) {C_SU2:.7f}, SU(3) {C_SU3:.7f}; ||v||^2=1 (o05: v=e3)")

    # ---- [CHK] 与 o05 存档对账
    print("\n[CHK] V_inf 对账 (O05-FORCEVAR §3: theta=pi/2 -> 0.0278, 3pi/4 -> 0.1618)")
    for th in CHECK_THETAS:
        v, eh, cH, v400 = vinf(th)
        print(f"  theta={th:.6f}  V_inf={v:.6f}  (beta=400: {v400:.6f})  E[H]/beta={cH:.6f}")

    # ---- [P1] beta_0 主表
    print("\n[P1] beta_0(theta) 主表  (V_inf 见 [CHK]; beta_0 = 6*V_inf for c=N/12|SU(2))")
    print("  theta/pi   V_inf      E[H]/beta   b0(c=1/6)  b0(c=1/4)  b0(c=1/6+cH)   vs 2.5   vs 6")
    rows = []
    for th in THETAS:
        v, eh, cH, _ = vinf(th)
        b16 = v / C_SU2
        b14 = v / C_SU3
        bmix = v / (C_SU2 + cH)
        rows.append((th, v, cH, b16, b14, bmix))
        print(f"  {th/np.pi:8.4f} {v:10.6f}  {cH:10.6f}   {b16:10.5f}  {b14:9.5f}"
              f"   {bmix:11.5f}     {'<=' if b16 <= 2.5 else ' >'}    "
              f"{'<=' if b16 <= 6.0 else ' >'}")

    # ---- [P1b] c 的三种取法 (口径敏感性): N/12 作 beta-系数 | 实测 cH | 含 beta-无关底
    print("\n[P1b] c 口径敏感性: b0 = V_inf/c  (c_N12 = N/12 当 beta-系数;"
          " c_H = 实测 E[H]/beta; c_full: RHS=beta*c_H+N/12 => b0 = max(0,(V_inf-N/12)/c_H))")
    print("  theta/pi   V_inf      c_H        b0(N/12)   b0(c_H)    b0(c_full)   vs 2.5 (c_H)")
    for th in THETAS:
        v, _, cH, _ = vinf(th)
        bN, bH = v / C_SU2, (v / cH if cH > 1e-12 else float("inf"))
        bF = max(0.0, (v - C_SU2) / cH) if cH > 1e-12 else float("inf")
        print(f"  {th/np.pi:8.4f} {v:10.6f}  {cH:9.6f}   {bN:10.5f}  {bH:10.5f}   {bF:11.5f}"
              f"     {'<=' if bH <= 2.5 else ' >'}")
    for tgt in (2.5, 6.0):
        thc = crossing(tgt, 0.0, np.pi * 0.9995, C_SU2)
        # c_H 口径: 单独二分
        lo, hi = 0.0, np.pi * 0.9995
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            v, _, cH, _ = vinf(mid)
            if v / cH < tgt:
                lo = mid
            else:
                hi = mid
        print(f"  c_H 口径: beta_0={tgt:4.1f} <=> theta_c = {0.5*(lo+hi)/np.pi:.5f} pi"
              f"  (N/12 口径: {thc/np.pi:.5f} pi)")

    # ---- [P2] 阈值定位: beta_0 = 2.5 / 6 的临界面
    print("\n[P2] 临界面: beta_0(theta_c)=2.5 与 =6  (c=1/6; 扫描单调, 二分)")
    for tgt in (2.5, 6.0):
        thc = crossing(tgt, 0.0, np.pi * 0.9995, C_SU2)
        print(f"  beta_0 = {tgt:4.1f}  <=>  theta_c = {thc:.6f} = {thc/np.pi:.5f} pi"
              f"  (cos(th_c/2) = {np.cos(thc/2):.5f}, 1-cos = {1-np.cos(thc/2):.5f})")
    print(f"  已测点落位: pi/2 -> beta_0 = {rows[4][3]:.5f}; 3pi/4 -> {rows[6][3]:.5f}"
          f" (均 <= 2.5 => 判 (a))")

    # ---- [P3] 直接核: [2.5,6] 上不等式逐点成立吗
    print("\n[P3] 直接核 (c=1/6): beta in {2.5,3,4,6} 时 Var(beta) <= V_inf 且 V_inf <= c*beta")
    print("  theta/pi   beta    Var(beta)   V_inf      c*beta    余量 c*beta-Var   余量/Var")
    worst = (1e9, None)
    for th in [np.pi / 2, 3 * np.pi / 4, 7 * np.pi / 8, 15 * np.pi / 16]:
        v, _, _, _ = vinf(th)
        for b in (2.5, 3.0, 4.0, 6.0):
            var = slice_stats(b, th)["var"]
            slack = C_SU2 * b - var
            print(f"  {th/np.pi:8.4f} {b:6.2f}  {var:10.6f}  {v:10.6f}  {C_SU2*b:9.5f}"
                  f"   {slack:12.6f}   {slack/var:9.4f}")
            if slack / var < worst[0]:
                worst = (slack / var, (th, b))
    print(f"  最紧一格: theta={worst[1][0]/np.pi:.4f}pi, beta={worst[1][1]},"
          f" 余量/Var = {worst[0]:.3f}")

    # ---- [P4] theta=pi 退化格与 H1A 的 beta* 对账
    print("\n[P4] 退化格 theta=pi: V_inf 发散 (=beta^2 V0), E[H]==0 => 策略第二假设本身失效")
    for b in (2.5, 6.0, 200.0):
        d = slice_stats(b, np.pi)
        print(f"  theta=pi beta={b:6.2f}  Var={d['var']:11.6f}  Var/beta^2={d['var']/b**2:.7f}"
              f"  E[H]={d['EH']:.2e}  E[Phi'']={d['EP2']:.7f}")
    print("  => beta_0(theta=pi) = +inf; 该格是零测退化切片 (H1A beta*=sqrt(0.2051153/0.0110593)"
          f" = {np.sqrt(0.2051153/0.0110593):.5f})")


if __name__ == "__main__":
    main()
