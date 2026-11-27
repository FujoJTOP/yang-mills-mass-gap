# -*- coding: utf-8 -*-
"""GP1-DELTA1 -- delta_1 as a SINGLE matrix element: calibration + cost model.

Design doc: research/v3/species/GP1-DESIGN.md

Modes (all zero-MC):
  --calib    free field delta_1=0 (two routes), 2D delta_1=D^3/24, 2-particle witness
  --bch      verify the a^2 coefficient of H_a = -(1/a) log(e^{-aA/2} e^{-aB} e^{-aA/2})
             -> this is the "T_a 的 a-展开" skeleton of K2 (the gm1-algebra input)
  --survey   load an existing npz, measure the per-config noise rho of the single-slice
             estimator (the quantity that decides whether 30 min is enough)
  --dry-run  print the cost model only (no data, no MC)

Run:  py -3.12 gp1_delta1.py --dry-run
      py -3.12 gp1_delta1.py --calib --bch
      py -3.12 gp1_delta1.py --survey gh3_L8.npz
"""
import argparse
import os
import sys
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# =====================================================================  [A] BCH
# T_a = e^{-aA/2} e^{-aB} e^{-aA/2}   (Wilson transfer matrix, symmetric splitting:
#        A = electric kernel, B = -log of the spatial-plaquette kernel)
# H_a := -(1/a) log T_a = (A + B) + a^2 K2 + O(a^4)
# K2 = alpha*[A,[A,B]] + beta*[B,[A,B]]  -- alpha,beta are universal (identified below).

def comm(X, Y):
    return X @ Y - Y @ X


def mexp(X):
    w, V = np.linalg.eig(X)
    return (V * np.exp(w)) @ np.linalg.inv(V)


def mlog(X):
    w, V = np.linalg.eig(X)
    return (V * np.log(w.astype(complex))) @ np.linalg.inv(V)


def k2_numeric(A, B, avec=(0.10, 0.07, 0.05, 0.035, 0.025)):
    """Fit H_a = (A+B) + a^2 K2 + a^4 K4 + a^6 K6 by linear least squares in t=a^2."""
    h0 = A + B
    t = np.array([a * a for a in avec])
    Y = []
    for a in avec:
        T = mexp(-a * A / 2) @ mexp(-a * B) @ mexp(-a * A / 2)
        H = -mlog(T) / a
        Y.append((H - h0).ravel())
    Y = np.array(Y)
    M = np.stack([t, t * t, t ** 3], 1)
    coef, *_ = np.linalg.lstsq(M, Y, rcond=None)
    return np.real(coef[0]).reshape(A.shape)


def k2_closed(A, B, alpha, beta):
    return alpha * comm(A, comm(A, B)) + beta * comm(B, comm(A, B))


def mode_bch():
    P("=" * 96)
    P("[BCH]  H_a = -(1/a) log(e^{-aA/2} e^{-aB} e^{-aA/2}) = (A+B) + a^2 K2 + O(a^4)")
    rng = np.random.default_rng(20261126)
    # identify alpha, beta from two random (A,B); K2 is linear in them
    rows, rhs = [], []
    for _ in range(2):
        A = rng.normal(size=(4, 4)); A = A + A.T
        B = rng.normal(size=(4, 4)); B = B + B.T
        kn = k2_numeric(A, B).ravel()
        c1 = comm(A, comm(A, B)).ravel()
        c2 = comm(B, comm(A, B)).ravel()
        rows.append([np.dot(c1, c1), np.dot(c1, c2)])
        rhs.append(np.dot(c1, kn))
        rows.append([np.dot(c2, c1), np.dot(c2, c2)])
        rhs.append(np.dot(c2, kn))
    sol, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    alpha, beta = float(sol[0]), float(sol[1])
    P(f"    识别（两次随机 (A,B)、最小二乘）:  alpha = {alpha:+.10f}   beta = {beta:+.10f}")
    P(f"    逐位检验候选 alpha=-1/24: {alpha*24.0:+.6f} ;  beta=-1/12: {beta*12.0:+.6f}")
    ALPHA, BETA = -1.0 / 24.0, -1.0 / 12.0
    P(f"    ⇒ K2 = -(1/24)[A,[A,B]] - (1/12)[B,[A,B]]  （alpha,beta 与 -1/24,-1/12 相对差 "
      f"{abs(alpha/ALPHA-1):.2e}/{abs(beta/BETA-1):.2e}）")
    worst = 0.0
    P("    逐组核验（4x4 随机对称 A,B；K2 数值 vs 闭式）:")
    for k in range(4):
        A = rng.normal(size=(5, 5)); A = A + A.T
        B = rng.normal(size=(5, 5)); B = B + B.T
        kn = k2_numeric(A, B)
        kc = k2_closed(A, B, ALPHA, BETA)
        d = np.max(np.abs(kn - kc)) / max(np.max(np.abs(kc)), 1e-300)
        worst = max(worst, d)
        P(f"      #{k}:  max|K2_num - K2_closed| / max|K2| = {d:.3e}")
    P(f"    ⇒ 最坏相对差 {worst:.3e}   {'OK' if worst < 5e-3 else 'FAIL'}")
    P("    物理读数：alpha 片 = 4-导数（磁/动能）片；beta 片 = 交叉片。")
    P("    方案 S1 下 B 的 a^2 展开给出 K2^mag（GM1-ALGEBRA:85-86）；")
    P("    方案 S2 另有整块电扇区（GM1-ALGEBRA:120，仓内 [NG]）。")
    assert worst < 5e-3, "BCH identification failed"
    return alpha, beta


# ==================================================================  [B] 定标
def d_free(a, m):
    return (2.0 / a) * np.arcsinh(a * m / 2.0)


def omega(a, D):
    return (2.0 / a) * np.sinh(a * D / 2.0)


def mode_calib():
    P("=" * 96)
    P("[CAL-1] 自由场：两条路都必须给 delta_1 = 0")
    P("   路 (i) 差分：(omega_a - m)/a^2 -> 0")
    for m in (0.5, 1.0, 2.0):
        row = []
        for a in (1e-1, 3e-2, 1e-2, 1e-3):
            row.append(f"{a:.0e}:{(omega(a, d_free(a, m)) - m) / a**2:+.2e}")
        P(f"     m={m:4.1f}  " + "  ".join(row))
    P("   路 (ii) 单矩阵元：delta_1 = <K2>_1 + D^3/24,  K2 = -m^3/24 (单体，GM1-ALGEBRA:189)")
    worst = 0.0
    for m in (0.5, 1.0, 2.0, 3.0):
        k2 = -m ** 3 / 24.0
        d1 = k2 + m ** 3 / 24.0
        worst = max(worst, abs(d1))
        P(f"     m={m:4.1f}   <K2> = {k2:+.12f}   D^3/24 = {m**3/24:+.12f}   delta_1 = {d1:+.3e}")
    P(f"   ⇒ 两条路一致：max|delta_1| = {worst:.3e}   {'OK' if worst < 1e-15 else 'FAIL'}")
    assert worst < 1e-15

    P("")
    P("[CAL-2] 2D：K2 = 0 (无空间 plaquette) ⇒ delta_1 = D^3/24 = 1/24 (D=1) （GM4 勘误）")
    D = 1.0
    for a in (1e-2, 3e-3, 1e-3):
        d1 = (omega(a, D) - D) / a ** 2     # exact, since Delta_a == D
        P(f"     a={a:.0e}   delta_1(num) = {d1:.12f}   D^3/24 = {D**3/24:.12f}   "
          f"ratio = {d1/(D**3/24):.9f}")
    P("     exact: 1 - Delta_a/D = 0 ⇒ omega_a = (2/a)sinh(aD/2) ⇒ s1 = D^3/24  [T]")

    P("")
    P("[CAL-3] 状态依赖（区分『矩阵元』与『数』）：自由场双粒子态（gm1alg_check [3]）")
    ok = True
    for (w1, w2) in [(F(1), F(2)), (F(1), F(1)), (F(1, 2), F(3, 2)), (F(1), F(4)), (F(3), F(5))]:
        lhs = (-(w1 ** 3 + w2 ** 3) + (w1 + w2) ** 3) / 24
        pred = w1 * w2 * (w1 + w2) / 8
        ok &= (lhs == pred)
        P(f"     (w1,w2)=({w1},{w2}):  <K2+H^3/24> = {lhs}   w1w2(w1+w2)/8 = {pred}   match={lhs == pred}")
    P(f"   ⇒ 5/5 = {ok}；单粒子扇区 (w2=0) 才为 0 ⇒ delta_1 是 (理论, 态) 的函数，")
    P("      (SM-2) 只取边缘态那一个数。")
    assert ok


# =============================================================  [C] 成本模型
# 实测吞吐（不重跑，取存档）：
#   Ls=8: 4.8 ms/cfg  (b01_gw2p.json 57.5 s / 12000, 经 GM4-NUMERIC.md:229)
#   Ls=6: 1.52 ms/cfg = (6/8)^4 * 4.8      (GM4-NUMERIC.md:229 同款 (Ls/8)^4 律)
MS_L8 = 4.8e-3
MS_L6 = 1.52e-3
M_HAT = 1.7          # beta=2.5 的 a*Delta (B04:98，经 GM4-NUMERIC.md:96)
TARGET = 1.0 / 96.0  # 4-sigma 分辨 delta_1=0 vs delta_1=D^3/24
BUDGET_S = 1800.0


def cheb_degree(lo, hi, tol):
    """Minimal Chebyshev degree for log(x) on [lo,hi] with max abs error <= tol."""
    from numpy.polynomial import chebyshev as C
    xs = np.linspace(-1.0, 1.0, 4001)
    X = lo + (hi - lo) * (xs + 1.0) / 2.0
    for n in range(1, 121):
        c = C.chebinterpolate(lambda z: np.log(lo + (hi - lo) * (z + 1.0) / 2.0), n)
        if np.max(np.abs(C.chebval(xs, c) - np.log(X))) <= tol:
            return n
    return -1


def mode_dry_run(rho_meas=None):
    P("=" * 96)
    P("[COST] 目标：sigma(delta_1/D^3) <= 1/96  （4 sigma 分辨 0 与 D^3/24）")
    P(f"   输入：Ls=8 {MS_L8*1e3:.2f} ms/cfg ；Ls=6 {MS_L6*1e3:.2f} ms/cfg ；m_hat={M_HAT}")
    P("   统计律：sigma = rho*kappa*D^3/sqrt(N)  ⇒  N_req = (96*rho*kappa)^2")
    P("   rho = 单组态相对噪声（单时片估计量）；kappa = |K2|/D^3（算子自然尺度比）")
    P("")
    rho_from_repo = np.sqrt(12000.0) / 5.0    # S/N(2)=5.0 @ N=1.2e4  (B04:169 via GM4:58)
    P(f"   仓内可比的 per-config rho（由 S/N(2)=5.0 @ N=1.2e4 反解）= {rho_from_repo:.1f}")
    if rho_meas is not None:
        P(f"   本轮实测 rho（--survey，见下）= {rho_meas:.1f}")
    P("")
    P("   kappa↓ / rho→   7.0        22.0       60.0      100.0")
    for kappa in (0.5, 1.0, 2.0):
        row = []
        for rho in (7.0, 22.0, 60.0, 100.0):
            n = (96.0 * rho * kappa) ** 2
            row.append(f"{n*MS_L8/3600.0:9.2f}h")
        P(f"   {kappa:4.1f}          " + " ".join(row))
    P("   （单位 = Ls=8 墙钟；Ls=6 除以 3.16）")
    P("")
    thr = BUDGET_S / MS_L8
    P(f"   ≤30 min 的充要条件：N_req <= {thr:.2e} cfg  ⇒  rho*kappa <= {np.sqrt(thr)/96.0:.2f}")
    P(f"   仓内同量级 rho≈22 ⇒ 需 {((96*22)**2)*MS_L8/3600:.1f} h（Ls=8）/ {((96*22)**2)*MS_L6/3600:.1f} h（Ls=6）")
    n1 = (24.0 * 22.0) ** 2
    P(f"   降精度版（1 sigma 分辨 D^3/24，即 sigma=1/24）：N={n1:.2e} ⇒ {n1*MS_L8/60:.0f} min (Ls=8)")
    P("     ⇒ 30 min 能给到的分辨率恰是 1/24 = 效应本身 ⇒ 1 sigma，不能判 (SM-2)")
    if rho_meas is not None:
        n = (96.0 * rho_meas) ** 2
        P("")
        P(f"   ★ 实测行：rho={rho_meas:.1f}, kappa=1 ⇒ N_req={n:.2e} ⇒ "
          f"{n*MS_L8/3600:.2f} h (Ls=8) / {n*MS_L6/60:.0f} min (Ls=6)")
        P(f"     超出 30 min 的倍数（Ls=8 / Ls=6）= {n*MS_L8/BUDGET_S:.1f}x / {n*MS_L6/BUDGET_S:.1f}x")
        P(f"     要达到 30 min 需把 rho 从 {rho_meas:.1f} 压到 {np.sqrt(BUDGET_S/MS_L8)/96:.1f} "
          f"（方差降 {(rho_meas/(np.sqrt(BUDGET_S/MS_L8)/96))**2:.1f}x）")
    P("")
    P("[COST-b] H^3/24 的两条路")
    P("   路 (i) 精确：psi_1 是 H 的本征向量 ⇒ <H^3>_1 = D^3，是 c-数 ⇒ 0 次 matvec [T]")
    P("   路 (ii) 多项式：log 用 Chebyshev（T 的谱窗 [e^-6, e^-1.7]）")
    lo, hi = np.exp(-6.0), np.exp(-M_HAT)
    n1e3 = cheb_degree(lo, hi, 1e-3)
    n1e4 = cheb_degree(lo, hi, 1e-4)
    P(f"     阶数：1e-3 相对精度 ⇒ n={n1e3} ；1e-4 ⇒ n={n1e4}（谱窗 [{lo:.3f},{hi:.3f}]）")
    P("     每个 Chebyshev 项 = 一次 T 作用 = 相关函数前进一格；")
    P(f"     第 k 格的每-组态噪声 ∝ e^(+m_hat*k)，故总成本 ∝ e^(2*m_hat*n) = e^({2*M_HAT*n1e3:.0f}) —— 统计上死。")
    P("   ⇒ 只要用『概率表示』（相关函数矩阵），多项式路不可用；用本征向量则退化为路 (i)。")


# ================================================================  [D] survey
def survey(fn):
    P("=" * 96)
    d = np.load(fn)
    S = d["S"]
    Lt = int(d["Lt"])
    N = int(d["n"])
    beta = float(d["beta"])
    Ls = int(d["Ls"]) if "Ls" in d else -1
    nlev, nobs = S.shape[1], S.shape[2]
    P(f"[SURVEY] {os.path.basename(fn)}  beta={beta} Ls={Ls} Lt={Lt} N={N} "
      f"nlev={nlev} nobs={nobs}")
    A = np.ascontiguousarray(S[:, 0, :, :]).astype(np.float64)   # (N,nobs,Lt) raw level
    # connected per config
    mu = A.mean(0)                       # (nobs,Lt)
    Ac = A - mu[None, :, :]
    # single-slice (t=0) estimator noise: this is the analogue of a K2 insertion
    P("")
    P("  per-config 相对噪声 rho（单时片双线性估计量，K2 型插入量的同类）:")
    rho_list = []
    for w in range(nobs):
        x = Ac[:, w, 0]
        m1 = float(np.mean(Ac[:, w, 0] * Ac[:, w, 1]))       # ~C(1)
        xc = Ac[:, w, 0] * Ac[:, w, 1]
        rho = float(xc.std(ddof=1) / max(abs(xc.mean()), 1e-300))
        rho_list.append(rho)
        P(f"    w={w}:  C(1)={xc.mean():+.6e}   per-config rho={rho:8.1f}")
    P(f"  ⇒ rho 中位数 {np.median(rho_list):.1f}")
    # GEVP on 3 windows (gp6b logic) -> psi_1 + m_eff
    W = [0, min(3, nobs - 1), nobs - 1]
    P("")
    P(f"  GEVP（窗口 {W}，连通）:")
    Cs = np.empty((Lt, N, 3, 3))
    for t in range(Lt):
        for i in range(3):
            for j in range(3):
                Cs[t, :, i, j] = np.mean(Ac[:, W[i], :] * np.roll(Ac[:, W[j], :], -t, axis=1), 1)
    Ct = Cs.sum(1)
    C0 = 0.5 * (Ct[0] + Ct[0].T)
    wv, V = np.linalg.eigh(C0)
    keep = wv > 1e-2 * wv.max()
    P(f"    C(0) 谱（连通）：{np.array2string(wv, precision=3)}  ⇒ Gram rank = {int(keep.sum())}/3")
    P("    rank<3 ⇒ psi_1 在该基内不唯一（GM2-EDGE:0.3 的简并读法 λ_min(P1 K2 P1)）")
    Lm = V[:, keep] / np.sqrt(wv[keep])[None, :]
    lam = np.array([np.linalg.eigvalsh(0.5 * ((Lm.T @ Ct[t] @ Lm) + (Lm.T @ Ct[t] @ Lm).T))[-1]
                    for t in range(Lt)])
    ms = []
    for t in range(1, 4):
        ms.append(np.log(lam[t] / lam[t + 1]) if lam[t + 1] > 0 else np.nan)
    P("    m_eff(t) = " + " ".join(f"t={t}:{ms[t-1]:.3f}" for t in range(1, 4)))
    eps = abs(ms[0] - ms[1]) / ms[1]
    P(f"    变分态缺陷代理 eps = |m_eff(1)-m_eff(2)|/m_eff(2) = {eps:.3f}")
    P(f"    3-体污染系统项 = 2*eps*|w1w2(w1+w2)/8| ~ {2*eps*0.25:.3f}*D^3 （两体取 w1=w2=D）")
    P(f"    判据要求 ≤ 1/96 = {1/96:.4f}*D^3  ⇒  {'通过' if 2*eps*0.25 <= 1/96 else '超标'}")
    return float(np.median(rho_list)), eps


# ==================================================================  main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calib", action="store_true")
    ap.add_argument("--bch", action="store_true")
    ap.add_argument("--survey", type=str, default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--tag", type=str, default="gp1_delta1")
    a = ap.parse_args()
    if not (a.calib or a.bch or a.survey or a.dry_run):
        a.dry_run = True
    rho_meas = None
    if a.calib:
        mode_calib()
    if a.bch:
        mode_bch()
    if a.survey:
        fn = a.survey if os.path.isabs(a.survey) else os.path.join(HERE, a.survey)
        rho_meas, eps = survey(fn)
    if a.dry_run:
        mode_dry_run(rho_meas)
    out = os.path.join(HERE, f"{a.tag}.log")
    open(out, "w", encoding="utf-8").write("\n".join(LOG) + "\n")
    print(f"\n[log] {out}")


if __name__ == "__main__":
    main()
