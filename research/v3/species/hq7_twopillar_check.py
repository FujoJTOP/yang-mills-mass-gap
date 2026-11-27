# HQ7 · 柱 2「P 的数据侧检验」零 MC 复核
# 输入：仓内既有 npz（research/v3/species/w7/），不产生新 MC，不写 Lean。
# 检验：完全单调性两点式  C(2t)/C(0) >= (C(t)/C(0))^2   （log C 凸 / Cauchy--Schwarz）
# 另附：二进位传播  rho -> rho^2（HZ1 §2.2 的下一步），看它能否支撑/否证 P。
import json, math, sys
import numpy as np

W = r"research/v3/species/w7"


def corr_periodic(S, Lt):
    """C(t) = (1/Lt) sum_{t0} <S(t0+t) S(t0)>_c, 周期时间。S: (N, Lt)"""
    S = np.asarray(S, dtype=float)
    S = S - S.mean()
    N = S.shape[0]
    out = np.empty(Lt)
    for t in range(Lt):
        prod = (S * np.roll(S, -t, axis=1)).mean(axis=0)  # per t0
        out[t] = prod.mean()
    return out


def jk_err(S, Lt, nblk=20):
    """delete-one-block jackknife on C(t)。"""
    S = np.asarray(S, dtype=float)
    N = S.shape[0]
    S = S - S.mean()
    idx = np.array_split(np.arange(N), nblk)
    vals = np.empty((nblk, Lt))
    for b, ii in enumerate(idx):
        mask = np.ones(N, bool)
        mask[ii] = False
        Sb = S[mask]
        for t in range(Lt):
            vals[b, t] = (Sb * np.roll(Sb, -t, axis=1)).mean(axis=0).mean()
    full = np.array([(S * np.roll(S, -t, axis=1)).mean(axis=0).mean() for t in range(Lt)])
    jk = np.sqrt((nblk - 1) / nblk * ((vals - vals.mean(0)) ** 2).sum(0))
    return full, jk


def report(tag, C, sC, tmax_note):
    C0 = C[0]
    rho = C / C0
    srho = sC / C0
    print(f"\n=== {tag} : C(t)/C(0)  ({tmax_note}) ===")
    print("  t   rho(t)      +-        S/N     rho(t)^2     gap=rho(2t)-rho(t)^2")
    for t in range(len(C)):
        r2 = rho[t] ** 2 if t + 1 < len(C) else float("nan")
        gap = rho[2 * t] - r2 if 2 * t < len(C) else float("nan")
        snr = abs(rho[t]) / srho[t] if srho[t] > 0 else float("inf")
        print(f"  {t:2d}  {rho[t]: .6e}  {srho[t]:.3e}  {snr:6.2f}  {r2: .6e}   {gap: .4e}")
    print("  --- 两点式判读（t -> 2t）---")
    rows = []
    for t in range(1, len(C)):
        if 2 * t >= len(C):
            break
        lhs = rho[2 * t]
        rhs = rho[t] ** 2
        sl = srho[2 * t]
        sr = 2 * abs(rho[t]) * srho[t]
        d = lhs - rhs
        sd = math.hypot(sl, sr)
        verdict = "PASS" if d >= 0 else ("FAIL(noise)" if abs(d) < 2 * sd else "FAIL")
        print(f"  t={t}: rho({2*t})={lhs:.6e}+-{sl:.2e} ; rho({t})^2={rhs:.6e}+-{sr:.2e}"
              f" ; d={d:.3e}+-{sd:.2e} ; z={d/sd if sd>0 else float('nan'):+.2f} -> {verdict}")
        rows.append(dict(t=t, lhs=lhs, rhs=rhs, d=d, sd=sd,
                         z=(d / sd if sd > 0 else None), verdict=verdict))
    return rho, srho, rows


def propagate(rho0, t0, horizons):
    print(f"  二进位传播（从 rho({t0})={rho0:.6e} 出发，rho -> rho^2）：")
    r = rho0
    for h in horizons:
        print(f"    t={h:3d} : rho >= {r:.6e}")
        r = r * r
    return r


def do_npzgp3(fn, tag):
    d = np.load(f"{W}/{fn}", allow_pickle=True)
    C = d["C"]
    sC = d["sigC"]
    print(f"\n########## {fn}  beta={float(d['beta'])} Ls={int(d['Ls'])} Lt={int(d['Lt'])} N={int(d['N'])}")
    rho, srho, rows = report(tag, C, sC, "periodic Lt, published C/sigC")
    return rho, srho, rows


def do_b01(fn):
    d = np.load(f"{W}/{fn}", allow_pickle=True)
    S = d["S"][:, 0, :]
    Lt = int(d["Lt"])
    C, sC = jk_err(S, Lt)
    print(f"\n########## {fn}  beta={float(d['beta'])} Lt={Lt} N={S.shape[0]}  (重算, jackknife 20 blk)")
    rho, srho, rows = report(f"b01 beta={float(d['beta'])}", C, sC, "jackknife 重算")
    return rho, srho, rows


if __name__ == "__main__":
    out = {}
    for fn in ["gp3_b25_Ls8_Lt16_N24k.npz", "gp3_b3_Ls8_Lt16_N100k.npz"]:
        rho, srho, rows = do_npzgp3(fn, fn)
        out[fn] = dict(rho=rho.tolist(), srho=srho.tolist(), rows=rows)
        # 传播：取"最大且 S/N>=3"的 t
        ok = [t for t in range(len(rho)) if srho[t] > 0 and abs(rho[t]) / srho[t] >= 3]
        if ok:
            t0 = max(ok)
            propagate(rho[0] if t0 == 0 else rho[t0], t0, [t0 * 2, t0 * 4, t0 * 8, t0 * 16])

    for fn in ["b01_b2.2.npz", "b01_b2.75.npz", "b01_b3.npz", "b01_b4.npz", "b01_b6.npz", "b01_b10.npz"]:
        rho, srho, rows = do_b01(fn)
        out[fn] = dict(rho=rho.tolist(), srho=srho.tolist(), rows=rows)

    # [w1] 原子地板的数据侧证人核：b01 的 w1 列（w1=f1^2/M0 必在 [0,1]）
    print("\n########## [w1] b01 的 w1 列（原子地板数据侧证人核；越界即不可用）")
    bad = []
    for fn in ["b01_b2.2.npz", "b01_b2.75.npz", "b01_b3.npz", "b01_b4.npz", "b01_b6.npz", "b01_b10.npz"]:
        d = np.load(f"{W}/{fn}", allow_pickle=True)
        w1 = float(d["w1"]); e1 = float(d["w1_err"])
        flag = "OK" if 0.0 <= w1 <= 1.0 else "OUT-OF-RANGE（w1∉[0,1]）"
        if flag != "OK":
            bad.append((fn, w1))
        print(f"  beta={float(d['beta']):5.2f}  w1={w1: .6g}  w1_err={e1:.3g}  -> {flag}")
    print(f"  越界点：{bad}；GP3-PLATEAU.md:154 已登记 w1_err 的轴 bug 产物（4.15e-15）"
          f" ⇒ w1 列不可用作地板证人")

    with open(r"research/v3/species/hq7_twopillar_check.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n[done] json written")
