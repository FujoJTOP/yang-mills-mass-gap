# HQ7 · 柱 2「P 的匹配时刻读数」零 MC 复核（β=2.5 的 L-序列）
# 输入：research/v3/species/w7/gh2_b2.5_Ls{6,8}_Lt{12,20}.npz（既有配置，不产生新 MC）
# 读出：rho(τ) = C(τ)/C(0)，在匹配时刻 τ = L_s 上取值（P 的 LHS 形状：G_c(L)=e^{√c}C_L(L)）
# 纪律：τ 必须 <= L_t/2 才免周期时间回绕污染；两条 L_t 都算，用来显式暴露回绕。
import json, math
import numpy as np

W = r"research/v3/species/w7"


def corr(S, Lt):
    S = np.asarray(S, float)
    S = S - S.mean()
    return np.array([(S * np.roll(S, -t, axis=1)).mean(axis=0).mean() for t in range(Lt)])


def jk(S, Lt, nblk=20):
    S = np.asarray(S, float)
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
    full = corr(S, Lt)
    err = np.sqrt((nblk - 1) / nblk * ((vals - vals.mean(0)) ** 2).sum(0))
    return full, err


def run(fn):
    d = np.load(f"{W}/{fn}", allow_pickle=True)
    S = d["S"][:, 0, :]
    Lt = int(d["Lt"]); Ls = int(d["Ls"])
    C, sC = jk(S, Lt)
    rho = C / C[0]
    srho = sC / C[0]
    pure = Lt // 2
    print(f"\n[{fn}] beta={float(d['beta'])} Ls={Ls} Lt={Lt} N={S.shape[0]} clean_tau<={pure}")
    print("   t   rho(t)        +-          S/N")
    for t in range(Lt):
        print(f"  {t:2d}  {rho[t]: .6e}  {srho[t]:.2e}  {abs(rho[t])/srho[t] if srho[t]>0 else float('inf'):6.2f}")
    t = Ls
    ok = t <= pure
    print(f"  >> 匹配时刻 tau=Ls={t}: rho={rho[t]:.4e} +- {srho[t]:.2e} "
          f"(S/N={abs(rho[t])/srho[t] if srho[t]>0 else float('inf'):.2f}) "
          f"{'CLEAN' if ok else 'WRAPPED (tau>Lt/2) -> 作废'}")
    return dict(fn=fn, Ls=Ls, Lt=Lt, rho=rho.tolist(), srho=srho.tolist(),
                match_t=t, match_rho=float(rho[t]), match_srho=float(srho[t]), clean=bool(ok))


if __name__ == "__main__":
    out = []
    for fn in ["gh2_b2.5_Ls6_Lt20.npz", "gh2_b2.5_Ls8_Lt20.npz",
               "gh2_b2.5_Ls6_Lt12.npz", "gh2_b2.5_Ls8_Lt12.npz"]:
        out.append(run(fn))

    print("\n=== P 的 L-序列读数（只取 CLEAN 的两格，Lt=20）===")
    clean = [o for o in out if o["clean"]]
    for o in sorted(clean, key=lambda x: x["Ls"]):
        print(f"  L={o['Ls']}: rho(L)={o['match_rho']:.4e} +- {o['match_srho']:.2e}")
    if len(clean) == 2:
        a, b = sorted(clean, key=lambda x: x["Ls"])
        d = b["match_rho"] - a["match_rho"]
        sd = math.hypot(a["match_srho"], b["match_srho"])
        print(f"  Delta = rho(8)-rho(6) = {d:.3e} +- {sd:.2e}  z={d/sd if sd>0 else float('nan'):+.2f}")
        print("  判读：z<0 且 |d|>2sd => 匹配时刻读数随 L 下降（与无平台相容）；"
              "z>0 且 |d|>2sd => 上升（平台征兆）；否则无分辨率")
    # 用实测 m_eff 估 mL，判断"匹配时刻是否已被体积锁死"
    for o in clean:
        r = o["rho"]
        if r[1] > 0 and r[0] > 0:
            meff0 = math.log(r[0] / r[1])
            print(f"  L={o['Ls']}: 实测 m_eff(0)={meff0:.3f}（格点单位）=> m*L = {meff0*o['Ls']:.2f}")
    with open(r"research/v3/species/hq7_seq_check.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n[done]")
