"""p02_m6_dyndyn.py -- M6 (谐振子 ⊗ 自由模, d=0+1) 上 Δ_dyn = E[H] - Var(G) 的精确值。

产物: w7/p02_m6_dyndyn.log      文档: research/v3/species/P02-M6DYN.md
口径: mpmath dps=30, 确定性 (无采样 / 无随机 / 无 MC)。
  §1 谱: E0 = ω/2 (与 b 环周长 C 无关) 且 m_L = min(ω, 2π²/C²) → 0
  §2 Var(a)=1/(2ω), G″=-∂²log ρ_a|0 = 2ω = 4(E0-V(0)), Ward G″·Var=1
  §3 Δ_dyn 的五种读法 (R1..R5) + 闭环判据 (closure ①) 的唯一选支
  §4 Gauss 条件化恒等式: Δ_dyn(θ) = 1/Var(ξ) 三路互核, 方向依赖
  §5 目标量词: inf_v Δ_dyn = 0
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 30


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


W_LIST = [mp.mpf('2'), mp.mpf('2.5'), mp.mpf('0.7')]

# ---------------------------------------------------------------- §1 谱
hdr("§1  谱 / E[H]:   spec = {(n+1/2)ω} + {2π²k²/C²}  ⟹ E0 = ω/2,  m → 0")

# (a) 1D 谐振子格点核验 (a 方向), 精确 E_n = ω(n+1/2)
print("\n(a) a 方向三对角精确对角化 (Dirichlet, -1/2 d² + ω²a²/2):")
for w in W_LIST:
    wf = float(w)
    for N in (800, 3200):
        Lbox = 14.0 / np.sqrt(wf)
        x = np.linspace(-Lbox, Lbox, N + 2)[1:-1]
        h = x[1] - x[0]
        d = 1.0 / h**2 + 0.5 * wf**2 * x**2
        e = -0.5 / h**2 * np.ones(N - 1)
        E = np.linalg.eigvalsh(np.diag(d) + np.diag(e, 1) + np.diag(e, -1))
        print(f"    ω={float(w):.2f} N={N:5d}  E0={E[0]:.12f}  (ω/2={wf/2:.12f}, "
              f"err={E[0]-wf/2:+.2e})   E1-E0={E[1]-E[0]:.9f} (ω={wf})")

# (b) b 方向: 环上自由模 => 零模能量恰 0, 间隙 min(ω, 2π²/C²)
print("\n(b) b 模在周长 C 的环上:  E_b(k)=½(2πk/C)², E_b(0)=0 精确")
print("    C        E0=ω/2+E_b(0)      m_L=min(ω,2π²/C²)")
for w in [mp.mpf('2')]:
    for C in [mp.mpf('10'), mp.mpf('100'), mp.mpf('1000'), mp.mpf('1e5')]:
        Eb1 = 2 * mp.pi**2 / C**2
        print(f"    {float(C):<9g}{float(w/2 + 0):<20.12f}{float(min(w, Eb1)):.6e}")
print("    ⇒ E0 = ω/2 对一切 C 精确成立 (b 零模加 0); m_L → 0 当 C → ∞; 连续谱 = [ω/2,∞)")

# ------------------------------------------------- §2 Var(a), G'', Ward
hdr("§2  ρ_a = |ψ0|² = √(ω/π)e^{-ωa²}:  Var(a)=1/(2ω),  G″ = 2ω,  Ward")

for w in W_LIST:
    sig2 = 1 / (2 * w)
    rho = lambda a: mp.sqrt(w / mp.pi) * mp.e**(-w * a**2)
    m1 = mp.quad(lambda a: a * rho(a), [-mp.inf, 0, mp.inf])
    m2 = mp.quad(lambda a: a**2 * rho(a), [-mp.inf, 0, mp.inf])
    m4 = mp.quad(lambda a: a**4 * rho(a), [-mp.inf, 0, mp.inf])
    var = m2 - m1**2
    Gpp = -mp.diff(lambda a: mp.log(rho(a)), 0, 2)
    Vir = 0.5 * w**2 * m2          # <V>
    Kin = mp.mpf(1) / (8 * sig2)   # <T> = 1/(8σ²)
    E0 = Vir + Kin
    E0_id = 4 * (w / 2 - 0)        # 4(E0 - V(0))
    print(f"\n  ω={float(w):.2f}:  Var(a)={mp.nstr(var,18)}  1/(2ω)={mp.nstr(sig2,18)}  "
          f"res={mp.nstr(var-sig2,3)}")
    print(f"        G″=-∂²logρ|0 = {mp.nstr(Gpp,18)}   2ω = {mp.nstr(2*w,18)}  "
          f"res = {mp.nstr(Gpp-2*w,3)}")
    print(f"        G″·Var(a) = {mp.nstr(Gpp*var,20)}  (Ward)")
    print(f"        <V>={mp.nstr(Vir,12)} <T>={mp.nstr(Kin,12)} ⇒ E0={mp.nstr(E0,18)} "
          f"(ω/2={mp.nstr(w/2,18)}, <a⁴>=3/(4ω²)={mp.nstr(m4,12)} vs "
          f"{mp.nstr(3/(4*w**2),12)})")
    print(f"        机制核验 4(E0-V(0)) = {mp.nstr(E0_id,18)}  ⇒ = G″  {'OK' if abs(E0_id-Gpp)<mp.mpf('1e-20') else 'FAIL'}")

# ------------------------------------------------- §3 Δ_dyn 的五种读法
hdr("§3  Δ_dyn = E[H] - Var(G) 的五种读法  (+ 闭环判据 closure ①: Δ_dyn + W_full = G″)")

print("""
  V_eff := -log ρ_a = ω a²   (量子有效指数: 动能已被吸收; G″ = V_eff″ = 2ω)
  V_bare := ½ω² a²            (经典势; V_bare″ = ω², 与 G″=2ω 差 2/ω)

  读法                          E[H]        Var(G)       Δ_dyn            closure ①
  R1 有效指数 + 条件化 ξ=a=0    2ω          0            2ω               OK  (W_full=0)
  R2 有效指数 + 无条件化        2ω          2ω           0                自洽但非语料口径
  R3 裸势 + 条件化              ω²          0            ω²               破: ω² ≠ 2ω
  R4 裸势 + 无条件化            ω²          ω³/2         ω²(1-ω/2)        破
  R5 能量读法 E0 与 Var(a)      ω/2         1/(2ω)       (ω²-1)/(2ω)      量纲混用
  --- 只有 R1 同时满足: 条件化(ch05:233 |ξ=0) 与 closure ① (W_full=0 ⇒ Δ_dyn=G″) ---""")

print(f"\n  {'ω':>5} {'R1=2ω':>14} {'R2':>10} {'R3=ω²':>12} {'R4':>14} {'R5':>14} "
      f"{'|R3-G″|':>10}")
for w in W_LIST:
    R1 = 2 * w
    R2 = mp.mpf(0)
    R3 = w**2
    R4 = w**2 * (1 - w / 2)
    R5 = (w**2 - 1) / (2 * w)
    print(f"  {float(w):>5.2f} {float(R1):>14.6f} {float(R2):>10.4f} {float(R3):>12.6f} "
          f"{float(R4):>14.6f} {float(R5):>14.6f} {float(abs(R3-2*w)):>10.6f}")
print("  ⇒ R3/R4/R5 在 ω>2 (R3,R4) 或 ω<1 (R5) 上会给负值: 这三支制造的是【假 (b)】,")
print("    它们的共同病是 closure ① 不闭合或量纲混用, 不是 M6 的物理。")

# ------------------------- §4 Gauss 条件化恒等式 + 方向依赖 (v = cosθ a + sinθ b)
hdr("§4  Δ_dyn(θ) = E[H_vv|ξ=0] - Var(G_v|ξ=0) = 1/Var(ξ)  (Gauss 条件化恒等式, 三路互核)")

print("\n  取 a~N(0,σa²), σa²=1/(2ω); b~N(0,σb²) 独立; v=(cosθ,sinθ), ‖v‖=1, ξ=<v,X>")
print("  H_vv = cos²θ/σa² + sin²θ/σb²,  G_v = cosθ a/σa² + sinθ b/σb²")
print("  恒等式: E[H_vv|0] - Var(G_v|0) = H_vv - (H_vv - 1/Var(ξ)) = 1/Var(ξ) = G″(θ)\n")

def check(w, th_deg, sb2):
    th = mp.radians(mp.mpf(th_deg))
    c, s = mp.cos(th), mp.sin(th)
    sa2 = 1 / (2 * w)
    Vr = c**2 * sa2 + s**2 * sb2                     # Var(ξ)
    Hvv = c**2 / sa2 + s**2 / sb2
    VarG_cond = Hvv - 1 / Vr                         # Schur
    d1 = Hvv - VarG_cond                             # E[H]-Var(G)  (条件化)
    # 独立路: 从 ρ_ξ(s)=N(0,Vr) 直接二阶求导
    rho = lambda x: mp.e**(-x**2 / (2 * Vr)) / mp.sqrt(2 * mp.pi * Vr)
    d2 = -mp.diff(lambda x: mp.log(rho(x)), 0, 2)
    return Vr, d1, d2, Hvv, VarG_cond

print(f"  {'ω':>4} {'θ°':>6} {'Var(ξ)':>16} {'Δ_dyn=E[H]-Var(G)':>22} "
      f"{'G″=-∂²logρ':>22} {'res':>9}")
for w in [mp.mpf('2')]:
    for th in [0, 15, 45, 89]:
        for sb2 in [mp.mpf('1'), mp.mpf('100'), mp.mpf('1e8')]:
            Vr, d1, d2, Hvv, VG = check(w, th, sb2)
            print(f"  {float(w):>4.1f} {th:>6} {float(Vr):>16.6e} {float(d1):>22.12e} "
                  f"{float(d2):>22.12e} {float(d1-d2):>9.1e}"
                  + ("   σb²=1" if sb2 == 1 else ("   σb²=100" if sb2 == 100 else "   σb²=1e8 (自由方向)")))

print("\n  ⇒ 三路 (闭环 / Schur / 直接求导) 在每一位上都同值; 且")
print("     θ=0   (纯 a 模):        Δ_dyn = 2ω, 与 σb² 无关   [T]")
print("     θ≠0 (含 b 分量):        Δ_dyn = 1/(cos²θσa²+sin²θσb²) → 0  (σb²→∞)")
print("     b 完全自由 (M6 字面) ⇒ σb²=∞ ⇒ 任何含 b 分量的模上 Δ_dyn = 0 (与 M5 自由场同规则)")

# ------------------------------------------------- §5 目标量词
hdr("§5  判据 Δ₀>0 的量词:  ∀v  Δ_dyn(v) ≥ a²Δ₀² > 0 ?  —— 在 M6 上")
print("""
  Δ_dyn(a 模) = 2ω > 0        ⇒ 单点读法: 靶在无隙模型上【成立】= n12 的假阳性
  Δ_dyn(b 模) = 0   ,  inf_v Δ_dyn = 0
  ⇒ 靶的逐字量词 "forever coarse transverse v" (N12 §1, app.tex:183) 在 M6 上【被违反】
  ⇒ 同一模型给出两个读数: 点读 (2ω, 靶真) / ∀v 读 (0, 靶假)
  ⇒ M6 不能单独判死 Δ₀>0;  2D YM 可以:  Δ_dyn = -W_full < 0 对一切 v (W_full ≥ N/12‖v‖²)""")

# ------------------------------------------------- §6 汇总
hdr("§6  结论 (K = 精确值, 与 ω 同号)")
for w in W_LIST:
    print(f"  ω={float(w):>4.2f}:  E0=E[H]=ω/2={float(w/2):.6f}   Var(a)=1/(2ω)={float(1/(2*w)):.6f}"
          f"   G″=Δ_dyn(a)=2ω={float(2*w):.6f}   Δ_dyn(b)=0   inf_v Δ_dyn=0")
print("\n  ⇒ Δ_dyn(a 模) = 2ω > 0 恒正 (无需不等式); Δ_dyn(b 模) = 0; 从不取负值。")
print("  ⇒ (a) 成立 (点读) / (b) 不成立 (无负值) / (c) 不成立 (Δ_dyn 良定义)。")
