# -*- coding: utf-8 -*-
"""
GAP4 — 熵 vs 能量：失效点 beta_c = 8/lambda_* 的定位，以及 [1/48, 2.5] (=120x) 的分解.

口径 (不许混用):
  N=2, a=1, d=4,  u = c_f/N = beta/(2N^2) = beta/8          [T] X01:15 / BETA-CONVENTION:64
  Delta = 20 = 过定 plaquette 的相邻 plaquette 数 = 4(2(d-1)-1)   [T] d=4 计数
  lambda_* = 4D plaquette 动物生长常数: a_n = 1,20,458,11132,280175 (n=1..5)
             相邻比 20, 22.9, 24.3, 25.2 (仍在升)             [E] C03-BYPASS:111
             B_n 比  20, 25.70, 28.69, 30.53, 31.77         [E] UR4:14
             包络 lambda_T = e*Delta = 54.3656365692        [E] P4:153
  判据: 活动度 x 熵 = lambda_* u ;  = 1  <=>  beta = 8/lambda_*
        Penrose/KP 型充分判据:  K u < 1,  K = 2N^2 = 8,  K u = beta  [T] TN2:88 / ADJ2:56

无 git commit；只写 research/v3/species/GAP4.md + research/v3/w7/gap4_*.py|.log
"""
import mpmath as mp

mp.mp.dps = 40
P = print

N, d = 2, 4
DELTA = 4 * (2 * (d - 1) - 1)          # = 20
KHAT = 2 * N * N                        # = 8
B48 = mp.mpf(1) / 48
B25 = mp.mpf("2.5")
SPAN = B25 / B48                        # 120 exactly

A = [1, 20, 458, 11132, 280175]         # C03:111, n=1..5
RATIO_OBS = [mp.mpf(20), mp.mpf("22.9"), mp.mpf("24.3"), mp.mpf("25.2")]
B_RATIO = [mp.mpf(20), mp.mpf("25.70"), mp.mpf("28.69"), mp.mpf("30.53"), mp.mpf("31.7727")]

LAM_A = mp.mpf(A[4]) / A[3]             # a5/a4 = 280175/11132 (精确比)
LAM_B = B_RATIO[4]                      # 31.7727
LAM_ENV = mp.e * DELTA                  # 54.3656365692

P("=" * 78)
P("GAP4  entropy-vs-energy  beta_c = 8/lambda_*   [mpmath dps=40]")
P("=" * 78)
P(f"  N={N} d={d} Delta={DELTA} K_hat=2N^2={KHAT}  u=beta/8  1/48={mp.nstr(B48,12)}")
P(f"  span 2.5/(1/48) = {mp.nstr(SPAN,15)}   (=120 exact)")

# ---------------------------------------------------------------- 0  sanity
P("\n[0] 口径自检")
P(f"  K_hat*u(beta) = 2N^2*beta/(2N^2) = beta   (N-无关)                     [T]")
for b in ("0.018394", "0.25", "2.5"):
    bb = mp.mpf(b)
    P(f"    beta={mp.nstr(bb,6):>8}  K_hat*u={mp.nstr(bb,10):>12}  u={mp.nstr(bb/8,10):>12}")

# ------------------------------------------------- 1  lambda_* 的三个读数
P("\n[1] lambda_* 的三个读数（含 root test 与 ratio test 的差别）")
P("  n   a_n        a_n^(1/n)      a_n/a_(n-1)")
prev = None
for n, an in enumerate(A, start=1):
    root = mp.mpf(an) ** (mp.mpf(1) / n)
    rat = "   -" if prev is None else mp.nstr(mp.mpf(an) / prev, 8)
    P(f"  {n}   {an:<10} {mp.nstr(root,10):>12}   {rat:>10}")
    prev = an
P("  => ratio test: lambda_* >= 25.17 (n<=5, 仍在升)                     [E]")
P("  => root  test 只给 lambda_* >= 12.29 (280175^(1/5))                  [E]")
P("     (root<ratio: a_n^(1/n) 仍在大幅上升; 25.17 的正确用法是相邻比的下界,")
P("      需要比值单调性才升级为定理 —— 语料无. 故 lambda_*>=25.17 是 [E] 不是 [T].)")
P(f"  lambda_A (a5/a4)      = {mp.nstr(LAM_A,20)}")
P(f"  lambda_B (B6/B5)      = {mp.nstr(LAM_B,20)}")
P(f"  lambda_env = e*Delta  = {mp.nstr(LAM_ENV,20)}")

# ------------------------------------------------------- 2  失效点 beta_c
P("\n[2] 失效点: lambda_* u = 1  <=>  u = 1/lambda_*  <=>  beta_c = 8/lambda_*")
P("    (题面写 1/(8 lambda_*); 按 beta=8u 的正确闭式是 8/lambda_* = 2N^2/lambda_*.)")
P("    lambda_*        u_c=1/lam       beta_c=8/lam    /(1/48)     2.5/beta_c   log-frac(above)")
rows = [("A  a5/a4", LAM_A), ("B  B6/B5", LAM_B), ("env e*Delta", LAM_ENV)]
for nm, lam in rows:
    uc = 1 / lam
    bc = 8 * uc                       # = 8/lam
    r48 = bc / B48
    r25 = B25 / bc
    frac = mp.log(r25) / mp.log(SPAN)
    P(f"  {nm:<12} {mp.nstr(uc,10):>12} {mp.nstr(bc,12):>14} {mp.nstr(r48,10):>12}"
      f" {mp.nstr(r25,10):>12} {mp.nstr(frac*100,6):>10}%")

P("\n  恒等式（可核验, [T] 算术）:")
P("    beta_c/(1/48) = 1/(lambda_* u)|_(1/48) ;   2.5/beta_c = (lambda_* u)|_(2.5)")
for nm, lam in rows:
    bc = 8 / lam
    P(f"   {nm:<12} lam*u(1/48)={mp.nstr(lam*B48/8,10):>12}  -> recip {mp.nstr(1/(lam*B48/8),10):>10}"
      f"   lam*u(2.5)={mp.nstr(lam*B25/8,10):>10}   1/48<bc<2.5 ? "
      f"{'yes' if B48 < bc < B25 else 'NO'}")

# ------------------------------------------------- 3  120x 的三因子分解
P("\n[3] 120x 的三因子分解（本路唯一的产出）")
P("    1/48  --f1-->  1/lambda_*  --f2-->  8/lambda_*  --f3-->  2.5")
P("    f1 = 48/lambda_*  (已证窗边 1/48  vs 半径边 1/lam; KP 型充分边)")
P("    f2 = 8 = K_hat = 2N^2  (Penrose/KP 判据常数 vs 真半径 —— 纯松弛)")
P("    f3 = 2.5*lambda_*/8 = lambda_* u(2.5)  (半径之外: 熵必然压过能量)")
for nm, lam in rows:
    f1 = 48 / lam
    f2 = KHAT
    f3 = B25 * lam / 8
    P(f"  {nm:<12} f1={mp.nstr(f1,10):>10}  f2={f2}  f3={mp.nstr(f3,10):>10}"
      f"   prod={mp.nstr(f1*f2*f3,14):>14}   lnf3/ln120={mp.nstr(mp.log(f3)/mp.log(SPAN)*100,6)}%")

# --------------------------------------- 4  两个端点的熵/能量积（题面 9.74）
P("\n[4] 端点处的熵/能量积 lambda_* u")
for nm, lam in rows:
    P(f"  {nm:<12} (1/48): {mp.nstr(lam*B48/8,8):>12}   2.5: {mp.nstr(lam*B25/8,8):>12}"
      f"   [u=beta/4 口径: {mp.nstr(lam*B25/4,8):>10}]")
P("  题面 beta=2.5 处 '9.74' 落在 [7.87, 9.93] 内 => 对应 lambda_* ~ 31.2  [E];")
P("  语料 beta=2.5 的 lambda_s u = 15.75 是 u=beta/4=0.625 口径 (RS3:111).")

# ------------------------------------------------- 5  (b) 范围 + 稳健性
P("\n[5] (b) lambda_* 的不确定性 -> beta_c 范围")
lo = 8 / LAM_ENV                      # lambda 最大 -> beta_c 最小
hi = 8 / LAM_A
P(f"  lambda_* in [{mp.nstr(LAM_A,8)}, {mp.nstr(LAM_ENV,8)}]  =>  beta_c in"
  f" [{mp.nstr(lo,8)}, {mp.nstr(hi,8)}]")
P(f"  相对移动 beta_c(env)/beta_c(A) = {mp.nstr(lo/hi,8)}x ;  lambda_B 给"
  f" {mp.nstr(8/LAM_B,8)}")
P(f"  夹逼稳健性: 1/48 < beta_c < 2.5  <=>  3.2 < lambda_* < 384 ; 语料区间"
  f" [{mp.nstr(LAM_A,6)}, {mp.nstr(LAM_ENV,6)}] 完全落在 (3.2, 384) 内  [T]+[C]")
P(f"  熵压过能量的对数占比范围: "
  f"{mp.nstr(mp.log(B25*LAM_ENV/8)/mp.log(SPAN)*100,6)}% ~ "
  f"{mp.nstr(mp.log(B25*LAM_A/8)/mp.log(SPAN)*100,6)}%")

# ------------------------------------------------- 6  (c) 严格判据的检查
P("\n[6] (c) 严格判据")
P("  (c-1) Cauchy-Hadamard [T]: 正系数级数 sum A(alpha) u^alpha 的收敛半径")
P("        u_c = 1/limsup|A(alpha)|^(1/|alpha|) = 1/lambda_*.")
P("        => u>u_c 时 不存在任何绝对收敛的重排/重求和 (正系数). 这就是'原理性'.")
P("        出处: 语料自己的写法 RS3:119.")
P("  (c-2) 反向限定 [NG]: 半径外'必然发散'只对绝对收敛成立; 有符号抵消(条件收敛)")
P("        不被排除. 物理 CPS/Weingarten 展开的系数带符号")
P("        (PROOF-CONFINEMENT-RANDOM-SURFACE:5 'Weingarten function - signed'),")
P("        故对物理量只能宣称'诊断', 不能宣称定理.")
P("  (c-3) 判据松弛 [T]: Penrose/Kotecky-Preiss 型 K*u<1 是充分条件,")
P("        取 K=2N^2 => K u = beta => 边 beta<1/lambda_*, 比真半径边 8/lambda_*")
P("        低 8x = K_hat. 这 8x 是判据的松弛, 不是物理.")
P("  (c-4) 包络层对照 [T]+[E]: lambda_T=e*Delta => beta<1/eDelta = "
  + mp.nstr(1/LAM_ENV, 10))
P(f"        = {mp.nstr(1/LAM_ENV/B48,8)}x(1/48)  => 1/48 已在包络窗外 1.133x"
  f" (= eDelta/48 = {mp.nstr(LAM_ENV/48,8)})")
P(f"        2.5 在窗外 {mp.nstr(B25*LAM_ENV,8)}x (= ubar(beta=2.5), L1SIGN:54 的 16.989)")

# ------------------------------------------------- 7  判：多少是熵造成的
P("\n[7] 判：120x 里多少是'熵必然压过能量'")
P("  乘法分解 (lambda_A): 120 = 1.907 x 8.000 x 7.866")
P("   1.907x  = 1/48 -> 1/lambda_*   : 已证窗边与半径边之差(技术, 可用更好计数吃掉)")
P("   8.000x  = 1/lambda_* -> 8/lambda_*  : K_hat=2N^2 判据松弛(纯技术)")
P("   7.866x  = 8/lambda_* -> 2.5    : 半径之外(原理, 正系数下不可救)")
P(f"  对数测度: f3 占 ln120 的 {mp.nstr(mp.log(B25*LAM_A/8)/mp.log(SPAN)*100,6)}% ;"
  f" 技术部分占 {mp.nstr(100-mp.log(B25*LAM_A/8)/mp.log(SPAN)*100,6)}%")
P(f"  (lambda_B: 120 = {mp.nstr(48/LAM_B,8)} x 8 x {mp.nstr(B25*LAM_B/8,8)},"
  f" f3 占 {mp.nstr(mp.log(B25*LAM_B/8)/mp.log(SPAN)*100,6)}%)")
P(f"  (env:      120 = {mp.nstr(48/LAM_ENV,8)} x 8 x {mp.nstr(B25*LAM_ENV/8,8)},"
  f" f3 占 {mp.nstr(mp.log(B25*LAM_ENV/8)/mp.log(SPAN)*100,6)}%)")
P("\n  => 结论: 120x 中 >= 7.87x (对数 >= 43%) 的失效是原理性的(熵),")
P("             其余 <= 15.3x (<= 57%) 是技术性的(包络 1.9x + 判据常数 8x).")
P("\n[9] 坐标 B: 干净 Wilson 坐标 (GAP2:36  beta_W = N^2 beta_A = 4 beta_A => 1/48 读作 1/12)")
P("    同一 beta_c=8/lam, 分母换成 30x 的跨度; 绝对段 f3 不变, 只变 f1 与占比")
B12 = mp.mpf(1) / 12
P(f"    1/12 = {mp.nstr(B12,10)}   1/12 与 1/48 差 N^2 = {mp.nstr(B12/B48,6)}")
for nm, lam in rows:
    bc = 8 / lam
    f1 = (1 / lam) / B12               # (1/lam)/(1/12) = 12/lam
    f3 = B25 * lam / 8
    P(f"  {nm:<12} bc={mp.nstr(bc,8):>10}  f1={mp.nstr(f1,8):>9}  f2=8  f3={mp.nstr(f3,8):>9}"
      f"  prod={mp.nstr(f1*8*f3,12):>12}  lnf3/ln30={mp.nstr(mp.log(f3)/mp.log(30)*100,6)}%")
P("  题面 120x 的分母 vs 30x 的分母: 占比 43.1% -> 60.6% (lam_A); 48% -> 67.5% (lam_B);")
P("  59% -> 83% (env). 绝对量 f3 不变 => '越过 beta_c 的那 7.87x' 是坐标稳定量.")

P("\n[8] 反酸测 / 限定")
P("  - beta->0 (u->0): lambda_* u -> 0 < 1, 判据说'能量胜', 与强耦合收敛一致. 不空.")
P("  - 自由场: u 仍 ~ beta/8 但无可数聚合物 => lambda_* = 0 => 条件恒满足(退化),")
P("    这不是有效判据(酸的): 判据只对有非平凡动物族的模型有内容; 在 4D YM 上有内容.")
P("  - 限定: N=2, a=1, d=4, u=beta/8; lambda_* 真值 [E], lambda_*<=384 未证;")
P("    正系数假设只在'计数'层成立, 物理(签名)展开不继承 (c-2).")
P("=" * 78)
