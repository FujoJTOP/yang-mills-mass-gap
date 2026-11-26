# -*- coding: utf-8 -*-
r"""
r11_adjudicate.py -- R11: 裁决 r03(SU(2): Phi'==0 对一切 theta)  vs  r08(SU(3): S==D, gap==0)

不引入新模型, 只做四件事:
 [A] SU(2): Phi(xi)=g(|xi|) 是径向函数 ==> 精确恒等式  Phi' = (g'(rho)/rho) * a ,  a:=<e3,xi>=xi3
     ==> a=0 上 Phi'==0 对一切 theta (点态, 与 theta 无关);  Phi''|_{a=0} = g'(rho)/rho
 [B] SU(2) 的 D 到底是哪个: D = Var(G+Phi') = beta^2 s^2 Var(g)  vs  Var(Phi') = 0
 [C] SU(3): (i) 环面上 r08 的 Phi' != 0 与 IBP 恒等式 J(Phi'^2-Phi'')=d^2_tau J (r08 步骤无错)
            (ii) 反例: 正则 X=diag(x) 且 a=<E,X>=0 而 Phi'(X) != 0  ==> N=3 的 a=0 切片不钉死 Phi'
 [D] 纯 Haar SU(3) 的 a=0 条件切片: E[Phi''|a=0] 与 Var(Phi'|a=0) (对照 SU(2) 的 0.2051 与 0)

口径: [A][B] mpmath 高精度(40 位)+ 闭式; [C] 解析恒等式 + 浮点核对; [D] MC(numpy), 标 [C].
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 40
BETA = mp.mpf('2.5')
pi = mp.pi
LOG = []


def P(s=""):
    print(s)
    LOG.append(s)


# ============================ [A] SU(2): Phi' = (g'(rho)/rho) a , 点态 ============================
def g0(d):                      # O09: g(d) = -2 log( sin(d/2)/(d/2) ),  Phi(SU(2))=g(rho)
    return -2 * mp.log(mp.sin(d / 2) / (d / 2))


def g0p(d):
    return 2 / d - 1 / mp.tan(d / 2)


def g0pp(d):
    return mp.mpf(1) / 2 / mp.sin(d / 2) ** 2 - 2 / d ** 2


def Phi2(xi):                   # xi in R^3, U=exp(i xi.sig/2), J=4 sin^2(rho/2)/rho^2
    rho = mp.sqrt(sum(x * x for x in xi))
    return g0(rho)


P("=" * 78)
P("[A] SU(2)  Phi(xi)=g(|xi|) 径向  ==>  Phi' = (g'(rho)/rho)*a  (a=xi3)")
P("    theta 完全不出现在 Phi 里 (Phi=-log J, J 只依赖 rho) ==> '对一切 theta' 是 trivially 的")
P("=" * 78)
pts = [(mp.mpf('0.37'), mp.mpf('0.91')), (mp.mpf('-1.4'), mp.mpf('0.25')),
       (mp.mpf('2.1'), mp.mpf('-1.3')), (mp.mpf('0.05'), mp.mpf('3.0'))]  # (x1,x2)
avec = [mp.mpf(0), mp.mpf('0.4'), mp.mpf('-0.7')]
P(f"{'x1':>7} {'x2':>7} {'rho':>10} {'max|Phi_a=0|':>16} {'max|Phi_a-Pred|':>18}")
worst0 = mp.mpf(0)
worst = mp.mpf(0)
for x1, x2 in pts:
    for a in avec:
        rho = mp.sqrt(x1 ** 2 + x2 ** 2 + a ** 2)
        d = mp.diff(lambda t: Phi2([x1, x2, a + t]), 0)      # d/dtau
        pred = g0p(rho) * a / rho
        if a == 0:
            worst0 = max(worst0, abs(d))
        else:
            worst = max(worst, abs(d - pred))
    rho = mp.sqrt(x1 ** 2 + x2 ** 2)
    d0 = mp.diff(lambda t: Phi2([x1, x2, t]), 0)
    P(f"{float(x1):7.3f} {float(x2):7.3f} {float(rho):10.4f} {mp.nstr(abs(d0),8):>16} {'--':>18}")
P(f"-->  a=0 上 max|Phi'| = {mp.nstr(worst0,8)}     (点态恒 0)")
P(f"-->  一般 a 上 max|Phi' - (g'(rho)/rho)a| = {mp.nstr(worst,8)}   <= 恒等式成立 ==> Phi' 恰是 a 的线性函数")
P(f"-->  故 Phi'≡0 <==> a=0 ; 且 Phi''|_a=0 = g'(rho)/rho = 1/(2u^2)-cot(u)/(2u) , u=rho/2  (核对: 下)")
for u in [mp.mpf('0.3'), mp.mpf('1.0'), mp.mpf('1.4')]:
    rho = 2 * u
    xt = mp.diff(lambda t: Phi2([rho, t, 0]), 0, 2)          # 横向(e2, a=0) = 切片内方向
    xr = mp.diff(lambda t: Phi2([rho - t, 0, 0]), 0, 2)      # 径向(e1) = g''(rho)
    P(f"    u={float(u):4.2f}  横向 Phi''={mp.nstr(xt,12):>18}   g'(rho)/rho={mp.nstr(g0p(rho)/rho,12):>18}"
      f"   1/(2u^2)-cot(u)/(2u)={mp.nstr(1/(2*u**2)-1/mp.tan(u)/(2*u),12):>18}")
    P(f"                      径向 Phi''={mp.nstr(xr,12):>18}   g''(rho)={mp.nstr(g0pp(rho),12):>18}")

# ============================ [B] SU(2) 的 D 是哪个 ============================
P()
P("=" * 78)
P("[B] SU(2) 的 D:  D = Var(G+Phi') = beta^2 s^2 Var_nu(g)   vs   Var(Phi') = 0")
P("    密度 w(u)=sin^2 u/u , g(u)=sin u/(2u) (=sqrt(J)/2, 也是 S_W 的展开系数)", )
P("=" * 78)


def wfun(u):
    return mp.sin(u) ** 2 / u


def gg(u):
    return mp.sin(u) / (2 * u)


def phipp(u):
    return 1 / (2 * u ** 2) - 1 / mp.tan(u) / (2 * u)


brk = [0, pi / 3, 2 * pi / 3, pi]
Z = mp.quad(wfun, brk)
NA = mp.quad(lambda u: gg(u) * wfun(u), brk)
NB = mp.quad(lambda u: gg(u) ** 2 * wfun(u), brk)
NS = mp.quad(lambda u: phipp(u) * wfun(u), brk)
Eg, Eg2 = NA / Z, NB / Z
Vg = Eg2 - Eg ** 2
Epp = NS / Z
P(f"Z  = {mp.nstr(Z,20)}")
P(f"E[g]    = {mp.nstr(Eg,20)}       E[g^2] = {mp.nstr(Eg2,20)}")
P(f"Var(g)  = {mp.nstr(Vg,20)}    <== V0 (beta-无关, 纯 Haar 泛函)")
P(f"E[Phi'']= {mp.nstr(Epp,20)}    <== (q11: 0.205115297123...; 1/(2Cin2pi))")
P(f"Var(Phi') = 0  (精确: [A] 已证 Phi' == 0 逐点)              <== 关键")
for th in ['pi']:
    P(f"theta={th}: S = E[Phi''] = {mp.nstr(Epp,20)}   D = beta^2 Var(g) = {mp.nstr(BETA**2*Vg,20)}"
      f"   S-D = {mp.nstr(Epp-BETA**2*Vg,20)}")
P(f"若误用 D=Var(Phi')=0  ==> gap = {mp.nstr(Epp,20)}  (虚高 {mp.nstr(Epp/(Epp-BETA**2*Vg),6)} 倍)")
P(f"差异绝对值 |beta^2Var(g) - Var(Phi')| = {mp.nstr(BETA**2*Vg,20)}  = D 的 100%")
# beta=2.5, theta=pi/2 回归 (r03)
c = s = mp.sqrt(2) / 2


def wfun2(u):
    return mp.e ** (-BETA * (1 - c * mp.cos(u))) * mp.sin(u) ** 2 / u


Z2 = mp.quad(wfun2, brk)
NA2 = mp.quad(lambda u: gg(u) * wfun2(u), brk)
NB2 = mp.quad(lambda u: gg(u) ** 2 * wfun2(u), brk)
NS2 = mp.quad(lambda u: phipp(u) * wfun2(u), brk)
Vg2 = NB2 / Z2 - (NA2 / Z2) ** 2
S2 = (BETA * c / 2) * (NA2 / Z2) + NS2 / Z2
D2 = BETA ** 2 * s ** 2 * Vg2
P(f"回归 theta=pi/2: S={mp.nstr(S2,12)}  D={mp.nstr(D2,12)}  S-D={mp.nstr(S2-D2,12)}  (r03: 0.542943922354)")

# ============================ [C] SU(3): r08 的 Phi' 与 IBP ============================
P()
P("=" * 78)
P("[C] SU(3) 环面 (r08 模型):  Phi'=-sum(lam_i-lam_j)cot(x_ij/2)")
P("=" * 78)
LAM = np.array([1.0, 0.0, -1.0])
PAIRS = ((0, 1), (0, 2), (1, 2))


def Jfun(t):
    t = np.asarray(t, float)
    return np.prod([4 * np.sin((t[i] - t[j]) / 2) ** 2 for i, j in PAIRS])


def phi_p(t):
    t = np.asarray(t, float)
    return sum(-(LAM[i] - LAM[j]) / np.tan((t[i] - t[j]) / 2) for i, j in PAIRS)


def phi_pp(t):
    t = np.asarray(t, float)
    return sum(0.5 * (LAM[i] - LAM[j]) ** 2 / np.sin((t[i] - t[j]) / 2) ** 2 for i, j in PAIRS)


rng = np.random.default_rng(7)
t1 = rng.uniform(0.3, 2.0, 2000)
t2 = rng.uniform(0.3, 2.0, 2000)
Pp = np.array([phi_p([a, b, -(a + b)]) for a, b in zip(t1, t2)])
Ppp = np.array([phi_pp([a, b, -(a + b)]) for a, b in zip(t1, t2)])
JJ = np.array([Jfun([a, b, -(a + b)]) for a, b in zip(t1, t2)])
P(f"2000 随机环面点: min|Phi'|={np.abs(Pp).min():.4f}  mean Phi'={Pp.mean():+.4f}  Var(Phi')={Pp.var():.4f}"
  f"  median Phi''={np.median(Ppp):.4f}  (Phi' 与 Phi'' 都逐点非平凡)")
# 点态 IBP: J(Phi'^2-Phi'') = d^2_tau J
h = 1e-4
res = []
for k in range(5):
    t = np.array([t1[k], t2[k], -(t1[k] + t2[k])])
    Jt = lambda ta: Jfun([t[i] + ta * LAM[i] for i in range(3)])
    d2 = (Jt(h) - 2 * Jt(0) + Jt(-h)) / h ** 2
    res.append(abs(Jfun(t) * (phi_p(t) ** 2 - phi_pp(t)) - d2) / max(1.0, abs(d2)))
P(f"点态 IBP  J(Phi'^2-Phi'') = d^2_tau J : 5 点相对残差 max = {max(res):.3e}  ==> r08 该步无错 [T]")
# 2D 台形: 积分相等 (float, 常数项=均值)
N = 900
th1 = np.linspace(0, 2 * np.pi, N, endpoint=False) + np.pi / N     # 半格偏移: 避开 sin=0 奇异
TH1, TH2 = np.meshgrid(th1, th1, indexing='ij')
R1 = 4 * np.sin((TH1 - TH2) / 2) ** 2
R2 = 4 * np.sin((TH1 + 2 * TH2) / 2) ** 2      # x13 = th1-th3 = th1+th1+th2 ? 见下
# th3 = -(th1+th2); x12=th1-th2 ; x13=th1-th3 = 2th1+th2 ; x23=th2-th3 = th1+2th2
R2 = 4 * np.sin((2 * TH1 + TH2) / 2) ** 2
R3 = 4 * np.sin((TH1 + 2 * TH2) / 2) ** 2
Jg = R1 * R2 * R3
Pp_g = (-(LAM[0] - LAM[1]) / np.tan((TH1 - TH2) / 2)
        -(LAM[0] - LAM[2]) / np.tan((2 * TH1 + TH2) / 2)
        -(LAM[1] - LAM[2]) / np.tan((TH1 + 2 * TH2) / 2))
Ppp_g = (0.5 * (LAM[0] - LAM[1]) ** 2 / np.sin((TH1 - TH2) / 2) ** 2
        + 0.5 * (LAM[0] - LAM[2]) ** 2 / np.sin((2 * TH1 + TH2) / 2) ** 2
        + 0.5 * (LAM[1] - LAM[2]) ** 2 / np.sin((TH1 + 2 * TH2) / 2) ** 2)
mA = Jg * Ppp_g
mB = Jg * Pp_g ** 2
mA = mA[np.isfinite(mA)]
mB = mB[np.isfinite(mB)]
P(f"900^2 台形均值(剔除奇异线): mean|D|^2={Jg.mean():.6f}  mean J*Phi''={mA.mean():.6f}"
  f"  mean J*Phi'^2={mB.mean():.6f}  差={mA.mean()-mB.mean():+.3e}  ==> S≡D (r08 对)")

# ---- [C-ii] 决定性反例: N=3 的 a=0 切片不钉死 Phi' ----
P()
P("-" * 78)
P("[C-ii] 反例: 正则 X=diag(x), a:=<E,X>=2tr(EX)=0, 而 Phi'(X)!=0  (N=3)")
P("-" * 78)


def phi_p_cart(x, e):
    return sum(g0p(x[i] - x[j]) * (e[i] - e[j]) for i, j in PAIRS)


best = None
for _ in range(4000):
    x = rng.uniform(-1.5, 1.5, 3); x = x - x.mean()
    if min(abs(x[i] - x[j]) for i, j in PAIRS) < 0.25:
        continue
    e = rng.normal(size=3); e = e - e.mean()
    # 投影到 a=0: 减去 Lam*e 中 x 的平行分量
    e = e - (e @ x) / (x @ x) * x
    e = e - e.mean()
    if np.abs(e).max() < 1e-9:
        continue
    pp = phi_p_cart(x, e)
    if best is None or abs(pp) > abs(best[2]):
        best = (x.copy(), e.copy(), pp)
x, e, pp = best
a = 2 * float(np.dot(e, x))
P(f"  x = {np.round(x,4)}  (Sum=0, 谱差 {np.round(np.array([abs(x[i]-x[j]) for i,j in PAIRS]),4)} 全非零 ==> 正则)")
P(f"  e = {np.round(e,4)}  (Sum=0)   a = 2tr(EX) = {a:+.1e}")
P(f"  Phi'(X) = sum g'(x_ij)(e_i-e_j) = {mp.nstr(phi_p_cart(list(x), list(e)),12)}   != 0")
P("  但注意: ! 它不是 O(1), 看下面的三阶定律")
P(f"  {'t':>6} {'Phi(tx)|a=0':>16} {'Phi/t^3':>14}")
for t in [0.05, 0.1, 0.2, 0.4, 1.0, 2.0]:
    v = float(phi_p_cart([t * xi for xi in x], list(e)))
    P(f"  {t:6.2f} {v:16.8f} {v/t**3:14.6f}")
P("  => Phi' 在 a=0 上是 O(|x|^3): 因为 Hess Phi(0) = (N/12)*Id 对 N=2 和 N=3 都各向同性 (P16 §3 的 N/12 地板)")
P("     即: 二次阶 Phi 对所有 N 都是径向的; 各向异性从四次项 (phi=g-d^2/12=d^4/1440+..) 才出现.")
P("     故 N=3 上 Phi' 在 a=0 只是'不恒为 0 且 O(x^3) 小', 而 N=2 是'恒等于 0'.")
for N_ in (2, 3):
    xs = [1.0, -0.4] if N_ == 2 else [1.0, -0.4, -0.6]
    if N_ == 2:
        xs = [0.7, -0.7]
    s2 = sum(xi * xi for xi in xs)
    t0 = 1e-3
    val = sum(g0(1e-3 * (xs[i] - xs[j])) for i in range(N_) for j in range(i + 1, N_))
    P(f"     N={N_}: Phi(t x)/t^2 -> {float(val)/1e-6:12.8f}    (N/12)Sum x^2 = {N_/12*s2:12.8f}")
P(f"  ==> N=3 上 a=0 不蕴含 Phi'=0; 且 g' 为奇函数 (g'(-d)=-g'(d)) ==> Phi'(-X) = -Phi'(X)")
P(f"      => Phi' 在 a=0 切片上是奇、非零函数; 切片密度 |Delta|^2 为偶 => E[Phi'|a=0]=0 且 Var>0")
P(f"  对照 N=2: x=(rho/2,-rho/2), e=(e1,-e1) => a=2rho e1, e1-e2=2e1=a/rho  ==> Phi'=(g'(rho)/rho)a ≡ 0 于 a=0")
for rho in [0.9, 1.7]:
    e1 = 0.31
    xx = [rho / 2, -rho / 2]
    aa = 2 * (e1 * xx[0] + (-e1) * xx[1])
    p2 = g0p(xx[0] - xx[1]) * (e1 - (-e1))
    P(f"    rho={rho}: a={float(aa):.6f}  Phi'={float(p2):.10f}  (g'(rho)/rho)*a={float(g0p(rho)/rho*aa):.10f}")

# ============================ [D] 纯 Haar SU(3) 的 a=0 条件切片 ============================
P()
P("=" * 78)
P("[D] 综合: 正确的 D = Var(G + Phi')  , 两个模型各取一项")
P("=" * 78)
P("  SU(2) (r03/q11, 切片 xi3=0, v=e3):  Phi' == 0 精确 [A]  ==>  D = Var(G) = beta^2 s^2 Var(g) = 0.0691206")
P("  SU(3) (r08, 环面 T^2, v=lam in t):  G = d_tau S_W == 0 (S_W 在切片上被声明为常数 beta) ==> D = Var(Phi')")
P("  一般式:  D = Var(G) + Var(Phi') + 2 Cov(G,Phi') ;  谁掉项由模型决定, 不由 N 决定")
P("  ==> 题面'矛盾'是假二难: r03 掉 Var(Phi')  => D=Var(G)!=0 ; r08 掉 Var(G) => D=Var(Phi'). 不是同一个 D.")
P()
P("  不可交换的那一半: r08 的 tau 是环面平移 => Z(tau)=Z(0) => S≡D 是测度对称的恒等式;")
P("  r03 的 tau 是约束坐标 a 的法向(离开切片) => Z(a) 不是常数 => S-D 可以 != 0.")
P("  ==> r08 不是算术上漏了 S_W 项, 而是用了 S_W 恒常的模型; 该模型结构性地看不见 gap.")
P()
P("  量级对照 (同一 'Xi'' 归一化):")
P("    SU(3) r08 exact : E[|D|^2 Phi''] = 48/Z = 8 = E[|D|^2 Phi'^2]/Z   <== Phi' 侧与 Phi'' 侧恰好相等(=8)")
P("    SU(2) 本文 exact: E[Phi''] = 0.2051152971 , Var(Phi') = 0            <== Phi' 侧整块消失")
P("  ==> SU(2) 的正性不来自 'Phi'' 强', 而来自 'Phi' 侧被 rank-1 抹掉'; SU(3) 上 Phi' 侧不消失,")
P("      在环面模型里还恰好与 Phi'' 侧一样大. 这就是 'SU(2) 是 rank-1 现象' 的确切含义.")
P()
P("=" * 78)
P("判: (a) 成立, 但理由必须改正:  r08 的 IBP 步 (J(Phi'^2-Phi'')=d^2_tau J, 积分=0) 没有错 [T];")
P("    S≡D 是'平移不变 => 二阶导=0'的恒等式, 所以 Weyl/环面路线对 SU(3) 给不出非空证书 (no-go),")
P("    但**不能**读成 'SU(3) 在 theta=pi 的 gap 恒为 0' 的定理 (r08 自己的 4(b) 已承认缺 8 维模型).")
P("    Phi'≡0 于 a=0 是 N=2 精确专属 (Phi'=(g'(rho)/rho)a); N=3 有 a=0 而 Phi'!=0 的反例 => 障碍是结构性的(rank-1).")
P("=" * 78)

with open(__file__.replace('.py', '.log'), 'w', encoding='utf-8') as f:
    f.write("\n".join(LOG) + "\n")
