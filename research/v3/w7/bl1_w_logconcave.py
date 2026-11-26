# BL1 — 核 SOS1 §4 的「IBP NS=1/4 + 1D Brascamp-Lieb ⇒ Var<=0.018832, beta<3.3003」
# 纪律：先核 w 的 log-凹性（决定 BL 能否用），再复算 Var/margin/window。
# 只读不写其它文件。mpmath。py -3.12.
from mpmath import mp, mpf, sin, cos, quad, pi, sqrt, log, findroot, diff

mp.dps = 45

# ---------- q11 的模型（Q11-CERT §1.1-§1.3）----------
# u = rho/2 in (0,pi);  theta=pi 切片上 S_W == beta, 故权重是常数
# 测度 = 纯 Haar 切片对 phi 积分后的 u-边缘:  w(u) du
def w(u):    return sin(u)**2/u
def g(u):    return sin(u)/(2*u)                  # G = beta*s*g(u),  theta=pi: s=1
def gp(u):   return (u*cos(u) - sin(u))/(2*u**2)  # g'
# (BL-1) 参考密度 w 的 log-凹模量:  phi := -log w
def phi(u):  return log(u) - 2*log(sin(u))
def phi2(u): return 2/sin(u)**2 - 1/u**2          # phi'' = -(log w)''
# q11 的 Phi''(= xi3 方向的 Hessian 分量, 不是 u 径向!)  —— 用来演示"别拿错曲率"
def Phi2_xi3(u): return 1/(2*u**2) - cos(u)/(sin(u)*2*u)

PTS = [mpf(0), mpf(1), mpf(2), mpf(3), pi]
def Q(f): return quad(f, PTS)

print("== [A] q11 锚点自核（用纯 mp.quad；对照 Q11-CERT §2.2 的区间值）==")
Z  = Q(w)
NA = Q(lambda u: g(u)*w(u))
NB = Q(lambda u: g(u)**2*w(u))
Eg, Eg2 = NA/Z, NB/Z
V0 = Eg2 - Eg**2
print(f"Z      = {mp.nstr(Z,30)}   (q11: 1.218826696528612205905355171)")
print(f"E[g]   = {mp.nstr(Eg,30)}   (q11: 0.3574113533914995816)")
print(f"E[g^2] = {mp.nstr(Eg2,30)}   (q11: 0.1388021693061169894)")
print(f"Var(g) = {mp.nstr(V0,30)}   (q11: 0.0110592937729735899)")
S = 1/(4*Z)
NS_num = Q(lambda u: Phi2_xi3(u)*w(u))
print(f"NS(=∫Phi''w) = {mp.nstr(NS_num,30)}   (解析 1/4)")
print(f"S = NS/Z = 1/(4Z) = {mp.nstr(S,30)}   (q11: 0.205115297123073148)")

print()
print("== [B] w 的 log-凹性：phi'' = 2/sin^2 u - 1/u^2 的符号  ==")
# 网格扫描（含端点附近）
N = 20000
m = (None, None)
for k in range(1, N):
    u = pi*k/N
    v = phi2(u)
    if m[0] is None or v < m[0]:
        m = (v, u)
print(f"网格 min phi'' = {mp.nstr(m[0],20)}  at u = {mp.nstr(m[1],12)}   (N={N})")
# 精确驻点: d/du phi'' = -4cos u/sin^3 u + 2/u^3 = 0
F  = lambda u: -4*cos(u)/sin(u)**3 + 2/u**3
ustar = findroot(F, mpf('1.39'))
rmin  = phi2(ustar)
print(f"驻点 u* = {mp.nstr(ustar,18)}   phi''(u*) = {mp.nstr(rmin,25)}")
print(f"phi''>0 全程?  {phi2(mpf('0.5'))>0 and phi2(mpf('2'))>0 and phi2(mpf('3'))>0 and rmin>0}")
print(f"对照：q11 的 Phi''(xi3 分量) 的 min = {mp.nstr(min(Phi2_xi3(pi*k/N) for k in range(1,N)),12)}  <- 拿这个当 rho 就把 BL 算废")
print(f"对照：Hess S_W 的 min eig(θ=0,β=2.5) = -β/4 = {mp.nstr(-mpf('2.5')/4,12)}  <- f03/sos1 的 -0.625，与 w 无关")

print()
print("== [C] BL: Var(f) <= E[f'^2]/rho_min  (log-凹 ⇒ 1D Brascamp-Lieb) ==")
Egp2 = Q(lambda u: gp(u)**2*w(u))/Z
print(f"E[g'^2] = {mp.nstr(Egp2,30)}   (sos1: 0.02917444834)")
BLc   = Egp2/rmin
print(f"rho_min = {mp.nstr(rmin,20)}")
print(f"Var <= E[g'^2]/rho_min = {mp.nstr(BLc,30)}   (sos1: 0.01883152401)")
print(f"真 Var(g) = {mp.nstr(V0,20)}  <= 上界? {V0 < BLc}   ratio={mp.nstr(BLc/V0,8)}")
# sharpened: 用函数曲率 1/phi''
Esh = Q(lambda u: gp(u)**2/phi2(u)*w(u))/Z
print(f"sharpened E[g'^2/phi''] = {mp.nstr(Esh,30)}   (sos1: 0.01487231)")

print()
print("== [D] 边距与窗口（beta=2.5 与临界 beta）==")
b = mpf('2.5')
for name, vb in (("rho_min 版", BLc), ("sharpened 版", Esh)):
    LB = S - b**2*vb
    bt = sqrt(S/vb)
    print(f"{name}: beta^2*Var <= {mp.nstr(b**2*vb,20)};  G'' >= S - that = {mp.nstr(LB,20)};  窗口 beta < {mp.nstr(bt,12)}")
print(f"sos1 报: +0.0874183 / 3.3003218 ; q11 精确 beta* = 4.30660785256")

print()
print("== [E] 酸测：同一 rho_min 对多个 f 是否都成立（BL 不能只对 g 成立）==")
tests = [("u", lambda u: u, lambda u: mpf(1)),
         ("u^2", lambda u: u**2, lambda u: 2*u),
         ("sin u", sin, cos),
         ("cos u", cos, lambda u: -sin(u)),
         ("g", g, gp),
         ("sin u/u", lambda u: sin(u)/u, lambda u: (u*cos(u)-sin(u))/u**2)]
for nm, f, fp in tests:
    Ef  = Q(lambda u: f(u)*w(u))/Z
    Ef2 = Q(lambda u: f(u)**2*w(u))/Z
    Vf  = Ef2 - Ef**2
    Ef2p= Q(lambda u: fp(u)**2*w(u))/Z
    print(f"f={nm:9s} Var={mp.nstr(Vf,12):>16s}  E[f'^2]/rho_min={mp.nstr(Ef2p/rmin,12):>16s}  ok={Vf <= Ef2p/rmin}")

print()
print("== [F] 窗口纪律：只有 theta=pi 时权重才与 beta 无关 ==")
print("  theta=pi: c=cos(pi/2)=0 ⇒ S_W≡beta ⇒ 测度 = 纯 Haar 切片 w/Z（与 beta 无关）⇒ rho_min=phi'' 最小 = 上式")
print("  theta=0 (c=1): 测度 ∝ exp(-beta(1-cos u)) w du ⇒ 曲率 = beta*cos u + phi''，min 会随 beta 降：")
for bb in ('2.5', '3.0', '3.5', '3.8'):
    bq = mpf(bb)
    mm = min(bq*cos(pi*k/20000) + phi2(pi*k/20000) for k in range(1, 20000))
    print(f"    beta={bb}: min(beta cos u + phi'') = {mp.nstr(mm,10)}")
# beta 临界: min_u [beta cos u + phi''] = 0
h = lambda bb: min(mpf(bb)*cos(pi*k/20000) + phi2(pi*k/20000) for k in range(1, 20000))
lo, hi = mpf('3.0'), mpf('6.0')
for _ in range(40):
    mid = (lo+hi)/2
    if h(mid) > 0: lo = mid
    else: hi = mid
print(f"   beta*(theta=0, log-凹失效) ≈ {mp.nstr((lo+hi)/2,10)}   （只作纪律演示，不是本条的窗口）")

print()
print("== [G] 判 ==")
print("  w = sin^2 u / u  log-凹（phi''>=1.5492346>0）⇒ 1D BL 可用 ⇒ Var<=0.0188315, margin +0.08742, 窗口 beta<3.30032")
print("  -0.625 = min eig Hess S_W (xi-空间的 3x3)，不是 u-线曲率，不进入 BL 假设")
