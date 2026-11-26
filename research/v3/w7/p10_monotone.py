"""P10 — Delta_dyn(beta) = E[H] - Var(G_v) 的单调性：单 plaquette SU(2) 精确一维求积。

约定逐字同 O05-FORCEVAR.md / w7/o05_var_scaling.py / H1A-EXACT.md §2.1:
  xi in R^3, |xi|<2pi; S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2)*u), u = xi_3/rho
  v = e_3, 条件化 xi_3 = 0 => 切片上 rho=r, u=0, Phi' == 0 (逐点精确) => X = 0
  G   = d_3 S   = beta*sin(th/2)*sin(r/2)/r            = beta*s*g(r),  g(r)=sin(r/2)/r
  H   = d_3^2 S = beta*cos(th/2)*sin(r/2)/(2r)         = beta*c*h(r),  h(r)=sin(r/2)/(2r)
  P'' = d_3^2 Phi = 2/r^2 - cot(r/2)/r  -> 1/6 (r->0)
  切片测度 w(r)dr = 2*pi*sin^2(r/2)/r * exp(-beta*S0(r)) dr,  S0(r) = 1 - c*cos(r/2)
  =>  Delta_dyn = beta*c*E[h] - beta^2*s^2*Var[g]                     [T]
  =>  d/dbeta E[f] = -Cov(f,S0),  d/dbeta Var[g] = -[Cov(g^2,S0) - 2E[g]*Cov(g,S0)]   [T]
只写 p10_*.py|.log（禁 git commit）。
"""
import numpy as np

BETAS = [0.5, 1.0, 2.5, 4.0, 6.0, 10.0]
THETAS = [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, 0.9 * np.pi, np.pi]
NFINE = 400001
NSCAN = 20001


def grid(n):
    r = np.linspace(1e-9, 2 * np.pi - 1e-9, n)
    Sr, Cr = np.sin(r / 2), np.cos(r / 2)
    return r, Sr, Cr, Sr / r, Sr / (2 * r), 2 * np.pi * Sr**2 / r, 2.0 / r**2 - (Cr / Sr) / r


def stats(beta, theta, n=NFINE, g=None):
    if g is None:
        g = grid(n)
    r, Sr, Cr, GG, HH, BASE, P2 = g
    return {k: v[0] for k, v in _stats(g, beta, np.array([theta])).items()}


def _stats(g, beta, thetas):
    """beta 固定, thetas 向量化; 返回 dict of arrays (nt,)"""
    r, Sr, Cr, GG, HH, BASE, P2 = g
    c = np.cos(thetas / 2.0)[:, None]
    s = np.sin(thetas / 2.0)[:, None]
    S0 = 1.0 - c * Cr
    w = BASE * np.exp(-beta * S0)
    Z = np.trapezoid(w, r, axis=1)[:, None]
    P = w / Z
    Eg = np.trapezoid(P * GG, r, axis=1)
    Eh = np.trapezoid(P * HH, r, axis=1)
    Eg2 = np.trapezoid(P * GG**2, r, axis=1)
    Vg = Eg2 - Eg**2
    W = np.trapezoid(P * P2, r, axis=1)
    ES0 = np.trapezoid(P * S0, r, axis=1)
    Cg = np.trapezoid(P * GG * S0, r, axis=1) - Eg * ES0
    Cg2 = np.trapezoid(P * GG**2 * S0, r, axis=1) - Eg2 * ES0
    Ch = np.trapezoid(P * HH * S0, r, axis=1) - Eh * ES0
    c1, s1 = c[:, 0], s[:, 0]
    dEh = -Ch
    dVg = -(Cg2 - 2 * Eg * Cg)
    D = beta * c1 * Eh - beta**2 * s1**2 * Vg
    dD = c1 * (Eh + beta * dEh) - s1**2 * (2 * beta * Vg + beta**2 * dVg)
    return dict(c=c1, s=s1, Eh=Eh, Vg=Vg, W=W, D=D, dD=dD,
                dVarG=2 * beta * Vg + beta**2 * dVg, EH=beta * s1 * Eg,
                EhH=beta * c1 * Eh, VG=beta**2 * s1**2 * Vg, dEh=dEh, dVg=dVg)


def D_fd(beta, theta, n=NSCAN, h=1e-5, g=None):
    if g is None:
        g = grid(n)
    a = _stats(g, beta * (1 + h), np.array([theta]))['D'][0]
    b = _stats(g, beta * (1 - h), np.array([theta]))['D'][0]
    return (a - b) / (2 * beta * h)


gF, gS = grid(NFINE), grid(NSCAN)

print("=" * 82)
print("[0] 对账 (theta=pi 切片 = Haar 冻结)")
st = stats(1.0, np.pi, g=gF)
V0, W0 = st['Vg'], st['W']
print(f"  Var_Haar[g] = {V0:.15f}   (存档 V0 = 0.011059293772974)")
print(f"  E[Phi'']    = {W0:.15f}   (存档 W0 = 0.205115297123)")
print(f"  sqrt(W0/V0) = {np.sqrt(W0 / V0):.8f}   (存档 beta* = 4.30660785)")
print(f"  |D + beta^2 V0| @theta=pi,beta=2.5: {abs(stats(2.5, np.pi, g=gF)['D'] + 6.25 * V0):.3e}")


print("=" * 82)
print("[1] Delta_dyn(beta; theta) = E[H] - Var(G|0) = beta*c*E[h] - beta^2*s^2*Var[g]   [T]式 / [C]值")
print("  theta/pi " + "".join(f"   b={b:<8g}" for b in BETAS) + " |  dD/db(b=6)  FD核对")
for th in THETAS:
    row = f"  {th / np.pi:6.3f}  "
    for b in BETAS:
        row += f"  {stats(b, th, g=gF)['D']:+.6f} "
    d6 = stats(6.0, th, g=gF)['dD']
    fd = D_fd(6.0, th, g=gS)
    row += f" | {d6:+.6f}  {fd:+.6f}"
    print(row)

print("=" * 82)
print("[2] 分项与 beta-导数")
print("  beta theta/pi      E[H]         Var(G)       dD/db      dVar(G)/db   dVar[g]/db")
for b in [2.5, 4.0, 6.0]:
    for th in [np.pi / 2, 0.9 * np.pi, np.pi]:
        st = stats(b, th, g=gF)
        print(f"  {b:4g}  {th / np.pi:.2f}   {st['EhH']:+.9f}  {st['VG']:+.9f}  "
              f"{st['dD']:+.9f}  {st['dVarG']:+.9f}  {st['dVg']:+.9f}")

print("=" * 82)
print("[3] 单调性判定: 网格 theta in [0,pi] x beta in [0.5,10]")
bg = np.linspace(0.5, 10.0, 61)
tg = np.linspace(0.0, np.pi, 181)
Dm = np.zeros((len(bg), len(tg)))
dDm = np.zeros((len(bg), len(tg)))
for j, b in enumerate(bg):
    o = _stats(gS, b, tg)
    Dm[j] = o['D']
    dDm[j] = o['dD']
inc = dDm > 0
print(f"  dD/dbeta>0 比例 {inc.mean():.4f} ; <0 比例 {(dDm < 0).mean():.4f}")
flip = np.where(np.sign(dDm[:, :-1]) * np.sign(dDm[:, 1:]) < 0)
if len(flip[0]):
    print(f"  dD/dbeta 沿 theta 变号: theta/pi in [{tg[flip[1]].min()/np.pi:.3f}, "
          f"{tg[flip[1]+1].max()/np.pi:.3f}], 覆盖 beta in [{bg[flip[0]].min():.2f}, {bg[flip[0]].max():.2f}]")
print("  每 beta: min_theta D / argmintheta ; (D+W)@theta=pi ; min_theta(D+W) / argmin")
for j, b in enumerate([0.5, 1.0, 2.5, 4.0, 6.0, 10.0]):
    jj = int(np.argmin(np.abs(bg - b)))
    k = int(np.argmin(Dm[jj]))
    HW = Dm[jj] + _stats(gS, b, tg)['W']
    k2 = int(np.argmin(HW))
    print(f"  {b:4g}  {Dm[jj, k]:+.6f}  {tg[k] / np.pi:.3f}      "
          f"{stats(b, np.pi, g=gF)['D'] + stats(b, np.pi, g=gF)['W']:+.6f}   "
          f"{HW[k2]:+.6f}  {tg[k2] / np.pi:.3f}")

print("=" * 82)
print("[4] Delta_dyn(beta) 在 [0.5,10] 上是否有内点极小 (每个 theta)")
print("  theta/pi   min_b D      argmin beta   D(0.5)      D(6)        D(10)      dD/db(b=10)")
bgf = np.linspace(0.5, 10.0, 800)
for th in THETAS:
    pass
    Dv = np.array([stats(b, th, g=gS)['D'] for b in bgf])
    k = int(np.argmin(Dv))
    print(f"  {th / np.pi:.3f}    {Dv.min():+.6f}   {bgf[k]:.3f}       "
          f"{Dv[0]:+.6f}  {stats(6.0, th, g=gF)['D']:+.6f}  {Dv[-1]:+.6f}  "
          f"{stats(10.0, th, g=gF)['dD']:+.6f}")

print("=" * 82)
print("[5] H1 判据 G''=D+W 的穿越 beta* (theta 扫描)")
for th in [0.0, np.pi / 4, np.pi / 2, 0.75 * np.pi, 0.9 * np.pi, np.pi]:
    b0, b1 = 0.2, 14.0
    f = lambda b: stats(b, th, g=gS)['D'] + stats(b, th, g=gS)['W']
    if f(b0) * f(b1) > 0:
        print(f"  theta={th / np.pi:.2f}pi : 端点同号 f({b0})={f(b0):+.4f}, f({b1})={f(b1):+.4f}")
        continue
    for _ in range(70):
        m = 0.5 * (b0 + b1)
        if f(b0) * f(m) <= 0:
            b1 = m
        else:
            b0 = m
    print(f"  theta={th / np.pi:.2f}pi : beta* = {0.5 * (b0 + b1):.6f}")

print("=" * 82)
print("[6] 充分条件 (任务 4): dD/db = c*(E[h]+beta E'[h]) - s^2*(2beta Var[g]+beta^2 Var'[g]) >= 0")
print("  theta/pi beta   A=c*(...)    B=s^2*(...)   A-B")
for th in [np.pi / 4, np.pi / 2, 0.75 * np.pi]:
    for b in [2.5, 4.0, 6.0]:
        st = stats(b, th, g=gF)
        A = st['c'] * (st['Eh'] + b * st['dEh'])
        B = st['s']**2 * (2 * b * st['Vg'] + b**2 * st['dVg'])
        print(f"  {th / np.pi:5.3f}  {b:4g}  {A:+.6f}   {B:+.6f}   {A - B:+.6f}")
print("=" * 82)
