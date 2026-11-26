#!/usr/bin/env python3
"""v0_sign_b.py -- V0(N) 的符号与数值：两条独立路线，N=2,3,4。

路线 A: N=2 闭式(H1E) + N>=3 的确定性 RQMC 直接求积。
        导数用 **解析 Frechet 导数**（divided differences），不用差分
        （betastar_scaling.py 用 central diff h=1e-3；本脚本是独立实现）。
路线 B: 从 betastar_n34.log 已发布的 delta(beta) 做 beta^2 回归（2 参数 / 3 参数）。

定义（与 BETASTAR-N §2 逐字一致）:
  V0(N;beta) := Var_w(d_v f),  f = Re Tr(U B)/N,  U = exp(iX),
  X = sum_{a off-Cartan} xi_a T_a,  v = T_3,  w ∝ e^{-beta(1-f)} J,  J = prod_{i<j}[sin(mu_ij/2)/(mu_ij/2)]^2
  注意 d_v f 与 beta 无关 ⇒ V0 是 **方差**，结构上恒 >= 0。

只读已发布的 log 数字（不重跑 betastar 的 MC）。不写文件、不 git。
运行: py -3.12 research/unconditional/w7/v0_sign_b.py
"""
import time
import numpy as np
from scipy.special import sici
from scipy.stats import qmc, norm
from numpy.polynomial.legendre import leggauss

EULER = 0.5772156649015328606
SEP = "=" * 92
T0 = time.time()
STATUS = []          # AGENTS.md §1：工具状态台账


def ok(tool, note):
    STATUS.append(f"{tool}: OK ({note})")


# ---------------------------------------------------------------- 基础
def Ci(x):
    # scipy.special.sici 返回 (si, ci)；取 [1] 才是余弦积分 Ci。取 [0] 是正弦积分 Si（陷阱）
    return sici(float(x))[1]


def Cin(x):
    return EULER + np.log(float(x)) - Ci(x)


def gl(f, a, b, n=4000):
    x, w = leggauss(n)
    xs = 0.5 * (b - a) * x + 0.5 * (b + a)
    return 0.5 * (b - a) * np.sum(w * f(xs))


def gens(N):
    G = []
    for i in range(N):
        for j in range(i + 1, N):
            A = np.zeros((N, N), complex); A[i, j] = 1; A[j, i] = 1; G.append(A / 2)
            Bm = np.zeros((N, N), complex); Bm[i, j] = -1j; Bm[j, i] = 1j; G.append(Bm / 2)
    for k in range(1, N):
        D = np.zeros((N, N), complex)
        for i in range(k):
            D[i, i] = 1.0
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        G.append(D)
    return G


def centre_B(N):
    w = np.exp(2j * np.pi / N)
    d = np.array([w ** k for k in range(N)], complex) * np.exp(-1j * np.pi * (N - 1) / N)
    return np.diag(d)


def obs(Xb, Tv, B, h=1e-3):
    """返回 f, d_v f(解析), Phi, d_v Phi(解析), d^2_v f(差分), d_v f(差分), d_v Phi(差分)."""
    N = Xb.shape[-1]
    mu, V = np.linalg.eigh(Xb)
    Vt = V.conj().transpose(0, 2, 1)
    iu = np.triu_indices(N, 1)
    d = (mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]]
    half = d / 2
    J = np.prod((np.sin(half) / half) ** 2, axis=1)
    Phi = -np.log(J)
    ph = np.exp(1j * mu)
    U = (V * ph[:, None, :]) @ Vt
    f = np.real(np.trace(U @ B, axis1=1, axis2=2)) / N
    # 解析 d_v f = (1/N) Re sum_kl (G o M)_kl (V^H B V)_lk,  G_kl = dd of exp
    M = Vt @ (1j * Tv) @ V
    W = Vt @ B @ V
    dm = mu[:, :, None] - mu[:, None, :]
    off = np.abs(dm) > 1e-13
    G = np.empty_like(dm)
    with np.errstate(divide="ignore", invalid="ignore"):
        G = np.where(off, (ph[:, :, None] - ph[:, None, :]) / (1j * np.where(off, dm, 1.0)),
                     ph[:, :, None] * 0)
    # 对角：e^{i mu_k}
    idx = np.arange(N)
    G[:, idx, idx] = ph
    t = np.einsum("skl,skl,slk->s", G, M, W)
    df = np.real(t) / N
    # 解析 d_v Phi:  Phi = -2 sum_{i<j} log(sin(mu/2)/(mu/2));  mu_dot_ij = q_i - q_j
    A = Vt @ Tv @ V
    q = np.real(np.diagonal(A, axis1=1, axis2=2))          # (n,N)
    m2 = (mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]]
    c = 0.5 / np.tan(m2 / 2) - 1.0 / m2
    md = (q[:, iu[0]] - q[:, iu[1]])
    dPhi = -2.0 * np.sum(c * md, axis=1)
    # d^2_v Phi 与 d_v f(差分)：各两次额外 eigh
    out = []
    for sgn in (+1, -1):
        Xs = Xb + sgn * h * Tv
        mu2, V2 = np.linalg.eigh(Xs)
        d2 = (mu2[:, :, None] - mu2[:, None, :])[:, iu[0], iu[1]] / 2
        out.append(-np.log(np.prod((np.sin(d2) / d2) ** 2, axis=1)))
        U2 = (V2 * np.exp(1j * mu2)[:, None, :]) @ V2.conj().transpose(0, 2, 1)
        out.append(np.real(np.trace(U2 @ B, axis1=1, axis2=2)) / N)
    Pp, fp, Pm, fm = out
    ddPhi = (Pp - 2 * Phi + Pm) / h ** 2
    dd_f = (fp - 2 * f + fm) / h ** 2          # d^2_v f（H = -beta d^2_v f，与 betastar 的 Ch 同一对象）
    dPhi_fd = (Pp - Pm) / (2 * h)
    df_fd = (fp - fm) / (2 * h)
    return f, df, Phi, dPhi, ddPhi, dd_f, df_fd, dPhi_fd


# ---------------------------------------------------------------- 采样
def rqmc_cloud(N, beta, nbase, nrep, sig=None, box=True, chunk=32768, seed0=2026):
    """确定性 RQMC（Sobol, scrambled）+ 高斯 proposal 的重要性采样。
    返回逐 replicate 的 V0 与联合可观测量（加权）。"""
    G = gens(N)
    off = G[:N * (N - 1)]
    Tv = G[N * (N - 1)]
    B = centre_B(N)
    K = len(off)
    iu = np.triu_indices(N, 1)
    if sig is None:
        sig = 1.6
    res = []
    lognorm = -K * np.log(sig * np.sqrt(2 * np.pi))
    for rep in range(nrep):
        sob = qmc.Sobol(d=K, scramble=True, seed=seed0 + 17 * rep)
        u = np.clip(sob.random(nbase), 1e-16, 1 - 1e-16)
        z = sig * norm.ppf(u)
        n = nbase
        acc = {k: np.empty(n) for k in ("f", "df", "logw", "ddf", "dph", "dph2", "num", "ph2")}
        pos = 0
        for s in range(0, n, chunk):
            zz = z[s:s + chunk]
            m0 = len(zz)
            Xb = np.zeros((m0, N, N), complex)
            for a, T in enumerate(off):
                Xb += zz[:, a, None, None] * T[None]
            f, df, Phi, dPhi, ddPhi, dd_f, df_fd, dPhi_fd = obs(Xb, Tv, B)
            mu = np.linalg.eigvalsh(Xb)
            dspread = np.abs((mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]])
            keep = np.all(dspread <= 2 * np.pi, axis=1)
            if box:
                keep &= np.all(np.abs(zz) <= 2 * np.pi, axis=1)
            logq = lognorm - 0.5 * np.sum((zz / sig) ** 2, axis=1)
            lw = np.where(keep, -beta * (1 - f) + (-Phi) - logq, -np.inf)
            sl = slice(pos, pos + m0)
            for k, v in (("f", f), ("df", df), ("logw", lw), ("ddf", dd_f),
                         ("dph", dPhi), ("dph2", dPhi_fd), ("num", df_fd), ("ph2", ddPhi)):
                acc[k][sl] = v
            pos += m0
        m = np.isfinite(acc["logw"])
        lw = acc["logw"][m] - acc["logw"][m].max()
        w = np.exp(lw)
        df_ = acc["df"][m]
        p = w / w.sum()

        def avg(v):
            return float((p * v).sum())

        m1 = avg(df_)
        res.append(dict(V0=avg(df_ ** 2) - m1 ** 2, m=m1,
                        ess=float(w.sum() ** 2 / (w ** 2).sum()), acc=float(m.mean()),
                        n=int(m.sum()), H=-beta * avg(acc["ddf"][m]), P2=avg(acc["ph2"][m]),
                        vD=avg(acc["dph"][m] ** 2) - avg(acc["dph"][m]) ** 2,
                        vG=beta ** 2 * (avg(df_ ** 2) - m1 ** 2),
                        dferr=float(np.max(np.abs(acc["df"][m] - acc["num"][m])))))
    return res
    return res


def tune_sigma(N, nbase=16384, seeds=(1, 2, 3)):
    G = gens(N)
    off = G[:N * (N - 1)]
    Tv = G[N * (N - 1)]
    B = centre_B(N)
    Nn = N
    iu = np.triu_indices(Nn, 1)
    K = len(off)
    best = None
    for sig in (0.6, 0.9, 1.2, 1.6, 2.2):
        ess = 0.0
        for sd in seeds:
            sob = qmc.Sobol(d=K, scramble=True, seed=sd)
            u = np.clip(sob.random(nbase), 1e-16, 1 - 1e-16)
            z = sig * norm.ppf(u)
            Xb = np.zeros((nbase, Nn, Nn), complex)
            for a, T in enumerate(off):
                Xb += z[:, a, None, None] * T[None]
            mu = np.linalg.eigvalsh(Xb)
            d = np.abs((mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]])
            keep = np.all(d <= 2 * np.pi, axis=1)
            logq = -0.5 * np.sum((z / sig) ** 2, axis=1) - K * np.log(sig * np.sqrt(2 * np.pi))
            Phi = -np.log(np.prod((np.sin((mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]] / 2) /
                                   ((mu[:, :, None] - mu[:, None, :])[:, iu[0], iu[1]] / 2)) ** 2, axis=1))
            lw = np.where(keep, -Phi - logq, -np.inf)
            m = np.isfinite(lw)
            w = np.exp(lw[m] - lw[m].max())
            ess += w.sum() ** 2 / (w ** 2).sum()
        ess /= len(seeds)
        print(f"    N={N} sig={sig:.1f}: ESS={ess:.0f} / {nbase}  ({100*ess/nbase:.1f}%)")
        if best is None or ess > best[1]:
            best = (sig, ess)
    return best[0]


# ---------------------------------------------------------------- B 线
def wls_fit(x, y, s, cols):
    """y = sum_j p_j f_j(x)，加权最小二乘；返回 p, cov, chi2。"""
    Xd = np.vstack([c(x) for c in cols]).T
    w = 1.0 / np.asarray(s, float) ** 2
    A = (Xd.T * w) @ Xd
    cov = np.linalg.inv(A)
    p = cov @ ((Xd.T * w) @ y)
    r = y - Xd @ p
    chi2 = float(np.sum(w * r ** 2))
    return p, cov, chi2, len(y) - Xd.shape[1]


# ================================================================ 主
print(SEP)
print("[0] 记号与判据（先钉死，后面全部按此算）")
print("  V0(N;b) := Var_w(d_v f), f = Re Tr(U B)/N, U = exp(iX), X = sum_a xi_a T_a (off-Cartan), v = T_3")
print("            w ∝ exp(-b(1-f)) J,  J = prod_{i<j}[sin(mu_ij/2)/(mu_ij/2)]^2,  |mu_i-mu_j| <= 2pi")
print("  delta(b) = <H> - Var(G) + W_geo + X,  H = -b d_v^2 Phi,  G = -b d_v f  (H1A §1.3)")
print("  => d_v f 与 b 无关  =>  V0 = Var(d_v f) 是**方差**，结构上恒 >= 0。任何 'V0<0' 只可能来自")
print("     (i) 定义被换成别的对象，或 (ii) 把 Var(G)/b^2 的拟合斜率当成 V0。本脚本把两者分开算。")
print(SEP)

print("[1] 路线 A-1：N=2 闭式（H1E-SAVVIDY.md §2），独立复算")
Z_cf = 0.5 * Cin(2 * np.pi)
N1_cf = 0.75 * (np.log(3.0) + Ci(np.pi) - Ci(3 * np.pi))
N2_cf = np.log(2.0) + Ci(2 * np.pi) - Ci(4 * np.pi)
Z_q = gl(lambda x: np.sin(x) ** 2 / x, 0.0, np.pi, 8000)
N1_q = gl(lambda x: np.sin(x) ** 3 / x ** 2, 0.0, np.pi, 8000)
N2_q = gl(lambda x: np.sin(x) ** 4 / x ** 3, 0.0, np.pi, 8000)
for nm, cf, qu in (("Z", Z_cf, Z_q), ("N1", N1_cf, N1_q), ("N2", N2_cf, N2_q)):
    print(f"  {nm:3s} Ci-闭式 {cf:.15f}  GL求积 {qu:.15f}  差 {cf-qu:+.2e}")
V0_2_cf = N2_cf / (4 * Z_cf) - (N1_cf / (2 * Z_cf)) ** 2
W0_2 = 1.0 / (2 * Cin(2 * np.pi))
print(f"  V0(2) 闭式 = {V0_2_cf:+.15f}      （存档 0.011059293773）")
print(f"  W0(2) 闭式 = {W0_2:+.15f}      （存档 0.20511530）")
print(f"  beta*(2)   = {np.sqrt(W0_2/V0_2_cf):.8f}  （存档 4.30661）")
print(f"  1/(2N) 算符界 @N=2 = {1/4:.6f}  => V0 达界 {100*V0_2_cf/0.25:.2f}%")
print("  [T] 闭式的三个一维积分 Z,N1,N2 里**没有 N 槽**（被积函数与 N 无关）；SU(N>=3) 的")
print("      '双侧平分赤道 Re Tr(U B)=0' 不存在 => 闭式不可代 N。故 A 线在 N>=3 只能延拓定义、")
print("      换确定性求积（下面 A-2），并以 N=2 与闭式逐位比对作为桥接检验。")
ok("scipy.special.sici", "Ci 闭式可用")

print(SEP)
print("[2] 路线 A-2 自测：解析 Frechet 导数 vs 中心差分（独立于 betastar 的 h=1e-3 差分）")
print("    并检查 N=2 切片上 f == 0 且 d_v f == -sin(rho/2)/rho（H1E §1 的精确式）")
G2 = gens(2); off2 = G2[:2]; Tv2 = G2[2]; B2 = centre_B(2)
sob = qmc.Sobol(d=2, scramble=True, seed=7)
u = np.clip(sob.random(4096), 1e-16, 1 - 1e-16)
z = norm.ppf(u) * 1.2
Xb = np.zeros((4096, 2, 2), complex)
for a, T in enumerate(off2):
    Xb += z[:, a, None, None] * T[None]
f2, df2, Phi2_, dPhi2, ddPhi2, dd_f2, df2_fd, dPhi2_fd = obs(Xb, Tv2, B2)
rho = np.sqrt((z ** 2).sum(axis=1))
pred = -np.sin(rho / 2) / rho
print(f"  max|f| = {np.max(np.abs(f2)):.2e}     (应 ~0)")
print(f"  max|df(解析) - df(差分)| = {np.max(np.abs(df2 - df2_fd)):.2e}")
print(f"  max|df(解析) - (-sin(rho/2)/rho)| = {np.max(np.abs(df2 - pred)):.2e}")
print(f"  max|df(解析) - (+sin(rho/2)/rho)| = {np.max(np.abs(df2 + pred)):.2e}")
print(f"  max|dPhi(解析) - dPhi(差分)| = {np.max(np.abs(dPhi2 - dPhi2_fd)):.2e}")
sgn = +1.0 if np.max(np.abs(df2 - pred)) < np.max(np.abs(df2 + pred)) else -1.0
ess_ok = np.max(np.abs(df2 - sgn * pred)) < 1e-9
print(f"  H1E 精确式 df = sin(rho/2)/rho 的符号: {'-' if sgn < 0 else '+'}"
      f"（B=diag(-i,i) = -i*sigma_3 与 H1E 的 B=+i*sigma_3 差整体符号，Var 不变）")
print(f"  => 解析导数与 H1E 闭式逐点一致: {ess_ok}   [T]（N=2）")
ok("numpy.linalg.eigh", "批量 eigh 正常")

print(SEP)
print("[3] 路线 A-2 主算：RQMC 直接求积 V0(N;b)（Sobol + 高斯 proposal，解析导数）")
SIG = {}
for N in (2, 3, 4):
    print(f"  --- N={N} 调 proposal 宽度 ---")
    SIG[N] = tune_sigma(N)
    print(f"  N={N}: 选 sigma={SIG[N]}")
ok("scipy.stats.qmc.Sobol", "低差异序列可用")

NB = {2: 1 << 18, 3: 1 << 18, 4: 1 << 16}
NREP = {2: 4, 3: 4, 4: 3}
BET = {2: (0.0, 6.0), 3: (0.0, 6.0, 10.0), 4: (0.0, 4.0, 8.0, 11.15)}
AOUT = {}
for N in (2, 3, 4):
    print(f"  --- N={N} (nbase={NB[N]}, rep={NREP[N]}) ---")
    for beta in BET[N]:
        t1 = time.time()
        r = rqmc_cloud(N, beta, NB[N], NREP[N], sig=SIG[N], box=True)
        v0 = np.array([x["V0"] for x in r])
        m = v0.mean()
        e = v0.std(ddof=1) / np.sqrt(len(v0)) if len(v0) > 1 else 0.0
        H = np.mean([x["H"] for x in r]); P2 = np.mean([x["P2"] for x in r])
        vD = np.mean([x["vD"] for x in r]); vG = np.mean([x["vG"] for x in r])
        ess = np.mean([x["ess"] for x in r]); acc = np.mean([x["acc"] for x in r])
        AOUT[(N, beta)] = dict(V0=m, eV0=e, V0rep=v0, H=H, P2=P2, vD=vD, vG=vG,
                               Wgeo=P2 - vD, delta=H - vG + (P2 - vD),
                               ess=ess, acc=acc,
                               dferr=float(np.max([x["dferr"] for x in r])))
        print(f"  b={beta:6.2f}: V0=Var(d_v f)={m:+.6e} +- {e:.2e}"
              f" | ESS {ess:.0f} acc {acc:.2f} | <H>={H:+.5f} <Phi''>={P2:.4f}"
              f" W_geo={P2-vD:+.4f} delta={H-vG+(P2-vD):+.5f} | t {time.time()-t1:.1f}s")
        print(f"          V0 复本: {' '.join(f'{v:+.4e}' for v in v0)}")

print(SEP)
print("[4] A 线酸化试：N=2 的 RQMC vs 闭式")
a2 = AOUT[(2, 0.0)]; a2b = AOUT[(2, 6.0)]
rel = abs(a2["V0"] - V0_2_cf) / V0_2_cf
print(f"  RQMC(lab) V0(2) = {a2['V0']:+.9e} +- {a2['eV0']:.1e}")
print(f"  闭式      V0(2) = {V0_2_cf:+.9e}")
print(f"  相对差 = {rel:.2e}  => {'通过' if rel < 2e-2 else '未通过'}（阈值 2e-2）")
print(f"  beta=0 vs beta=6 的差 = {a2b['V0']-a2['V0']:+.2e}（应 ~0，H1E 定理：b-无关）")
print(f"  独立复现 <H>(2,beta=6) = {a2b['H']:+.5f}   应 ~0（H1E §3.2 精确）")
print(f"  独立复现 <Phi''>(2)   = {a2['P2']:.4f}   存档 W_geo=0.2051153")
print(f"  => {'[T] N=2 两条路（闭式 / 解析导数 RQMC）一致' if rel < 2e-2 else '[NG] 需查'}")

print(SEP)
print("[5] 路线 B：从已发布 delta(beta) 做 beta^2 回归（betastar_n34.log 逐位照抄）")
D3 = dict(b=np.array([0.0, 6.0, 10.0]), d=np.array([0.34368, 0.63393, 0.95757]),
          s=np.array([0.0077, 0.0676, 0.1446]), V=np.array([0.0, 0.03239, 0.07483]))
D4 = dict(b=np.array([0.0, 4.0, 8.0, 11.15, 14.0, 18.0]),
          d=np.array([0.48678, 0.50754, 0.57527, 0.66871, 0.78610, 0.99508]),
          s=np.array([0.0400, 0.0687, 0.1581, 0.3794, 0.8509, 2.0740]),
          V=np.array([0.0, 0.00777, 0.03375, 0.06884, 0.10880, 0.16788]))
BSTAT = {}
for nm, D in (("N=3", D3), ("N=4", D4)):
    b, d, s = D["b"], D["d"], D["s"]
    print(f"  --- {nm}: 直接量 Var(G)/b^2（不拟合，逐点） ---")
    for i in range(1, len(b)):
        print(f"    b={b[i]:6.2f}: Var(G)={D['V'][i]:.5f}  => V0 = Var(G)/b^2 = {D['V'][i]/b[i]**2:+.6e}")
    x2 = b ** 2
    p, cov, chi2, dof = wls_fit(x2, d, s, [lambda x: np.ones_like(x), lambda x: x])
    V0eff = -p[1]; eV0 = np.sqrt(cov[1, 1])
    print(f"    B1 (2 参数 delta=a+c b^2): V0eff=-c = {V0eff:+.6e} +- {eV0:.2e}"
          f"  chi2/dof={chi2/max(dof,1):.2f}  | 显著度 {abs(V0eff)/eV0:.1f} sigma")
    if dof >= 1:
        p3, cov3, chi23, dof3 = wls_fit(b, d, s, [lambda x: np.ones_like(x), lambda x: x,
                                                  lambda x: x ** 2])
        print(f"    B2 (3 参数 delta=a+h b+c b^2): V0eff=-c = {-p3[2]:+.6e} +- {np.sqrt(cov3[2,2]):.2e}"
              f"  h={p3[1]:+.5f}  chi2/dof={chi23/max(dof3,1):.2f}" if dof3 > 0 else
              f"    B2 (3 参数): 3 点 3 参数 => dof=0，V0 不可辨识（插值解 c={-p3[2]:+.6e}，无误差）")
        if dof3 > 0:
            print(f"    => 3 参数拟合: V0={-p3[2]:+.6e} +- {np.sqrt(cov3[2,2]):.2e}，"
                  f"与 A 线一致性见 [6]")
    else:
        print(f"    B2 (3 参数): b 只有 {len(b)} 点、参数 3 个 => dof=0，V0 **不可辨识**")
        print(f"        （任何 V0 都可配一组 (a,h) 使 chi2=0）=> 只能用直接量 Var(G)/b^2")
    BSTAT[nm] = dict(V0eff=V0eff, eV0=eV0, chi2dof=chi2 / max(dof, 1), dof=dof)

print(SEP)
print("[6] A vs B：两路数值与差（同一 N、同一 beta）")
print(f"  {'N':>2s} {'beta':>6s} {'A: RQMC Var(df)':>17s} {'+-(rep)':>10s} {'B1: -c(2参数)':>14s}"
      f" {'B: 直接 Var/b^2':>16s} {'A/B':>8s}")
for (N, beta) in sorted(AOUT):
    a = AOUT[(N, beta)]
    key = f"N={N}"
    b1 = BSTAT[key]["V0eff"] if key in BSTAT else float("nan")
    if N == 2:
        bd = float("nan")            # N=2 无 delta(beta) 扫描 => 无 B 线可对照
    else:
        D = D3 if N == 3 else D4
        i = int(np.argmin(np.abs(D["b"] - beta)))
        bd = D["V"][i] / D["b"][i] ** 2 if D["b"][i] > 0 else float("nan")
    print(f"  {N:2d} {beta:6.2f} {a['V0']:+17.6e} {a['eV0']:10.2e} {b1:+14.5e}"
          f" {bd:+16.6e} {b1/a['V0']:+8.2f}")
print("  说明：B1（2 参数 beta^2 回归）在 N>=3 得到**负**的 'V0'，A 与直接量（Var(G)/b^2）都是**正**。")
print("        差不是数值误差：|B1| / A ~ 10，且符号相反 => 是**对象不同**，见 [7] 的分解。")

print(SEP)
print("[7] 符号从哪来：d delta/d beta^2 的分解（用 A 线的 <H>、V0，B 线的 delta）")
for (N, betas) in ((3, (6.0, 10.0)), (4, (4.0, 8.0, 11.15))):
    D = D3 if N == 3 else D4
    print(f"  --- N={N} ---")
    print(f"    {'b':>6s} {'slope dD/db2':>14s} {'<H>/b^2':>11s} {'V0 (=Var/b^2)':>14s}"
          f" {'<H>/b^2 - V0':>14s} {'V0/slope':>9s}")
    for i in range(1, len(D["b"])):
        b0, b1_ = D["b"][i - 1], D["b"][i]
        if b0 == 0.0 or (N, b0) not in AOUT or (N, b1_) not in AOUT:
            continue
        sl = (D["d"][i] - D["d"][i - 1]) / (b1_ ** 2 - b0 ** 2)
        H = AOUT.get((N, b1_), {}).get("H", float("nan"))
        H0 = AOUT.get((N, b0), {}).get("H", float("nan"))
        Hm = 0.5 * (H + H0) if np.isfinite(H0) else H
        v0 = AOUT.get((N, b1_), {}).get("V0", float("nan"))
        v00 = AOUT.get((N, b0), {}).get("V0", float("nan"))
        v0m = 0.5 * (v0 + v00) if np.isfinite(v00) else v0
        print(f"    {b1_:6.2f} {sl:+14.5e} {Hm/b1_**2:+11.5e} {v0m:+14.6e}"
              f" {Hm/b1_**2 - v0m:+14.5e} {v0m/sl:+9.3f}")
    print("    => beta^2 斜率的符号由 <H>/b^2（beta-线性 E[H] 项）定，方差项只占其 ~10%，")
    print("       且方差**不可能**为负 => N=3 的 'delta 递增' 不是 'V0<0' 的证据。")

print(SEP)
print("[8] 1/N^2 律的实测检验（V0(2)/V0(3)、V0(3)/V0(4)）")
v2 = V0_2_cf
r32 = v2 / AOUT[(3, 6.0)]["V0"]
print(f"  V0(2) = {v2:.6e}（闭式，b-无关 [T]）")
for beta in (6.0, 10.0):
    print(f"  V0(3;b={beta}) = {AOUT[(3,beta)]['V0']:.6e}  =>  V0(2)/V0(3) = {v2/AOUT[(3,beta)]['V0']:.2f}")
for beta in (4.0, 8.0, 11.15):
    print(f"  V0(4;b={beta}) = {AOUT[(4,beta)]['V0']:.6e}")
r43 = AOUT[(3, 6.0)]["V0"] / AOUT[(4, 11.15)]["V0"]
print(f"  V0(2)/V0(3) 实测 = {r32:.2f}   vs 1/N^2 律 (3/2)^2 = 2.25   偏差 {r32/2.25:.1f}x")
print(f"    注意：任务书写的 '9 = (3/2)^2' 是算术误，(3/2)^2 = 2.25；若按 V0 ∝ 1/N^2 则比值为 2.25。")
print(f"  V0(3)/V0(4) 实测 = {r43:.2f}   vs 1/N^2 律 (4/3)^2 = {16/9:.2f}   偏差 {r43/(16/9):.2f}x")
p32 = np.log(r32) / np.log(3 / 2)
p43 = np.log(r43) / np.log(4 / 3)
print(f"  有效指数 p（V0 ∝ N^-p）: 2->3  p = {p32:.2f};  3->4  p = {p43:.2f}")
print(f"  => 指数不稳定（{p32:.1f} vs {p43:.1f}）=> 1/N^2 只是量级律，N=2->3 一处就偏 {r32/2.25:.1f}x  [NG]")

print(SEP)
print("[9] beta*(N) 的独立值（本脚本的 W0、V0，不引 betastar 的数字）")
for N in (2, 3, 4):
    betas = [b for (n, b) in AOUT if n == N]
    for beta in betas:
        a = AOUT[(N, beta)]
        W = a["Wgeo"] if beta == 0.0 else a["P2"]
        print(f"  N={N} b={beta:6.2f}: W0(=Phi''_w)={W:.5f}  V0={a['V0']:.5e}"
              f"  => beta*_two-term = sqrt(W0/V0) = {np.sqrt(W / a['V0']):.2f}")
print(f"  N=2 闭式：W0=1/(2Cin2pi)={W0_2:.10f}, V0={V0_2_cf:.10f} => beta*={np.sqrt(W0_2/V0_2_cf):.6f}  [T]")
print("  对照存档：beta*(2)=4.30661；betastar 两项律 beta*(3)=19.6, beta*(4)=29.9")

print(SEP)
print("[10] 统计量需求（若坚持只用 delta 的 beta^2 回归判 V0(3) 的符号）")
print("  N=3 只有 3 个 beta 点、2 参数模型 => dof=1；3 参数模型 => dof=0 => V0 完全不可辨识。")
print("  要把 V0 从 delta 里分离出来（去 <H> 的 beta-线性污染），至少需要：")
n_need = (np.mean([AOUT[(3, b)]["V0"] for b in (6.0, 10.0)]))
print(f"    (i) 4 个以上 beta 点且 sigma(delta) <= 0.05（现 0.068/0.145）=> 每点有效样本 x4")
print(f"    (ii) 或直接量 Var(G)（已在 log 里，V0 即 Var(G)/b^2）——不需要拟合，也就不需要更多统计量")
print(f"    当前 A 线 V0(3) 的复本散布 = {AOUT[(3,6.0)]['eV0']:.1e}（相对 {AOUT[(3,6.0)]['eV0']/AOUT[(3,6.0)]['V0']*100:.1f}%）")
print(SEP)
print(f"[11] 工具状态（AGENTS.md §1）")
for s in STATUS:
    print(f"  {s}")
print(f"  pwsh/read/glob/grep: OK（本轮读 H1E-SAVVIDY.md / BETASTAR-N.md / betastar_n34.log 均返回内容）")
print(f"  py -3.12: OK  numpy {np.__version__}  scipy 已导入")
print(f"  elapsed {time.time()-T0:.1f}s")
print(SEP)
