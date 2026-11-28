import numpy as np, io, sys

def minimg(L):
    i = np.arange(L); dd = np.minimum(i, L - i)
    sq = np.zeros((L, L, L, L))
    for mu in range(4):
        a = dd.reshape([L if k == mu else 1 for k in range(4)])
        sq = sq + a ** 2
    return np.sqrt(sq)

def lamgrid(L):
    lam = np.zeros((L, L, L, L))
    for mu in range(4):
        x = np.arange(L).reshape([L if k == mu else 1 for k in range(4)])
        lam = lam + 4 * np.sin(np.pi * x / L) ** 2
    return np.broadcast_to(lam, (L, L, L, L)) if lam.shape != (L, L, L, L) else lam

out = io.StringIO(); w = out.write
w("# 数据包：4D SU(2) Wilson，beta=2.5，L=6/8/10/12\n\n")
w("来源 `research/unconditional/w7/bmc_L{6,8,10,12}_b2.500.npz`。\n\n")
w("键：`Ck` = 完整 `C_hat_L(k)`（shape `L^4`）；`Cz` = 完整 `Gamma_L(0,z)`（shape `L^4`，最小镜像坐标）；\n")
w("`chi_q`、`varq`、`varbar`、`plaq`、`sweeps`、`M`、`N`、`acc`、`eps`、`seed`。\n\n")
w("**`q_fields` 为空**（`shape=(0,L,L,L,L)`）=> 无配置级数据 => 无 bootstrap、无自相关、无 `argmax` 选择偏差。\n\n")
w("定义：`q_x = (1/2) tr U_{P(x)}`；`Gamma_L(x,y)=Cov(q_x,q_y)`；`C_hat_L(k)=Sum_z Gamma_L(0,z) e^{ikz}`；\n")
w("`lambda(k)=4 Sum_mu sin^2(k_mu/2)`；`g_c(z)=(1+|z|_2)^{-4}`；`Lambda_L=1+sup_x Sum_y g_c(x-y)`。\n\n")

w("## 1. 标量汇总\n\n")
w("| L | sweeps | N | acc | seed | varq | varbar | chi_q | Cz(0) | Sum\\|Cz\\| | Lambda_L |\n")
w("|---|---|---|---|---|---|---|---|---|---|---|\n")
for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    r = minimg(L); Lam = 1.0 + float(((1.0 + r) ** -4.0).sum())
    w("| %d | %d | %d | %.4f | %d | %.8f | %.10e | %.8f | %.8f | %.6f | %.4f |\n" % (
        L, int(d["sweeps"].shape[0]), int(d["N"]), float(d["acc"]), int(d["seed"]),
        float(d["varq"]), float(d["varbar"]), float(d["chi_q"]),
        float(d["Cz"][(0,) * 4]), float(np.abs(d["Cz"]).sum()), Lam))
w("\n`Lambda_L` 含 `y=x` 项 `g_c(0)=1`。自检：`varbar * L^4 = chi_q`（FFT 归一化）。\n\n")

w("## 2. `Cz` 实空间壳层（最小镜像 `r=|z|_2`，按 `r^2` 分组）\n\n")
w("| L | r^2 | # | mean | min | max |\n|---|---|---|---|---|---|\n")
for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Cz = d["Cz"]; r = minimg(L); r2 = np.round(r ** 2, 9)
    for s in sorted(set(r2.ravel())):
        m = np.isclose(r2, s, atol=1e-9); v = Cz[m]
        w("| %d | %g | %d | %+.6e | %+.3e | %+.3e |\n" % (L, s, int(m.sum()), v.mean(), v.min(), v.max()))
w("\n")
w("## 3. `C_hat` 壳层（按 `lambda(k)` 分组）\n\n")
w("| L | lambda | # | mean | std | min | max |\n|---|---|---|---|---|---|---|\n")
for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Ck = d["Ck"]; xs = lamgrid(L).ravel(); ys = Ck.ravel()
    for s in sorted(set(np.round(xs, 9))):
        m = np.isclose(xs, s, atol=1e-9); v = ys[m]
        w("| %d | %.2f | %d | %.8f | %.2e | %.8f | %.8f |\n" % (
            L, s, int(m.sum()), v.mean(), v.std(ddof=1) if m.sum() > 1 else 0.0, v.min(), v.max()))

txt = out.getvalue()
io.open("docs/DATA-su2-wilson-b2.5.md", "w", encoding="utf-8").write(txt)
sys.stdout.write("written bytes=%d lines=%d\n" % (len(txt.encode()), txt.count("\n")))
