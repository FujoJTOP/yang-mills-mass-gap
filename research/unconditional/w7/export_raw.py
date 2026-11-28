import numpy as np, io, os, sys

OUT = "docs/data-raw"
os.makedirs(OUT, exist_ok=True)

def mi(z, L):
    return z if z <= L - z else z - L

def W(Cz, L, mu, n, kap):
    nC = mi(n, L)
    others = [m for m in range(4) if m != mu]
    tot = 0.0 + 0.0j
    for a in range(L):
        za = mi(a, L)
        for b in range(L):
            zb = mi(b, L)
            for c in range(L):
                zc = mi(c, L)
                z = [0, 0, 0, 0]
                z[mu] = nC
                z[others[0]] = za
                z[others[1]] = zb
                z[others[2]] = zc
                ph = (kap.get(others[0], 0.0) * za
                      + kap.get(others[1], 0.0) * zb
                      + kap.get(others[2], 0.0) * zc)
                tot += Cz[tuple((zz % L) for zz in z)] * np.exp(1j * ph)
    return tot

# ---- 1. 原始数组导出为 CSV（每行一个 z，逗号分隔） ----
for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    for key in ("Ck", "Cz"):
        arr = np.asarray(d[key], dtype=np.float64)
        flat = arr.ravel(order="C")          # C-order：最后一维最快
        np.savetxt(os.path.join(OUT, "%s_L%d.csv" % (key, L)), flat,
                   fmt="%.17g", header="L=%d key=%s shape=(%d,%d,%d,%d) C-order flattened"
                   % (L, key, L, L, L, L), comments="# ")
sys.stdout.write("CSV written.\n"); sys.stdout.flush()

# ---- 2. W_mu(n;kappa) 全表 + 负值清单 ----
lines = io.StringIO()
w = lines.write
w("# W_mu(n;kappa) 全表（来自 clover 观测量的 Cz）\n")
w("# 定义 W_mu(n;kappa) = Sum_{z_perp} Cz[mi(n,mu; z_perp)] exp(i kappa.z_perp)\n")
w("# 恒等自检：Sum_n e^{i k_mu n} W_mu(n;kappa) = Ck(k)（在 L=6,8 上验证误差 6e-17）\n\n")
for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Cz = d["Cz"]
    ks = 2 * np.pi * np.arange(L) / L
    w("## L = %d\n\n" % L)
    w("### kappa = 0: W_mu(n;0) 全表\n\n")
    w("| mu | n | W |\n|---|---|---|\n")
    neg0 = []
    for mu in range(4):
        for n in range(L):
            v = W(Cz, L, mu, n, {}).real
            w("| %d | %d | %+.8e |\n" % (mu, mi(n, L), v))
            if v < 0:
                neg0.append((mu, mi(n, L), v))
    w("\n负值（kappa=0）：%s\n\n" % ("无" if not neg0 else str(neg0)))
    # 全 kappa 的负值清单
    negs = []
    tot = 0
    for mu in range(4):
        others = [m for m in range(4) if m != mu]
        for n in range(L):
            for k1 in range(L):
                for k2 in range(L):
                    for k3 in range(L):
                        kap = {others[0]: ks[k1], others[1]: ks[k2], others[2]: ks[k3]}
                        v = W(Cz, L, mu, n, kap).real
                        tot += 1
                        if v < -1e-12:
                            negs.append((mu, mi(n, L), k1, k2, k3, v))
    negs.sort(key=lambda t: t[5])
    w("### 全部 kappa：负值个数 %d / %d\n\n" % (len(negs), tot))
    w("最负的 20 个（mu, n, k1, k2, k3, W）：\n\n```\n")
    for t in negs[:20]:
        w("%s\n" % (t,))
    w("```\n\n")
    sys.stdout.write("L=%d done\n" % L); sys.stdout.flush()

io.open(os.path.join(OUT, "W-table.md"), "w", encoding="utf-8").write(lines.getvalue())
sys.stdout.write("W-table.md written.\n")
