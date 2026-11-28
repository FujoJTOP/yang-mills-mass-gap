import numpy as np, sys, os

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

print("=" * 76)
print("决定性自检：Sum_n e^{i k_mu n} W_mu(n;kappa) 必须等于 C_hat(k)")
print("  若成立 => 我的 W 计算正确 => 负值是数据的真实性质")
print("  若不成立 => 我的 W 有 bug")
print("=" * 76)
for L in (6, 8):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Cz = d["Cz"]; Ck = d["Ck"]
    ks = 2 * np.pi * np.arange(L) / L
    worst = 0.0
    worst_info = None
    # 抽查若干 (mu, k) 组合
    for mu in range(4):
        others = [m for m in range(4) if m != mu]
        for k1 in range(L):
            for k2 in range(L):
                for k3 in range(L):
                    kap = {others[0]: ks[k1], others[1]: ks[k2], others[2]: ks[k3]}
                    for kmu in range(L):
                        s = 0.0 + 0.0j
                        for n in range(L):
                            s += np.exp(1j * ks[kmu] * mi(n, L)) * W(Cz, L, mu, n, kap)
                        kvec = [0, 0, 0, 0]
                        kvec[mu] = kmu
                        kvec[others[0]] = k1
                        kvec[others[1]] = k2
                        kvec[others[2]] = k3
                        ref = Ck[tuple(kvec)]
                        err = abs(s.real - ref)
                        if err > worst:
                            worst = err
                            worst_info = (mu, kmu, k1, k2, k3, s.real, ref)
        sys.stdout.write("  L=%d mu=%d done, running max err=%.3e\n" % (L, mu, worst)); sys.stdout.flush()
    print("  L=%d: max |Sum_n e^{ikn}W(n,kappa) - C_hat(k)| = %.3e" % (L, worst))
    print("         worst at %s" % (worst_info,))
    sys.stdout.flush()
