import numpy as np, sys

def mi(z, L):
    return z if z <= L - z else z - L

def W(Cz, L, mu, n, kap_vec):
    """kap_vec: dict {direction: momentum}.  transverse directions = all != mu."""
    nC = mi(n, L)
    tot = 0.0 + 0.0j
    others = [m for m in range(4) if m != mu]
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
                ph = (kap_vec.get(others[0], 0.0) * za
                      + kap_vec.get(others[1], 0.0) * zb
                      + kap_vec.get(others[2], 0.0) * zc)
                idx = tuple((zz % L) for zz in z)
                tot += Cz[idx] * np.exp(1j * ph)
    return tot

for L in (6, 8, 10, 12):
    d = np.load("research/unconditional/w7/bmc_L%d_b2.500.npz" % L)
    Cz = d["Cz"]
    ks = 2 * np.pi * np.arange(L) / L
    print("=" * 70)
    print("L =", L)
    print("  --- kappa = 0 的 W（应全部 >= 0） ---")
    mn = 1e30
    arg = None
    for mu in range(4):
        for n in range(L):
            v = W(Cz, L, mu, n, {})
            if v.real < mn:
                mn = v.real
                arg = (mu, mi(n, L))
    print("    min over (mu,n) of W_mu(n;0) = %+.6e at %s" % (mn, arg))
    print("  --- 全部 kappa ---")
    mn = 1e30
    arg = None
    nv = 0
    nt = 0
    for mu in range(4):
        others = [m for m in range(4) if m != mu]
        for n in range(L):
            for k1 in range(L):
                for k2 in range(L):
                    for k3 in range(L):
                        kap = {others[0]: ks[k1], others[1]: ks[k2], others[2]: ks[k3]}
                        v = W(Cz, L, mu, n, kap)
                        re = v.real
                        nt += 1
                        if re < mn:
                            mn = re
                            arg = (mu, mi(n, L), k1, k2, k3)
                        if re < -1e-12:
                            nv += 1
    print("    min over all (mu,n,kappa) = %+.6e at %s" % (mn, arg))
    print("    违反 -1e-12 的个数 = %d / %d" % (nv, nt))
    sys.stdout.flush()
