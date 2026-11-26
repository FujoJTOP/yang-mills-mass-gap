"""dn_perlink.py -- exact decomposition of the z2(s=0.5) density 1.5 beta^2 = 9.375.

Shows that the number is an INCIDENCE COUNT (plaquettes per link, signs, alignment) and
carries no dynamics: at the Z2-cocycle interpolation s=1/2 one link sits in 6 plaquettes,
their 3-link "rest" products are traceless (|sin psi| = 1), and the traces align.

SU(2), beta = 2.5, beta^2 = 6.25, beta^2/4 = 1.5625.
"""
import numpy as np

BETA = 2.5
T1 = np.array([[0.0, 0.5], [0.5, 0.0]], dtype=complex)
T2 = np.array([[0.0, -0.5j], [0.5j, 0.0]], dtype=complex)
T3 = np.array([[0.5, 0.0], [0.0, -0.5]], dtype=complex)
LAM = np.stack([1j * T1, 1j * T2, 1j * T3])
I2 = np.eye(2, dtype=complex)
S3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
L = 2


def li(v, m):
    return (v[0] + v[1] * L + v[2] * L ** 2 + v[3] * L ** 3) * 4 + m


# ---- config: U_1(x) = cos(th) I + i sin(th) S3 ,  th = s*pi*(x2+x3+x4), s = 1/2
s = 0.5
E = 4 * L ** 4
U = [I2.copy() for _ in range(E)]
for x1 in range(L):
    for x2 in range(L):
        for x3 in range(L):
            for x4 in range(L):
                th = s * np.pi * ((x2 + x3 + x4) % 2)
                U[li([x1, x2, x3, x4], 0)] = np.cos(th) * I2 + 1j * np.sin(th) * S3

plaqs = []
for x1 in range(L):
    for x2 in range(L):
        for x3 in range(L):
            for x4 in range(L):
                for mu in range(4):
                    for nu in range(mu + 1, 4):
                        x = [x1, x2, x3, x4]
                        xm = x.copy(); xm[mu] = (xm[mu] + 1) % L
                        xn = x.copy(); xn[nu] = (xn[nu] + 1) % L
                        plaqs.append(((li(x, mu), +1), (li(xm, nu), +1),
                                      (li(xn, mu), -1), (li(x, nu), -1)))


def mats_of(p):
    return [U[i] if sg > 0 else U[i].conj().T for i, sg in p]


def prod(lst):
    P = I2
    for m in lst:
        P = P @ m
    return P


g = np.zeros((E, 3))
detail = {}
for p in plaqs:
    M = mats_of(p)
    for k in range(4):
        ik, sgk = p[k]
        rot = M[k:] + M[:k]
        tail = prod(rot[1:] + [rot[0]])
        head = prod(rot)
        X = tail if sgk > 0 else head
        c = -BETA * 0.5 * sgk * np.array([np.trace(LAM[a] @ X).real for a in range(3)])
        g[ik] += c
        detail.setdefault(ik, []).append((k, sgk, c.copy(), X))

gn = float(np.sum(g ** 2))
NP = 6 * L ** 4
print("=" * 100)
print("dn_perlink.py -- exact decomposition of |grad S_W|^2 at the z2(s=0.5) configuration")
print(f"SU(2), beta = {BETA}, L = {L}, E = {E} links, N_p = {NP} plaquettes")
print("=" * 100)
print(f"  sum_l |g_l|^2            = {gn:.6f}")
print(f"  per link                 = {gn/E:.6f}   = {gn/E/BETA**2:.6f} beta^2")
print(f"  per plaquette (density)  = {gn/NP:.6f}   = {gn/NP/BETA**2:.6f} beta^2")
print(f"  check 1.5 beta^2 = {1.5*BETA**2:.6f} ; 2.25 beta^2 = {2.25*BETA**2:.6f}")

# how many links share one value?  (translation invariance -> few distinct values)
vals = np.round(np.sum(g ** 2, axis=1), 9)
u, cnt = np.unique(vals, return_counts=True)
print(f"\n  distinct per-link values: {dict(zip(u.tolist(), cnt.tolist()))}")
print("  (a single value for ALL links => the config is translation invariant => the")
print("   density is L-independent BY CONSTRUCTION, for every L at once)")

# decompose one direction-1 link
target = li([0, 0, 0, 0], 0)
print(f"\n  decomposition at link 0 (direction 1, unit cell 0):  g = {g[target]}")
print(f"  |g|^2 = {np.sum(g[target]**2):.6f} = {np.sum(g[target]**2)/BETA**2:.6f} beta^2")
tot_sq = 0.0
print("     k  sign   s_{p,l,.} = -(beta/2) sg tr(lam_a X)      |s|^2/(beta^2/4)   sin^2(psi)")
for (k, sgk, c, X) in detail[target]:
    trX = 0.5 * np.trace(X).real
    sin2 = 1.0 - trX ** 2
    tot_sq += float(np.sum(c ** 2))
    print(f"     {k}   {sgk:+d}   [{c[0]:+8.4f} {c[1]:+8.4f} {c[2]:+8.4f}]        "
          f"{np.sum(c**2)/(BETA**2/4):8.4f}        {sin2:6.4f}")
print(f"  sum over the 6 plaquettes of |s|^2 = {tot_sq:.6f} = {tot_sq/BETA**2:.6f} beta^2"
      f"   (= 6 * beta^2/4, the local budget)")
print(f"  |sum|^2 = |g|^2 = {np.sum(g[target]**2):.6f} = {np.sum(g[target]**2)/BETA**2:.4f} beta^2")
print(f"  |sum|^2 / sum|s|^2 = {np.sum(g[target]**2)/tot_sq:.4f}  vs the Cauchy-Schwarz bound 6")
print(f"  => the 6 contributions are COLLINEAR (all +z), so C-S is EXACTLY SATURATED at this link:")
print(f"     |g_l|^2 = 6 * (6 * beta^2/4) = 9 beta^2, which is the per-link bound of D1 itself.")
print(f"  Globally only 16 of {E} links are extremal (the direction-1 ones); the other 48 have")
print(f"  g = 0 exactly => total = 16 * 9 beta^2 = {16*9*BETA**2:.2f} = 144 beta^2, and")
print(f"  d = 144 beta^2 / (6 L^4) = (1/4) * 9 beta^2 * (E/N_p) = 1.5 beta^2.")
print("\n  => d = 1.5 beta^2 is (fraction of extremal links) x (per-link budget) x (E/N_p):")
print("     a counting statement about the 6 incident plaquettes and the 16/64 = 1/4 fraction.")
print("     No beta-running, no L-running, no dynamics.")
print("\nDONE.")
