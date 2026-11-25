# -*- coding: utf-8 -*-
"""
BC-BETHE / T0 : ferromagnetic Ising, Bethe lattice (z=4) vs exact Z^2 strip.

Bethe (z-regular tree, z=4, matched to the square lattice degree):
    cavity recursion  x -> ((t+x)/(1+t x))^{z-1},  t = tanh(beta J)
    disordered fixed point x=0 stable for (z-1) t < 1  =>  beta_c = atanh(1/(z-1))
    mass (inverse correlation length along a branch):
        m_tree = -log( (z-1) tanh(beta) )      [(z-1)t < 1]
               = 0                            [ordered]

Z^2 (exact, Onsager): row-to-row transfer matrix on a periodic strip of width W
    T = V^{1/2} (2 sinh 2K)^{W/2} exp(K* sum sigma^x_i) V^{1/2},
    V = exp(K sum sigma^z_i sigma^z_{i+1}),  K* = -1/2 log tanh K.
    lambda_0/lambda_1 gives the strip correlation length; W -> infinity is the
    bulk 2D mass.  Computed with a matrix-free eigensolver (no 2^W x 2^W matrix).

Purpose: fix the SIGN of the tree-vs-lattice mass comparison in the FKG /
positive-association class, to be contrasted with the Gaussian result
(bc_gaussian.py), where the sign is opposite.  A comparison whose sign depends
on the model class is not a geometric theorem.
"""
import json
import math

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


def m_tree(beta, z=4):
    t = math.tanh(beta)
    x = (z - 1) * t
    return -math.log(x) if x < 1.0 else 0.0


def beta_c_tree(z=4):
    return math.atanh(1.0 / (z - 1))


def ising_strip_mass(W, K, tol=1e-12):
    """log(lambda_0/lambda_1) for the periodic W x infinity Ising strip, exact."""
    if K <= 0:
        return float("inf")
    N = 1 << W
    idx = np.arange(N, dtype=np.int64)
    bits = (idx[:, None] >> np.arange(W, dtype=np.int64)) & 1
    s = 1.0 - 2.0 * bits                       # sigma^z eigenvalues +-1
    E = K * np.sum(s * np.roll(s, -1, axis=1), axis=1)
    Vdiag = np.exp(E)
    sqrtV = np.sqrt(Vdiag)
    Ks = -0.5 * math.log(math.tanh(K))
    c, sh = math.cosh(Ks), math.sinh(Ks)

    def apply_hx(v):
        """exp(Ks sum_i sigma^x_i) v  =  tensor product of 2x2 matrices."""
        a = v.reshape((2,) * W)
        for ax in range(W):
            a = np.moveaxis(a, ax, 0)
            new0 = c * a[0] + sh * a[1]
            new1 = sh * a[0] + c * a[1]
            a = np.stack([new0, new1], axis=0)
            a = np.moveaxis(a, 0, ax)
        return a.reshape(-1)

    def matvec(v):
        return sqrtV * apply_hx(sqrtV * v)

    op = LinearOperator((N, N), matvec=matvec, dtype=np.float64)
    k = 2 if N > 8 else min(N, 2)
    vals = eigsh(op, k=k, which="LA", tol=tol, return_eigenvectors=False,
                 maxiter=20000, ncv=min(N, 24))
    vals = np.sort(vals)[::-1]
    return float(math.log(vals[0] / vals[1]))


def main():
    out = {}
    z = 4
    bc = beta_c_tree(z)
    print("=" * 78)
    print("BC-BETHE T0  --  ferromagnetic Ising: Bethe(z=4) vs exact Z^2 strip")
    print("=" * 78)
    print(f"\nbeta_c^tree(z=4) = atanh(1/3) = {bc:.6f}   ;  beta_c(Z^2) = 0.440687")
    out["beta_c_tree"] = bc
    out["beta_c_square"] = 0.44068679350977147

    # strip convergence check at one beta
    print("\n[1] strip convergence (m_lat(W) should converge as W -> infinity)")
    print(f"{'beta':>6} " + " ".join(f"{'W='+str(w):>11}" for w in (8, 10, 12)))
    conv = {}
    for beta in (0.15, 0.25, 0.35):
        K = beta
        vals = [ising_strip_mass(w, K) for w in (8, 10, 12)]
        conv[beta] = vals
        print(f"{beta:>6.2f} " + " ".join(f"{v:>11.6f}" for v in vals))
    out["strip_convergence"] = {str(k): v for k, v in conv.items()}

    # main table, W = 12
    print("\n[2] masses, beta below both critical points")
    print(f"{'beta':>6} {'m_lat(Z^2,W=12)':>16} {'m_tree(z=4)':>12} {'tree/lat':>9} {'tree-lat':>9}")
    rows = []
    for beta in (0.10, 0.15, 0.20, 0.25, 0.30, 0.32, 0.34):
        ml = ising_strip_mass(12, beta)
        mt = m_tree(beta, z)
        rows.append({"beta": beta, "m_lat": ml, "m_tree": mt,
                     "ratio": mt / ml, "diff": mt - ml})
        print(f"{beta:>6.2f} {ml:>16.6f} {mt:>12.6f} {mt/ml:>9.4f} {mt-ml:>9.4f}")
    out["mass_table"] = rows

    # exact Onsager cross-check: axial mass m = 2|K - K*|, K* = -1/2 log tanh K
    print("\n[3] exact Onsager cross-check of the strip:  m_lat = 2|K - K*|,  K* = -1/2 log tanh K")
    chk = []
    for r in rows:
        b = r["beta"]
        Ks = -0.5 * math.log(math.tanh(b))
        pred = 2.0 * abs(Ks - b)
        chk.append({"beta": b, "m_lat_strip": r["m_lat"], "m_lat_onsager": pred,
                    "diff": r["m_lat"] - pred})
        print(f"    beta={b:>5.2f}  strip={r['m_lat']:.6f}  Onsager={pred:.6f}  diff={r['m_lat']-pred:+.2e}")
    out["onsager_check"] = chk
    print("    (agreement to ~1e-3: the W=12 strip is already the bulk mass at these beta)")

    print("\n[4] reading")
    print("    * in the whole disordered window  m_tree < m_lat  (tree is the WORST case:")
    print("      more paths, slower decay, smaller mass); the tree gap closes FIRST.")
    print("    * beta_c^tree = 0.3466 < beta_c^lat = 0.4407 : the tree orders EARLIER,")
    print("      i.e. at fixed beta it is CLOSER to criticality and therefore LESS massive.")
    print("    * this is the OPPOSITE sign to the Gaussian model (bc_gaussian.py: m_tree > m_lat).")
    print("    * consequence: any 'tree -> lattice' inequality is a statement about the FKG /")
    print("      positive-association class, not about connectivity; SU(N) Wilson has no FKG.")

    with open(__file__.replace(".py", ".json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
