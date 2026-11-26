"""
TA03-RHO2 : what is rho^2 at beta = 2.5 ?

Two blocks.

BLOCK 1 (exact, Fraction/mpmath): the four corpus anchors of rho^2
    rho^2 := mu2/mu1  (mu_n = 0++-projected coarse transfer generator levels)
           = Lambda*/m1^2  (LOOP-RUNBOOK 12.50 slot ratio)
  (1) 1b single plaquette SU(2) beta=0   : mu = l(l+2) = 3, 8        -> 8/3
  (2) 2D Yang-Mills SU(2)                 : C2(adj)/C2(fund)=2/(3/4) -> 8/3
  (3) 2D Yang-Mills SU(3)                 : 3/(4/3)                 -> 9/4
  (4) N=inf / free tower on T^4           : e1=1, e2=2              -> 2
  plus the Dirichlet string (q=0) -> 4, and the free-field-with-mass scan
  rho^2(mu)=(mu^2+2)/(mu^2+1) -> the acid test that rho^2>1 is NOT a gap criterion.

BLOCK 2 (MC, SU(2) Wilson, beta=2.5, L0=2, L=2): can rho^2 be measured?
  Only one nontrivial time separation exists (t=0,1), so a mass *ratio* needs a
  matrix: 3 0++ operators (plaq, clover, 1x2 rect), C_ij(0), C_ij(1), then
  GEVP  C(1) v = lambda C(0) v ,  lambda_n = 1/cosh(m_n)  (T=2 single-pole)
  => m_n = arccosh(1/lambda_n),  rho^2 = (m_2/m_1)^2 .
  Estimator only: same family as GEVP-4D.md, which is not one-sided.

Run: py -3.12 ta03_rho2_b25.py --n 40000 --gap 2
"""
import sys
import time
import argparse
from fractions import Fraction

import numpy as np

OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


# ---------------------------------------------------------------- BLOCK 1
def block_exact():
    P("=" * 78)
    P("BLOCK 1  exact anchors of rho^2 = mu2/mu1  (definition: ch05-theorems.tex:155,")
    P("         ch16-sixth-species.tex:527-532, app.tex:705; slot form LOOP-RUNBOOK:5545)")
    P("=" * 78)

    # (1) one plaquette, beta=0: rotor on S^3, eigenvalues l(l+2)
    mu = [l * (l + 2) for l in (1, 2, 3, 4)]
    r2 = Fraction(mu[1], mu[0])
    P(f"[T] (1) 1b single plaquette SU(2), beta=0: mu_l = l(l+2) = {mu}")
    P(f"        rho^2 = mu2/mu1 = {mu[1]}/{mu[0]} = {r2} = {float(r2):.10f}")

    # (2)(3) 2D Yang-Mills: quadratic Casimirs
    def c2_su2(j):
        return Fraction(int(4 * j * (j + 1)), 4)          # j(j+1), j = n/2

    c2_f2 = c2_su2(Fraction(1, 2))                        # 3/4
    c2_a2 = c2_su2(Fraction(1))                           # 2
    P(f"[T] (2) 2D YM SU(2): C2(fund)={c2_f2}, C2(adj)={c2_a2}, "
      f"rho^2 = {c2_a2 / c2_f2} = {float(c2_a2 / c2_f2):.10f}  (beta-independent)")
    c2_f3, c2_a3 = Fraction(4, 3), Fraction(3)            # SU(3)
    P(f"[T] (3) 2D YM SU(3): C2(fund)={c2_f3}, C2(adj)={c2_a3}, "
      f"rho^2 = {c2_a3 / c2_f3} = {float(c2_a3 / c2_f3):.10f}  (beta-independent)")

    # (4) free tower on T^4: n in Z^4, e = sum n_i^2
    n2 = []
    for a in range(-3, 4):
        for b in range(-3, 4):
            for c in range(-3, 4):
                for d in range(-3, 4):
                    n2.append(a * a + b * b + c * c + d * d)
    n2 = np.array(n2)
    P(f"[T] (4) free tower T^4: e1 = min>0 = {n2[n2 > 0].min()}, "
      f"e2 = next distinct > e1 = {np.unique(n2[n2 > 0])[1]}"
      f"  (degeneracies {int((n2 == 1).sum())}, {int((n2 == 2).sum())})")
    P(f"        rho^2 = 2/1 = {Fraction(2, 1)} = 2.0000000000   [no gap: m=0]")

    # Dirichlet string q = 0: mu_k = k^2 pi^2
    from mpmath import mp, mpf, pi
    mp.dps = 30
    r2d = (2 * pi) ** 2 / pi ** 2
    P(f"[T] extra: 1D Dirichlet string q=0: mu_k=k^2 pi^2, rho^2 = 4pi^2/pi^2 = "
      f"{mp.nstr(r2d, 12)}")

    # free field with mass mu: rho^2 = (mu^2+2)/(mu^2+1)  (B02-LARGEN.md:88)
    P("")
    P("[T] acid test (reproduce B02-LARGEN.md:88): free field with mass mu on T^4,")
    P("    e1 = mu^2+1, e2 = mu^2+2  =>  rho^2(mu) = (mu^2+2)/(mu^2+1) :")
    for m in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0):
        P(f"      mu = {m:>6} : rho^2 = {(m*m+2)/(m*m+1):.6f}"
          f"   gap m = {m}")
    P("    => rho^2 > 1 holds for EVERY mu >= 0, including mu = 0 (no gap at all).")
    P("       rho^2 > 1 is therefore not a gap criterion; at mu=0 it is 2 exactly,")
    P("       the same value as N=inf (B02-LARGEN.md:81-90).")
    P("")


# ---------------------------------------------------------------- BLOCK 2
def block_mc(args):
    sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
    import w4_blocked as W
    from su3 import su2_heatbath, SIG

    adj = W.adj

    def staple(U, mu):
        S = np.zeros_like(U[mu])
        Um = U[mu]
        for nu in range(4):
            if nu == mu:
                continue
            Un = U[nu]
            S += np.roll(Un, -1, axis=mu) @ adj(np.roll(Um, -1, axis=nu)) @ adj(Un)
            S += (adj(np.roll(np.roll(Un, 1, axis=nu), -1, axis=mu))
                  @ adj(np.roll(Um, 1, axis=nu)) @ np.roll(Un, 1, axis=nu))
        return S

    def plaq(U, mu, nu):
        return (U[mu] @ np.roll(U[nu], -1, axis=mu)
                @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))

    def update(U, beta, rng):
        kappa = beta / 2.0
        shp = U.shape[1:5]
        ns = int(np.prod(shp))
        par = sum(np.indices(shp)[i] for i in range(4)) % 2
        idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
        Uf = U.reshape(4, ns, 2, 2)
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = staple(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                Uf[mu, ids] = su2_heatbath(kappa, A @ S, 0, 1, rng) @ A
        return U

    # ---- three 0++ operators (spatial averages, one value per time slice)
    def ops(U):
        # o0 : S_geo  = 12 (1 - Pbar)   (plaquette)
        # o1 : 1/4 tr F^2 clover form   (u08_corr.py:50-60, corpus-certified)
        # o2 : 1x2 rectangular Wilson loop, real part summed
        tot = 0.0
        clov = 0.0
        rect = 0.0
        for mu in range(4):
            for nu in range(mu + 1, 4):
                Pp = plaq(U, mu, nu)
                tot = tot + 0.5 * np.real(np.trace(Pp, axis1=-2, axis2=-1))
                C = 0.25 * (Pp + adj(np.roll(Pp, 1, axis=nu)) + adj(np.roll(Pp, 1, axis=mu))
                            + np.roll(np.roll(Pp, 1, axis=mu), 1, axis=nu))
                ca = 0.5 * np.imag(np.einsum('aij,...ji->...a', SIG, C))
                clov = clov + np.sum(ca * ca, axis=-1)
                # 1x2 in the (mu,nu) plane: U_mu(x) U_mu(x+mu) U_nu(x+2mu) U_mu^dag ... closed
                U1 = U[mu]
                U2 = np.roll(U[mu], -1, axis=mu)
                V1 = np.roll(U[nu], -1, axis=mu)
                V2 = np.roll(U[nu], 1, axis=mu)
                Wr = (U1 @ U2 @ V1 @ adj(np.roll(U2, -1, axis=nu))
                      @ adj(np.roll(U1, -1, axis=nu)) @ adj(V2))
                rect = rect + np.real(np.trace(Wr, axis1=-2, axis2=-1)) / 2.0
        # per-time-slice spatial averages, shape (L0,)
        o0 = 12.0 * (1.0 - (tot / 6.0).mean(axis=(1, 2, 3)))
        o1 = clov.mean(axis=(1, 2, 3))
        o2 = rect.mean(axis=(1, 2, 3))
        return np.stack([o0, o1, o2], axis=1), float((tot / 6.0).mean())

    rng = np.random.default_rng(args.seed)
    U = W.random_su2(rng, (4, args.L0, args.L, args.L, args.L))
    t0 = time.time()
    for _ in range(args.therm):
        update(U, args.beta, rng)
    NO = 3
    vals = np.empty((args.n, args.L0, NO))
    plaqm = np.empty(args.n)
    for k in range(args.n):
        for _ in range(args.gap):
            update(U, args.beta, rng)
        v, pl = ops(U)
        vals[k] = v
        plaqm[k] = pl
    el = time.time() - t0
    P(f"# TA03 MC SU(2) beta={args.beta}  {args.L0}x{args.L}^3  ncfg={args.n} "
      f"therm={args.therm} gap={args.gap} seed={args.seed}  [{el:.1f}s]")
    P(f"# <1/2 Re tr U_p> = {plaqm.mean():.6f} +- {plaqm.std(ddof=1)/np.sqrt(args.n):.6f}")

    def cmat(v):
        """C_ij(t) = (1/L0) sum_t0 [<Oi(t0+t)Oj(t0)> - <Oi><Oj>], t = 0,1"""
        n, L0, k = v.shape
        X = v - v.mean(axis=0, keepdims=True)
        C0 = np.zeros((k, k))
        C1 = np.zeros((k, k))
        for t0 in range(L0):
            a = X[:, t0, :]
            b = X[:, t0, :]
            c = X[:, (t0 + 1) % L0, :]
            C0 += a.T @ b
            C1 += c.T @ b
        return C0 / (n * L0), C1 / (n * L0)

    def masses(C0, C1):
        ev, V = np.linalg.eigh(C0)
        if ev.min() <= 0:
            return None
        A = V @ np.diag(ev ** -0.5) @ V.T          # C0^{-1/2}
        M = A @ C1 @ A
        lam = np.linalg.eigvalsh(M)[::-1]           # descending
        ms = []
        for L_ in lam:
            if L_ <= 0 or L_ >= 1:
                ms.append(np.nan)
            else:
                ms.append(float(np.arccosh(1.0 / L_)))
        return lam, np.array(ms), ev

    names = ['S_geo=12(1-Pbar)', '1/4 trF^2 clover', '1x2 rect Wilson']
    C0, C1 = cmat(vals)
    np.set_printoptions(precision=6, suppress=True)
    P("")
    P("C(0) matrix (0++ channel, per time slice):")
    P(str(C0))
    P("C(1) matrix:")
    P(str(C1))
    P(f"eig C(0) = {np.linalg.eigvalsh(C0)}  (need all > 0)")
    P("")
    P("single-operator T=2 effective mass  m = arccosh(C0_ii/C1_ii),  rho2 = m^2/m^2 (n/a):")
    for i in range(NO):
        if C1[i, i] > 0:
            P(f"   {names[i]:<22} C0/C1 = {C0[i,i]/C1[i,i]:8.5f}  "
              f"m_eff = {np.arccosh(C0[i,i]/C1[i,i]):8.5f}")
        else:
            P(f"   {names[i]:<22} C1_ii <= 0 : undefined")
    P("")
    P("2x2 / 3x3 GEVP:  C(1)v = lambda C(0)v ,  lambda_n = 1/cosh(m_n)")
    for pair in [(0, 1), (0, 2), (1, 2), (0, 1, 2)]:
        idx = list(pair)
        out = masses(C0[np.ix_(idx, idx)], C1[np.ix_(idx, idx)])
        lbl = "+".join(names[i].split()[0] for i in idx)
        if out is None:
            P(f"   [{lbl}] C(0) not positive definite -> GEVP undefined")
            continue
        lam, ms, ev = out
        P(f"   [{lbl}]  lambda = {lam}")
        P(f"        m = {ms}")
        if len(ms) >= 2 and np.all(np.isfinite(ms[:2])) and ms[0] > 0:
            r2 = (ms[1] / ms[0]) ** 2
            P(f"        rho^2 = (m2/m1)^2 = {r2:.6f}")
        else:
            P("        rho^2 undefined (a non-positive or >=1 GEVP eigenvalue)")

    # jackknife on the 2x2 (0,1) and (0,1,2)
    def jack(nb, idx):
        n = vals.shape[0]
        b = n // nb
        acc = []
        for i in range(nb):
            mask = np.ones(n, bool)
            mask[i * b:(i + 1) * b] = False
            C0b, C1b = cmat(vals[mask])
            o = masses(C0b[np.ix_(idx, idx)], C1b[np.ix_(idx, idx)])
            if o is None:
                acc.append((np.nan, np.nan))
            else:
                acc.append((o[1][0], o[1][1]))
        a = np.array(acc, float)
        return a
    P("")
    P("jackknife (8 bins) on the two lowest GEVP masses:")
    for pair in [(0, 1), (0, 2), (1, 2), (0, 1, 2)]:
        idx = list(pair)
        a = jack(8, idx)
        lbl = "+".join(names[i].split()[0] for i in idx)
        m1, m2 = a[:, 0], a[:, 1]
        ok = np.isfinite(m1) & np.isfinite(m2)
        if ok.sum() < 4:
            P(f"   [{lbl}] only {int(ok.sum())}/8 bins usable -> no error estimate")
            continue
        nb = ok.sum()
        e1 = np.sqrt((nb - 1) / nb * np.sum((m1[ok] - m1[ok].mean()) ** 2))
        e2 = np.sqrt((nb - 1) / nb * np.sum((m2[ok] - m2[ok].mean()) ** 2))
        r = (m2[ok] / m1[ok]) ** 2
        er = np.sqrt((nb - 1) / nb * np.sum((r - r.mean()) ** 2))
        P(f"   [{lbl}] m1 = {m1[ok].mean():.4f} +- {e1:.4f} ; "
          f"m2 = {m2[ok].mean():.4f} +- {e2:.4f} ; "
          f"rho^2 = {r.mean():.4f} +- {er:.4f}  [{nb}/8 bins]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L0', type=int, default=2)
    ap.add_argument('--L', type=int, default=2)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--n', type=int, default=40000)
    ap.add_argument('--therm', type=int, default=1000)
    ap.add_argument('--gap', type=int, default=2)
    ap.add_argument('--seed', type=int, default=2026)
    a = ap.parse_args()

    block_exact()
    block_mc(a)

    with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("[log written]")


if __name__ == '__main__':
    main()
