"""
U08 -- connected correlator C(t) = <O(t)O(0)>_c of O = 1/4 tr F^2, SU(2) Wilson
action, beta = 2.5, and the E1 estimate from C(t).

O on the lattice -- CORPUS-CERTIFIED form (LOOP-RUNBOOK.md:5377, U04-MEASF1.md:52-60):
    clover  C_{uv}(x) = 1/4 [ U_uv(x) + U_uv(x-nu)^dag + U_uv(x-mu)^dag + U_uv(x-mu-nu) ]
    O(x) = sum_{mu<nu} sum_a ( Re c_a(C_{uv}(x)) )^2           (a = 1, lattice units)
  with U = c0 + i c.sigma ,  Re c_a(C) = 1/2 Im tr(sigma_a C).
  => O(x) >= 0 pointwise, bounded.  NOT 1 - 1/2 Re tr U_p  (that is S_geo, not 1/4 tr F^2),
  NOT (U_p - U_p^dag)/(2 i a^2 g)  (covariant tensor component, wrong channel).
  --op plaq switches to 12(1-Pbar) for an operator-definition cross-check only.

Geometry (L0,L1,L2,L3) periodic, time axis = 0:
    Obar(t) = (1/V_s) sum_xs O(xs,t)
    C(t)    = (1/L0) sum_t0 [ <Obar(t0+t) Obar(t0)> - <O>^2 ]      =>  C(0) = Var(Obar)
  [asserted below]

Run:  py -3.12 u08_corr.py --L0 2 --L 2 --beta 2.5 --n 30000 --gap 2
"""
import sys
import time
import argparse
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W          # noqa: E402
from su3 import su2_heatbath, SIG  # noqa: E402

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


def o_field_clover(U):
    """O(x) = sum_{mu<nu} sum_a (Re c_a(C_uv))^2 , clover, a=1, >= 0 pointwise."""
    O = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            P = plaq(U, mu, nu)
            C = 0.25 * (P + adj(np.roll(P, 1, axis=nu)) + adj(np.roll(P, 1, axis=mu))
                        + np.roll(np.roll(P, 1, axis=mu), 1, axis=nu))
            ca = 0.5 * np.imag(np.einsum('aij,...ji->...a', SIG, C))
            O = O + np.sum(ca * ca, axis=-1)
    return O


def o_field_plaq(U):
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            tot = tot + 0.5 * np.real(np.trace(plaq(U, mu, nu), axis1=-2, axis2=-1))
    return 12.0 * (1.0 - tot / 6.0)


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


def run(L0, L, beta, ncfg, ntherm, gap, seed, op):
    rng = np.random.default_rng(seed)
    of = o_field_clover if op == 'clover' else o_field_plaq
    U = W.random_su2(rng, (4, L0, L, L, L))
    t0 = time.time()
    for _ in range(ntherm):
        update(U, beta, rng)
    vals = np.empty((ncfg, L0))
    plaqm = np.empty(ncfg)
    ox, ox2 = np.empty(ncfg), np.empty(ncfg)
    for k in range(ncfg):
        for _ in range(gap):
            update(U, beta, rng)
        O = of(U)
        vals[k] = O.mean(axis=(1, 2, 3))
        plaqm[k] = (0.5 * np.real(np.trace(plaq(U, 0, 1), axis1=-2, axis2=-1))).mean()
        ox[k], ox2[k] = O.mean(), (O ** 2).mean()
    return dict(vals=vals, plaq=plaqm, ox=ox, ox2=ox2, el=time.time() - t0)


def corr(vals):
    n, L0 = vals.shape
    X = vals - vals.mean()
    return np.array([np.mean(X * np.roll(X, -t, axis=1)) for t in range(L0)])


def jk_err(vals, t, nbin):
    n = vals.shape[0]
    b = n // nbin
    cs = np.array([corr(vals[i * b:(i + 1) * b])[t] for i in range(nbin)])
    return float(np.sqrt((nbin - 1) / nbin * np.sum((cs - cs.mean()) ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L0', type=int, default=2)
    ap.add_argument('--L', type=int, default=2)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--n', type=int, default=30000)
    ap.add_argument('--therm', type=int, default=1000)
    ap.add_argument('--gap', type=int, default=1)
    ap.add_argument('--seed', type=int, default=2026)
    ap.add_argument('--op', default='clover')
    a = ap.parse_args()

    r = run(a.L0, a.L, a.beta, a.n, a.therm, a.gap, a.seed, a.op)
    vals, L0, T = r['vals'], a.L0, a.L0
    C = corr(vals)
    err = np.array([jk_err(vals, t, 8) for t in range(L0)])
    var1 = float(np.var(vals[:, 0]))
    var_pw = float(r['ox2'].mean() - r['ox'].mean() ** 2)
    x = vals[:, 0] - vals[:, 0].mean()
    n = len(x)
    ac = np.correlate(x, x, 'full')[n - 1:] / (np.arange(n, 0, -1) * x.var())
    tau = 0.5 + float(np.sum(ac[1:min(n // 2, 2000)]))

    print(f"# U08 SU(2) beta={a.beta}  {L0}x{a.L}^3  op={a.op}  ncfg={a.n} "
          f"therm={a.therm} gap={a.gap} seed={a.seed}  [{r['el']:.1f}s]")
    print(f"# <O> = {vals.mean():.6f}   <1/2 Re tr U_p> = {r['plaq'].mean():.6f} "
          f"+- {r['plaq'].std(ddof=1)/np.sqrt(a.n):.6f}")
    print(f"# [T] C(0) = Var(Obar):  C(0)={C[0]:.10e}  "
          f"mean-of-slice-Var={var1:.6f}  (diff {abs(C[0]-var1):.2e} = MC, "
          f"identical by construction)")
    print(f"#     Var_pointwise(O(x)) = {var_pw:.6e}  != C(0)  "
          f"(ratio {var_pw/C[0]:.2f}, V_s={a.L**3})")
    print(f"# [E] tau_int={tau:.2f} meas-units")
    print("\n  t     C(t)/C(0)        C(t)             err")
    for t in range(L0):
        print(f"  {t:<4} {C[t]/C[0]:.8f}   {C[t]:+.8e}   {err[t]:.2e}")
    print("\n  t   m_eff=ln[C(t-1)/C(t)]   single-pole finite-T solve"
          "  [C(t)=f^2(e^{-mt}+e^{-m(T-t)})/(1+e^{-mT})]")
    ms = np.linspace(1e-4, 9.0, 400000)
    for t in range(1, L0):
        if C[t] <= 0:
            print(f"  {t:<4} (C(t)<=0)")
            continue
        mod = (np.exp(-ms * t) + np.exp(-ms * (T - t))) / (1.0 + np.exp(-ms * T))
        mm = float(ms[np.argmin(np.abs(mod - C[t] / C[0]))])
        print(f"  {t:<4} {np.log(C[t-1]/C[t]):.6f}              {mm:.6f}")
    print()
    if C[1] > 0:
        print(f"[E] E1(naive, ln C0/C1) = {np.log(C[0]/C[1]):.6f}")
    if T == 2 and C[1] > 0:
        print(f"[T] L_t=2 single pole: C(1)=f^2/cosh(E1) => E1 = arccosh(C0/C1) "
              f"= {np.arccosh(C[0]/C[1]):.6f}   (U04-MEASF1.md:108)")
    print(f"[wall] {r['el']:.1f}s")


if __name__ == '__main__':
    main()
