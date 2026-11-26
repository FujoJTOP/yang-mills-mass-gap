"""
A1-RHO-L4 : can rho^2 = mu2/mu1 be *measured* on the L=4 lattice at beta=2.5 ?

rho^2 := mu2/mu1 = e2/e1 = Lambda*/m1^2   (ch05-theorems.tex:155; ch16:527-533; app:705)
criterion:  rho^2 <= 4  <=>  m2 <= 2 m1  <=>  the 2nd 0++ is below the 2-glueball cut
            <=>  the 2nd level is an ATOM       (atom3; DT-DIRICHLET:53)

ta03 (A03-RHO2.md) on L0=2: four legal pairings give rho^2 = 2.48 ... 4.45 (2.41 ... 5.79 across
seeds), i.e. the *scheme spread* (1.8x) dwarfs the statistical error bars.  On L=4 the corpus
(ch13-machine.tex:646-651) reports the single legal point (t0,t)=(1,2) gives rho^2 = 41.9 +- 46.0
with a 3-operator fit returning lambda2 <= 0 (-0.011) => rho^2 undefined.

This script does ONE thing: build a larger 0++ basis (6 operators, incl. two genuinely new
species -- adjoint-representation plaquette and magnetic-only plaquette) on L=4, run one MC,
and ask whether the 2nd GEVP eigenvalue is > 0, and how big the scheme spread is with (a) more
operators and (b) higher t / higher GEVP level.

BLOCK 1  arithmetic / structural (no MC):
   - if the 2nd level is the two-glueball cut, m2 = 2 m1 EXACTLY => rho^2 = 4 EXACTLY
     (= the criterion's own boundary)  [T]
   - resource scaling for the L=4 job (memory, sweeps, wall time per config)
BLOCK 2  MC SU(2) Wilson beta=2.5, L0 x L^3 = 4 x 4^3, 6 0++ operators, GEVP at (t0,t).

Run: py -3.12 a1_rho_l4.py --L0 4 --L 4 --n 1200 --gap 1 --therm 200 --budget 300
"""
import sys
import time
import argparse
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W                    # noqa: E402
from su3 import su2_heatbath, SIG         # noqa: E402

adj = W.adj
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


# ------------------------------------------------------------------ geometry
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


def path_prod(U, steps):
    """ordered product of links along `steps` = [(dir, +-1), ...], periodic b.c.
    returns array of shape U.shape[1:] + (2,2)."""
    disp = [0, 0, 0, 0]
    Wm = np.broadcast_to(np.eye(2, dtype=complex),
                         U.shape[1:5] + (2, 2)).copy()
    for (d, s) in steps:
        arr = U[d] if s > 0 else adj(U[d])
        for ax in range(4):
            if disp[ax]:
                arr = np.roll(arr, -disp[ax], axis=ax)
        Wm = Wm @ arr
        disp[d] += s
    return Wm


def loop(U, mu, nu, a, b):
    steps = ([(mu, 1)] * a + [(nu, 1)] * b + [(mu, -1)] * a + [(nu, -1)] * b)
    return path_prod(U, steps)


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


# ------------------------------------------------------------------ operators
NOPS = 6
NAMES = [
    'S_geo=12(1-Pbar)',      # 0 : plaquette (all 6 planes)
    'clover 1/4 trF^2',      # 1 : corpus-certified pointwise-positive form
    '1x2 rect Wilson',       # 2 : rectangular loop, all planes, both orientations
    '2x2 square Wilson',     # 3 : larger loop -> different overlap on excited states
    'adjoint plaquette',     # 4 : 1 - (1/3) tr_adj(P) = generic higher-rep 0++
    'magnetic-only (i<j)',   # 5 : space-space plaquettes only (tr B^2 channel)
]


def measure(U):
    """returns (L0, NOPS) array: spatial average of each operator on each time slice."""
    shp = U.shape[1:5]
    L0 = shp[0]
    acc = np.zeros((NOPS,) + tuple(shp))
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Pp = plaq(U, mu, nu)
            tr = 0.5 * np.real(np.trace(Pp, axis1=-2, axis2=-1))       # (L0,L,L,L)
            acc[0] += 12.0 * (1.0 - tr) / 6.0
            if mu >= 1 and nu >= 1:
                acc[5] += 6.0 * (1.0 - tr) / 3.0
            acc[4] += (1.0 - (np.abs(np.trace(Pp, axis1=-2, axis2=-1)) ** 2 - 1.0) / 3.0) / 6.0
            # clover
            C = 0.25 * (Pp + adj(np.roll(Pp, 1, axis=nu)) + adj(np.roll(Pp, 1, axis=mu))
                        + np.roll(np.roll(Pp, 1, axis=mu), 1, axis=nu))
            ca = 0.5 * np.imag(np.einsum('aij,...ji->...a', SIG, C))
            acc[1] += np.sum(ca * ca, axis=-1) / 6.0
            # loops (both orientations in the plane)
            for (a, b) in ((1, 2), (2, 2)):
                w = loop(U, mu, nu, a, b)
                r = np.real(np.trace(w, axis1=-2, axis2=-1)) / 2.0
                acc[2 if (a, b) == (1, 2) else 3] += r / 12.0
                w2 = loop(U, nu, mu, a, b)
                r2 = np.real(np.trace(w2, axis1=-2, axis2=-1)) / 2.0
                acc[2 if (a, b) == (1, 2) else 3] += r2 / 12.0
    return acc.reshape(NOPS, L0, -1).mean(axis=2).T


# ------------------------------------------------------------------ analysis
def cmat(v):
    """C_ij(t) = (1/L0) sum_t0 [<Oi(t0+t)Oj(t0)> - <Oi><Oj>], t=0..L0-1 (periodic)."""
    n, L0, k = v.shape
    X = v - v.mean(axis=0, keepdims=True)
    Cs = np.zeros((L0, k, k))
    for t in range(L0):
        acc = 0.0
        for t0 in range(L0):
            acc = acc + X[:, (t0 + t) % L0, :].T @ X[:, t0, :]
        Cs[t] = acc / (n * L0)
    return Cs


def gevp_masses(Cs, t0, t, T):
    """GEVP C(t) v = lam C(t0) v ; lam = cosh(m(T/2-t)) / cosh(m(T/2-t0)).
    returns (lam, m) sorted by ascending m, or None."""
    C0, C1 = Cs[t0], Cs[t]
    ev, V = np.linalg.eigh(C0)
    if ev.min() <= 0:
        return None
    A = V @ np.diag(ev ** -0.5) @ V.T
    M = A @ C1 @ A
    lam = np.linalg.eigvalsh(M)[::-1]
    ms = []
    xx = np.linspace(1e-5, 12.0, 200000)
    f = np.cosh(xx * (T / 2.0 - t)) / np.cosh(xx * (T / 2.0 - t0))
    for L_ in lam:
        if L_ <= 0:
            ms.append(np.nan)
            continue
        ms.append(float(xx[np.argmin(np.abs(f - L_))]))
    return lam, np.array(ms)


def rho2(Cs, idx, t0, t, T):
    sub = Cs[:, idx, :][:, :, idx]
    o = gevp_masses(sub, t0, t, T)
    if o is None:
        return None, None
    lam, ms = o
    if len(ms) >= 2 and np.isfinite(ms[:2]).all() and ms[0] > 0:
        return (ms[1] / ms[0]) ** 2, ms
    return None, ms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L0', type=int, default=4)
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--n', type=int, default=1200)
    ap.add_argument('--therm', type=int, default=200)
    ap.add_argument('--gap', type=int, default=1)
    ap.add_argument('--seed', type=int, default=77)
    ap.add_argument('--budget', type=float, default=300.0, help='wall seconds for the MC loop')
    a = ap.parse_args()
    L0, L, beta = a.L0, a.L, a.beta
    V = L0 * L ** 3

    P("=" * 84)
    P("BLOCK 1  structural / resource (no MC)")
    P("=" * 84)
    P(f"lattice {L0} x {L}^3 = {V} sites, {4*V} links, SU(2) 3 d.o.f./link = {3*4*V} real d.o.f.")
    P(f"basis: {NOPS} 0++ operators  ->  GEVP matrices {NOPS}x{NOPS} = "
      f"{NOPS*NOPS*8*L0*8/1e3:.1f} kB per configuration set (nothing)")
    P("")
    P("[T] criterion boundary is kinematic, not dynamic: if the 2nd 0++ level is the")
    P("    two-glueball cut, then m2 = 2 m1 EXACTLY, hence rho^2 = 4 EXACTLY.")
    P("    So rho^2 <= 4 is satisfied *with equality* by a NON-atomic second level:")
    P("    a measurement returning rho^2 = 4.0 +- s cannot separate 'atom just below")
    P("    the cut' from 'the cut'.  Only m2 < 2m1 by >> s is evidence for an atom.")
    P("    => the L=4 verdict hinges on resolving m2 vs 2m1, i.e. on Delta = (rho^2-1)mu1,")
    P("       not on rho^2 itself (atom3: Delta=(rho^2-1)mu1 is the load-bearing object).")
    P("")

    # ---------------- MC (one job, adaptive length)
    rng = np.random.default_rng(a.seed)
    U = W.random_su2(rng, (4, L0, L, L, L))
    t_start = time.time()
    for _ in range(a.therm):
        update(U, beta, rng)
    t_th = time.time() - t_start

    # benchmark 8 sweeps to pick n within budget
    tb = time.time()
    for _ in range(8):
        update(U, beta, rng)
    t_sw = (time.time() - tb) / 8.0
    tb = time.time()
    for _ in range(4):
        measure(U)
    t_ms = (time.time() - tb) / 4.0
    per_cfg = a.gap * t_sw + t_ms
    n_fit = max(200, int(min(a.n, a.budget / per_cfg)))
    P(f"[bench] therm {a.therm} sw = {t_th:.1f}s ; sweep = {t_sw*1e3:.1f} ms ; "
      f"measure = {t_ms*1e3:.1f} ms ; per config (gap={a.gap}) = {per_cfg*1e3:.0f} ms")
    P(f"[bench] budget {a.budget:.0f}s => n = {n_fit} configs (requested {a.n}), "
      f"estimated {per_cfg*n_fit:.0f}s")
    P("")

    vals = np.empty((n_fit, L0, NOPS))
    t0 = time.time()
    for k in range(n_fit):
        for _ in range(a.gap):
            update(U, beta, rng)
        vals[k] = measure(U)
        if (k + 1) % 100 == 0:
            el = time.time() - t0
            P(f"    ... {k+1}/{n_fit} configs  [{el:.0f}s, "
              f"{el/(k+1)*1000:.0f} ms/cfg]")
    el = time.time() - t0
    P(f"# A1 MC SU(2) beta={beta}  {L0}x{L}^3  ncfg={n_fit} therm={a.therm} "
      f"gap={a.gap} seed={a.seed}  [{el:.1f}s]")
    plq = 1.0 - vals[:, :, 0].mean() / 12.0
    P(f"# <1/2 Re tr U_p> = {plq:.6f}   (sanity; SU(2) beta=2.5 should be ~0.63-0.67)")
    P("")

    T = L0
    Cs = cmat(vals)
    np.set_printoptions(precision=5, suppress=True)
    P("C(0) (connected, per time slice, 6x6):")
    P(str(Cs[0]))
    P(f"eig C(0) = {np.linalg.eigvalsh(Cs[0])}")
    for t in range(1, L0):
        P(f"diag C({t})/diag C(0) = {np.diag(Cs[t])/np.diag(Cs[0])}")
    P("")

    # ---------------- GEVP ladder: basis level x (t0,t)
    # operator priority = decreasing C(0) signal; large loops are noise-dominated at L=4
    PRIOR = [0, 1, 5, 4, 2, 3]
    P("single-operator T=4 correlator shape, per operator:")
    for i in range(NOPS):
        P(f"   {NAMES[i]:<22} C0={Cs[0][i,i]:+.3e}  "
          f"diag C(t)/C0 = " + " ".join(f"{Cs[t][i,i]/Cs[0][i,i]:+.4f}" for t in (1, 2, 3)))
    P("")
    P(f"operator priority order (by C(0) signal): "
      f"{[NAMES[i] for i in PRIOR]}")
    P("")
    P("GEVP ladder:  C(t) v = lam C(t0) v ,  lam = cosh(m(T/2-t))/cosh(m(T/2-t0)) ,  T=L0")
    P(f"{'k':>2} {'(t0,t)':>8} {'eig C(t0) min':>13} {'lam':>40} {'m':>26} {'rho^2':>9}")
    rows = []
    for k in range(2, NOPS + 1):
        idx = list(PRIOR[:k])
        for (t0p, tp) in ((1, 2), (1, 3), (2, 3)):
            if tp >= L0:
                continue
            sub = Cs[:, idx, :][:, :, idx]
            ev0 = np.linalg.eigvalsh(sub[t0p])
            o = gevp_masses(sub, t0p, tp, T)
            if o is None:
                P(f"{k:>2} {str((t0p,tp)):>8} {ev0.min():>13.2e}  "
                  f"C(t0) not positive definite -> undefined")
                continue
            lam, ms = o
            r2, _ = rho2(Cs, idx, t0p, tp, T)
            ls = " ".join(f"{x:+.4f}" for x in lam)
            mstr = " ".join(("nan" if not np.isfinite(x) else f"{x:.4f}") for x in ms)
            P(f"{k:>2} {str((t0p,tp)):>8} {ev0.min():>13.2e} {ls:>40} {mstr:>26} "
              f"{'undef' if r2 is None else f'{r2:9.4f}'}")
            if r2 is not None:
                rows.append((k, t0p, tp, r2, lam[1] if len(lam) > 1 else np.nan))
    P("")

    # explicit ta03 pairings, for direct comparison with the L0=2 numbers
    P("ta03 pairings on L=4 (same three operators as A03-RHO2.md:104-107):")
    for pair in [(0, 1), (0, 2), (1, 2), (0, 1, 2)]:
        for (t0p, tp) in ((1, 2), (1, 3)):
            r2, ms = rho2(Cs, list(pair), t0p, tp, T)
            lbl = "+".join(NAMES[i].split()[0] for i in pair)
            P(f"   [{lbl:<28}] (t0,t)=({t0p},{tp})  "
              f"{'rho^2 = %.4f' % r2 if r2 is not None else 'rho^2 undefined (lambda2<=0)'}")
    P("")

    if rows:
        arr = np.array([r[3] for r in rows])
        P(f"[E] rho^2 over {len(rows)} defined (basis-level, t0, t) schemes: "
          f"min {arr.min():.4f}  max {arr.max():.4f}  "
          f"spread {arr.max()/max(arr.min(),1e-9):.2f}x  ;  >4: {int((arr>4).sum())}/{len(arr)}")
        for k in range(2, NOPS + 1):
            sub = np.array([r[3] for r in rows if r[0] == k])
            if len(sub) >= 2:
                P(f"    k={k}: rho^2 in [{sub.min():.3f}, {sub.max():.3f}] "
                  f"({sub.max()/sub.min():.2f}x)")
    else:
        P("[E] every scheme returned lambda2 <= 0 or a non-positive metric: rho^2 undefined "
          "on this lattice even with a 6-operator basis.")
    P("")

    # ---------------- jackknife on the largest basis
    P("jackknife (8 bins), largest basis (all 6):")
    idx = list(PRIOR)
    for (t0p, tp) in ((1, 2), (1, 3), (2, 3)):
        if tp >= L0:
            continue
        nb = 8
        n = vals.shape[0]
        b = n // nb
        acc = []
        for i in range(nb):
            mask = np.ones(n, bool)
            mask[i * b:(i + 1) * b] = False
            Cb = cmat(vals[mask])
            r2, ms = rho2(Cb, idx, t0p, tp, T)
            acc.append((np.nan if r2 is None else r2,
                        ms[0] if ms is not None else np.nan,
                        ms[1] if (ms is not None and len(ms) > 1) else np.nan))
        A_ = np.array(acc, float)
        ok = np.isfinite(A_[:, 0])
        if ok.sum() < 4:
            P(f"   (t0,t)=({t0p},{tp}): only {int(ok.sum())}/8 bins give a defined rho^2 "
              f"-> no error estimate (this IS the L=4 obstruction)")
            continue
        nb2 = int(ok.sum())
        r = A_[ok, 0]
        er = np.sqrt((nb2 - 1) / nb2 * np.sum((r - r.mean()) ** 2))
        m1 = A_[ok, 1]
        e1 = np.sqrt((nb2 - 1) / nb2 * np.sum((m1 - m1.mean()) ** 2))
        P(f"   (t0,t)=({t0p},{tp}) k={NOPS}: rho^2 = {r.mean():.4f} +- {er:.4f} ; "
          f"m1 = {m1.mean():.4f} +- {e1:.4f} ({nb2}/8 bins)")
    P("")
    P("[E] statistics projection (sigma ~ 1/sqrt(n)): to decide rho^2 <= 4 you need")
    P("    |rho^2-4| >> sigma(stat) AND spread(scheme) << 1.  Statistics shrink the bars")
    P("    around each scheme; they cannot shrink the spread between schemes.")
    P("")
    P(f"[wall] total {time.time() - t_start:.1f}s")
    with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("[log written]")


if __name__ == '__main__':
    main()
