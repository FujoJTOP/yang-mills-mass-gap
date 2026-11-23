"""
WEYL-sign: numerical sign test of the Weyl-modulus marginal curvature (IDEA-new).

SU(2) Wilson theory.  beta=2.5.  We test whether the Vandermonde/Jacobian
curvature of the eigenvalue marginal makes V_eff'' > 0 where the raw action
(link) Hessian is negative, and whether the discriminant locus dominates.

Stages
  A  exact single-plaquette SU(2) modulus marginal (analytic + numeric Hessian,
     radial vs transverse decomposition of the link Hessian).
  B  lattice SU(2) plaquette-angle marginal: empirical mu(theta), log-curvature,
     discriminant-locus weight; raw and APE-smeared; L=4 and L=6.
  C  raw link Hessian (per-link + full gauge-projected transverse) vs the
     pointwise modulus marginal curvature.
"""
import os, sys, time, argparse
import numpy as np
import scipy.ndimage as ndi
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
CPB = os.path.abspath(os.path.join(HERE, '..', 'lattice-cpb'))
sys.path.insert(0, CPB)
import w4_blocked as W
from su3 import su2_heatbath, ID2, SIG
T2 = SIG / 2.0


# ---------------------------------------------------------------- helpers
def plaquette_angle_field(U):
    """Return theta_p = arccos(0.5 Re tr U_p) for all (mu<nu) plaquettes;
    flat array of angles in [0,pi]."""
    th = []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ W.adj(np.roll(U[mu], -1, axis=nu)) @ W.adj(U[nu]))
            tr = np.real(np.trace(Up, axis1=-2, axis2=-1)) / 2.0
            th.append(np.arccos(np.clip(tr, -1.0, 1.0)).ravel())
    return np.concatenate(th)


def mean_plaq(U):
    L = U.shape[1]
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ W.adj(np.roll(U[mu], -1, axis=nu)) @ W.adj(U[nu]))
            tot += np.real(np.trace(Up, axis1=-2, axis2=-1)).sum() / 2.0
    return tot / (6 * L ** 4)


# ===================================================================
# STAGE A : exact single-plaquette modulus marginal + radial/transverse
# ===================================================================
def stageA(beta):
    print("\n==== STAGE A : exact single-plaquette SU(2) modulus marginal ====")
    th = np.linspace(1e-5, np.pi - 1e-5, 4_000_001)
    S = beta * (1 - np.cos(th))
    Spp = beta * np.cos(th)                    # raw action curvature (theta)
    J = (2.0 * np.sin(th)) ** 2                # Vandermonde |lambda1-lambda2|^2
    Veff = S - np.log(J)
    Veffpp = beta * np.cos(th) + 2.0 / np.sin(th) ** 2
    print(f"  beta = {beta}")
    print(f"  raw S''(theta) = beta cos(theta):  min = {Spp.min():+.4f} (at theta=pi), "
          f"max = {Spp.max():+.4f};  S''<0 for theta>pi/2")
    i = np.argmin(Veffpp)
    print(f"  V_eff'' = beta cos + 2 csc^2:  GLOBAL MIN = {Veffpp[i]:+.4f} "
          f"at theta={th[i]:.4f} rad ({np.degrees(th[i]):.1f} deg), "
          f"cos={np.cos(th[i]):+.4f}, spread s=2sin={2*np.sin(th[i]):.4f}")
    neg = Spp < 0
    print(f"  on the Savvidy branch S''<0:  min V_eff'' = {Veffpp[neg].min():+.4f}  "
          f"(i.e. V_eff''>0 EVERYWHERE S''<0)")
    j = np.argmin(np.abs(Spp + 0.165))
    print(f"  at the VAR-FINAL proxy S''=-0.165: theta={th[j]:.4f}, "
          f"V_eff'' = {Veffpp[j]:+.4f}")
    csp = (np.sqrt(1 + beta ** 2) - 1) / beta
    thp = np.arccos(csp)
    print(f"  marginal mu peak: cos(theta*)={csp:.4f}, theta*={thp:.4f}, "
          f"spread s*={2*np.sin(thp):.4f},  V_eff''(peak)={beta*csp+2/(1-csp**2):+.4f}")

    # ---- raw link Hessian of one plaquette in link-angle phi, along modulus ----
    # U_p = U4^dag with U1=U2=U3=I ; vary U4 = exp(i eps.t^a) U4^0.
    # U4^0 = diag(e^{i phi/2}, e^{-i phi/2}) -> field-strength angle THETA=phi/2.
    # Raw action curvature in link coords:  S''(phi)=beta/4 cos(THETA).
    # Modulus-marginal curvature in link coords: (1/4)(beta cos THETA + 2 csc^2 THETA).
    print("\n  -- raw single-plaquette link Hessian along the modulus (link angle phi, "
          "Theta=phi/2):")
    for phi0 in [0.8, 1.8, np.pi + 0.1, 4.0, 5.0, 2 * np.pi - 0.2]:
        TH = phi0 / 2.0
        U40 = np.array([[np.exp(1j * phi0 / 2), 0], [0, np.exp(-1j * phi0 / 2)]], dtype=complex)
        h = 1e-5
        def Sc(eps):
            U4 = _exp2(np.array(eps)) @ U40
            Up = W.adj(U4)                          # plaquette (others identity)
            return beta * (1 - 0.5 * np.real(np.trace(Up)))
        e0 = np.zeros(3); H = np.zeros((3, 3))
        for a in range(3):
            for b in range(3):
                ep = e0.copy(); ep[a] += h; ep[b] += h
                ea = e0.copy(); ea[a] += h
                eb = e0.copy(); eb[b] += h
                H[a, b] = (Sc(ep) - Sc(ea) - Sc(eb) + Sc(e0)) / h ** 2
        H = 0.5 * (H + H.T)
        w = np.linalg.eigvalsh(H)
        Smarg = 0.25 * (beta * np.cos(TH) + 2.0 / np.sin(TH) ** 2)  # link-coord marginal
        print(f"    phi={phi0:5.3f} Theta={TH:5.3f} cos={np.cos(TH):+.3f}  "
              f"RAW H_modulus(link)={H[2,2]:+.4f} [beta/4 cosTh={beta/4*np.cos(TH):+.4f}]  "
              f"MARGINAL V_eff''(link)={Smarg:+.4f}  "
              f"eig=({w[0]:+.3f},{w[1]:+.3f},{w[2]:+.3f})")
    return


def _exp2(v):
    """exp(i v.t^a) for SU(2), v in R^3 (t^a = sigma^a/2)."""
    n = np.linalg.norm(v)
    if n < 1e-15:
        return ID2.copy()
    vh = v / n
    c = np.cos(n / 2.0); s = np.sin(n / 2.0)
    return c * ID2 + 1j * s * (vh[0] * SIG[0] + vh[1] * SIG[1] + vh[2] * SIG[2])


# ===================================================================
# STAGE B : lattice plaquette-angle marginal
# ===================================================================
def collect_angles(L, beta, seed, ntherm, nkeep, gapevery, ape=False, nsm=1):
    rng = np.random.default_rng(seed)
    cfgs = W.sample2(L, beta, rng, ntherm, nkeep, gapevery)
    allth = []
    for U in cfgs:
        if ape:
            U = W.ape_all2(U, alpha=0.5, nsteps=nsm)
        allth.append(plaquette_angle_field(U))
    return np.concatenate(allth), cfgs


def marginal_curvature(th, logdrop=3.5, nband=6):
    """Empirical mu(theta) via KDE (no bin-edge artifacts) -> curvature of -ln mu
    on the well-sampled support (density > peak*e^{-logdrop}); plus coarse
    curvature profile in nband equal-count bands."""
    kde = gaussian_kde(th)
    g = np.linspace(0.02, np.pi - 0.02, 4000)
    dens = kde(g)
    good = dens > dens.max() * np.exp(-logdrop)
    dg = g[1] - g[0]
    neglog = -np.log(np.maximum(dens, 1e-300))
    curv = np.gradient(np.gradient(neglog, dg), dg)
    curv_g = curv[good]; gg = g[good]
    imin = np.argmin(curv_g)
    # coarse profile: curvature averaged within each band (by density-decent range)
    lo, hi = gg.min(), gg.max()
    edges = np.linspace(lo, hi, nband + 1)
    prof = []
    for k in range(nband):
        sel = (gg >= edges[k]) & (gg < edges[k + 1])
        if sel.sum() > 3:
            prof.append((0.5 * (edges[k] + edges[k + 1]), float(curv_g[sel].mean())))
    return dict(grid=g, mu=dens, neglog=neglog, curv=curv, good=good,
                min_curv=float(curv_g.min()), at_theta=float(gg[imin]),
                peak_theta=float(g[np.argmax(dens)]),
                frac_negcurv=float(np.mean(curv_g < 0)), ngood=int(good.sum()),
                profile=prof)


def disc_weights(th):
    s = 2.0 * np.sin(th)                       # eigenvalue spread in [0,2]
    return dict(
        frac_s_lt_0p3=float(np.mean(s < 0.3)),
        frac_s_lt_0p5=float(np.mean(s < 0.5)),
        frac_theta_lt_0p15=float(np.mean(th < 0.15)),
        frac_theta_gt_pi_m0p15=float(np.mean(th > np.pi - 0.15)),
        median_s=float(np.median(s)),
        peak_theta=float(th[np.argmax(np.histogram(th, bins=90, range=(0, np.pi))[0])] +
                         (np.pi / 90) / 2),
    )


def stageB(beta, cases):
    print("\n==== STAGE B : lattice plaquette-angle marginal (SU(2)) ====")
    out = {}
    for (L, ape, nsm, tag) in cases:
        t0 = time.time()
        th, cfgs = collect_angles(L, beta, seed=100 + L, ntherm=250,
                                  nkeep=24 if L == 4 else 10,
                                  gapevery=4, ape=ape, nsm=nsm)
        P = np.mean([mean_plaq(U) for U in cfgs])
        mc = marginal_curvature(th)
        dw = disc_weights(th)
        print(f"\n  [{tag}] L={L} n={th.size} meanP={P:.4f} peak_theta={mc['peak_theta']:.3f} "
              f"({time.time()-t0:.1f}s)")
        print(f"     -ln mu curvature (KDE): min={mc['min_curv']:+.3f} at theta={mc['at_theta']:.3f}; "
              f"frac(curv<0)={mc['frac_negcurv']:.3f} over support "
              f"({'LOG-CONCAVE' if mc['min_curv']>=0 else 'NON-CONCAVE in some band'})")
        print("       curvature profile (theta:mean curv): " +
              "  ".join(f"{t:.2f}:{c:+.2f}" for t, c in mc['profile']))
        print(f"     discriminant weight: frac(s<0.3)={dw['frac_s_lt_0p3']:.4f} "
              f"frac(s<0.5)={dw['frac_s_lt_0p5']:.4f} "
              f"frac(theta<0.15)={dw['frac_theta_lt_0p15']:.4f} "
              f"frac(theta>pi-0.15)={dw['frac_theta_gt_pi_m0p15']:.4f} "
              f"median_s={dw['median_s']:.3f}")
        # pointwise modulus-marginal curvature on THIS theta distribution
        c = np.cos(th); s2 = np.sin(th) ** 2
        Veffpp = beta * c + 2.0 / s2
        Veffpp = Veffpp[np.isfinite(Veffpp)]
        sav = c < 0
        print(f"     pointwise V_eff''=beta cos+2csc^2 over sampled theta: "
              f"min={Veffpp.min():+.3f} frac<0={np.mean(Veffpp<0):.4f}; "
              f"Savvidy(S''<0) frac={np.mean(sav):.4f}, "
              f"min V_eff'' on Savvidy={Veffpp[sav].min() if sav.any() else float('nan'):+.3f}")
        out[tag] = dict(P=P, mc=mc, dw=dw)
    return out


# ===================================================================
# STAGE C : raw link Hessian vs modulus marginal
# ===================================================================
def stageC(beta, ncfg=3, L=4):
    print("\n==== STAGE C : raw link Hessian vs modulus marginal (L=%d) ====" % L)
    rng = np.random.default_rng(202)
    cfgs = W.sample2(L, beta, rng, 250, ncfg, 5)
    for i, U in enumerate(cfgs):
        t0 = time.time()
        pl = W.perlink_min2(U, beta)
        # full intrinsic Hessian + gauge projection
        H = W.hessian2(U, beta)
        Q = W.gauge_basis2(U)
        Z = W.transverse_basis(Q)
        Ht = Z.T @ H @ Z
        Ht = 0.5 * (Ht + Ht.T)
        ev = np.linalg.eigvalsh(Ht)
        th = plaquette_angle_field(U)
        c = np.cos(th); s2 = np.sin(th) ** 2
        Vpp = beta * c + 2.0 / s2
        print(f"\n  cfg{i}  meanP={mean_plaq(U):.4f}  ({time.time()-t0:.0f}s)")
        print(f"    RAW per-link Hessian: neg_frac={np.mean(pl<0):.4f} "
              f"min/beta={pl.min()/beta:+.4f}")
        print(f"    RAW full transverse Hessian (gauge-projected): dim={ev.size} "
              f"neg_frac={np.mean(ev<0):.4f} min/beta={ev.min()/beta:+.4f} "
              f"min={ev.min():+.3f}")
        print(f"    modulus marginal on same config: pointwise V_eff'' "
              f"min={Vpp.min():+.3f} frac<0={np.mean(Vpp<0):.4f}; "
              f"plaquette theta mean={th.mean():.3f} "
              f"(Savvidy frac={np.mean(c<0):.4f})")


# ===================================================================
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', default='all')
    ap.add_argument('--beta', type=float, default=2.5)
    a = ap.parse_args()
    if a.stage in ('all', 'A'):
        stageA(a.beta)
    if a.stage in ('all', 'B'):
        cases = [(4, False, 0, 'L4 raw'),
                 (4, True, 1, 'L4 APE1'),
                 (6, False, 0, 'L6 raw'),
                 (6, True, 1, 'L6 APE1'),
                 (6, True, 4, 'L6 APE4')]
        stageB(a.beta, cases)
    if a.stage in ('all', 'C'):
        stageC(a.beta, ncfg=3, L=4)
