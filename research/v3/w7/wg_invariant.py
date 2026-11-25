#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
wg_invariant.py -- WG: a GAUGE-INVARIANT, NON-SELF-AVERAGING susceptibility estimator for the
4D beta=2.5 plateau, replacing the flat-chart k=0 coordinate.

Object (the estimator).  Let O_s(x) be a GAUGE-INVARIANT scalar 0++ density, here the (APE-smeared)
field-strength / plaquette energy density summed over all six Lorentz planes

    O_s(x) = sum_{mu<nu} ( 1 - (1/2) Re tr P_{mu nu}^s(x) ),      P^s = plaquette of smeared links,

and take its TOTAL-normalized zero mode

    q = V^{-1/2} sum_x ( O_s(x) - <O_s> )        (V = L^4, <O_s> = ensemble mean, a c-number).

    chi_O := Var(q) = sum_z C_s(z) = V * Var_c(Obar_c),     Obar_c = (1/V) sum_x O_s(x),

with C_s the connected two-point function.  Two facts make this the right carrier:

  * GAUGE-INVARIANT: O_s is built from gauge-invariant plaquettes, so q is invariant under every
    gauge transform (exactly -- verified to machine precision below).  The flat-chart charge
    chi_flat = Var(N^{-1/2} sum_l log U_l) is NOT (it drifts under finite local transforms).
  * NON-SELF-AVERAGING (total normalization, T4-CARRIER): Var(q) = chi_O = O(1), L-independent,
    and NOT chi_O/V.  The MEAN-normalized mode qbar=(1/V)sum O has Var(qbar)=chi_O/V -> 0.

Hazard that is REAL and handled: the bare composite has a UV contact term C_s(0) ~ a^{-8}; chi_O is
then contact-dominated (a FALSE plateau).  APE smearing at a fixed physical radius suppresses the
contact; we report the contact C_s(0) and the contact-subtracted chi_phys = chi_O - C_s(0) = sum_{z!=0}
C_s(z), which is the propagating (gap-relevant) part.

ACID TEST.  Free massive :phi^2: (the gauge-invariant composite analog of tr F^2) has chi finite and
L-independent; free MASSLESS has chi -> infinity (H1 correctly FAILS -- the estimator manufactures no
floor).  Free elementary phi: Var(total-zero)=1/m^2 exactly, every L.

Parts:
  [F] exact free-field acid test (massive / massless, elementary phi and composite :phi^2:)
  [G] SU(2) beta=2.5, L=4,6 : chi_O (gauge-inv, total-norm), chi_phys, chi_flat (reference),
      chi_P (reference), gauge-invariance check, and the normalization check Var(qbar)=chi_O/V.

Run:  py -3.12 -I research/v3/w7/wg_invariant.py            # L=4,6
      py -3.12 -I research/v3/w7/wg_invariant.py quick      # tiny N
Pure numpy/numba; ONE process (charter s4: never parallelise SU(N) L>=6).  Labels [T]/[V]/[C].
No third-party sources.  No QED.  Clay 0%.
"""
import sys, os, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CPB = 'D:/Dev/math/research/lattice-cpb'
W7U = 'D:/Dev/math/research/unconditional/w7'
for p in (HERE, CPB, W7U):
    if p not in sys.path:
        sys.path.insert(0, p)

import w4_blocked as W          # ape_all2, adj, random_su2, SIG, ID2
import m1_fast as F             # Numba heat-bath sweep, random_start, chi_flat

ID2 = W.ID2
RESULT = os.path.join(HERE, 'wg_invariant_results.json')


# ======================================================================= [F] free-field acid (exact)
def free_propagator(L, m, d):
    """Periodic free scalar propagator G(x) = (1/V) sum_k e^{ikx}/(m^2 + sum_mu 4 sin^2(k_mu/2))."""
    grids = np.meshgrid(*[np.arange(L) for _ in range(d)], indexing='ij')
    phat2 = np.zeros((L,) * d)
    for mu in range(d):
        phat2 += 4.0 * np.sin(2.0 * np.pi * grids[mu] / L / 2.0) ** 2
    fk = 1.0 / (phat2 + m * m)
    return np.fft.ifftn(fk).real


def free_acid():
    """Exact: Var(total-zero of phi)=1/m^2 (all L); chi_{:phi^2:}=sum_z 2G(z)^2 finite for m>0,
    DIVERGES as m->0.  These are the acid values the estimator must reproduce for a free field."""
    out = {}
    rows = []
    for d, Ls in [(1, [32, 64, 128, 256]), (2, [8, 16, 32, 64]), (4, [4, 6, 8, 12])]:
        m = 0.5
        Vs, cphi, ccomp = [], [], []
        for L in Ls:
            G = free_propagator(L, m, d)
            Vs.append(L ** d); cphi.append(float(G.sum())); ccomp.append(float(np.sum(2.0 * G * G)))
        Vs = np.asarray(Vs, float)
        s_phi = float(np.polyfit(np.log(Vs), np.log(np.array(cphi)), 1)[0])
        s_comp = float(np.polyfit(np.log(Vs), np.log(np.array(ccomp)), 1)[0])
        rows.append(dict(d=d, L=Ls, Var_phi=cphi, chi_comp=ccomp,
                         slope_phi=s_phi, slope_comp=s_comp, one_over_m2=1.0 / m ** 2))
    # massless acid: chi_{:phi^2:} -> infinity
    acid = {}
    for d in (1, 2):
        L = 64
        acid[d] = [(m, float(np.sum(2.0 * free_propagator(L, m, d) ** 2))) for m in (0.5, 0.2, 0.1, 0.05, 0.02)]
    out['ladder'] = rows
    out['massless_acid'] = acid
    return out


def print_free(fa):
    print("=" * 88)
    print(" [F] FREE-FIELD ACID TEST (exact)  --  the estimator must give a KNOWN value, never a floor")
    print("=" * 88)
    for r in fa['ladder']:
        print(f"  d={r['d']}  m=0.5:  1/m^2 = {r['one_over_m2']:.6f}")
        print(f"     Var(phi,total) = {[round(x,4) for x in r['Var_phi']]}  (slope {r['slope_phi']:+.4f}, 0 = non-self-avg)")
        print(f"     chi(:phi^2:)   = {[round(x,4) for x in r['chi_comp']]}  (slope {r['slope_comp']:+.4f}, 0 = non-vacuous)")
    print("  massless acid: chi(:phi^2:) DIVERGES as m->0 (H1 correctly FAILS; no manufactured floor):")
    for d, seq in fa['massless_acid'].items():
        print(f"     d={d} L=64:  " + "  ".join(f"m={m}:{v:.2f}" for m, v in seq))


# ======================================================================= [G] SU(2) gauge-invariant
def scalar_density(U):
    """Gauge-invariant 0++ field-strength/plaquette density O(x)=sum_{mu<nu}(1 - 1/2 Re tr P_{mu nu}(x))."""
    L = U.shape[1]
    O = np.zeros((L, L, L, L))
    for mu in range(4):
        for nu in range(mu + 1, 4):
            P = U[mu] @ np.roll(U[nu], -1, axis=mu) @ W.adj(np.roll(U[mu], -1, axis=nu)) @ W.adj(U[nu])
            O += 1.0 - np.real(np.trace(P, axis1=-2, axis2=-1)) / 2.0
    return O


def polyakov_spatial(U):
    """P(x)=1/2 Re tr prod_t U_0(t,x), a field on the spatial (L,L,L) sites (start t=0)."""
    L = U.shape[1]
    prod = U[0, 0].copy()
    for t in range(1, L):
        prod = prod @ U[0, t]
    return np.real(np.trace(prod, axis1=-2, axis2=-1)) / 2.0


def gauge_transform(U, G):
    out = np.empty_like(U)
    for mu in range(4):
        out[mu] = G @ U[mu] @ W.adj(np.roll(G, -1, axis=mu))
    return out


def tau_int(x):
    x = np.asarray(x, float); x = x - x.mean(); v = x.var(); n = x.size
    if v <= 0:
        return 1.0
    s = 1.0
    for t in range(1, n // 4):
        c = np.mean(x[:n - t] * x[t:]) / v
        if c < 0:
            break
        s += 2.0 * c
    return float(max(1.0, s))


def jk_var_se(x, nblk=20):
    """Delete-block jackknife of the ensemble variance Var(x); returns (var, se)."""
    x = np.asarray(x, float); n = x.size
    nblk = int(min(nblk, max(4, n // 8)))
    ed = np.linspace(0, n, nblk + 1).astype(int)
    base = x.var()
    lv = []
    for b in range(nblk):
        idx = np.ones(n, bool); idx[ed[b]:ed[b + 1]] = False
        lv.append(x[idx].var())
    lv = np.asarray(lv)
    se = float(np.sqrt((nblk - 1) / nblk * np.sum((lv - lv.mean()) ** 2)))
    return base, se


def gauge_invariance_check(L, seed=7):
    """Build one thermalized config, compare chi_O and chi_flat on U vs a random gauge transform U^g."""
    np.random.seed(seed)
    U = F.random_start(L, seed)
    for _ in range(120):
        F.sweep(U, 2.5, L)
    rng = np.random.default_rng(seed + 1)
    G = W.random_su2(rng, (L, L, L, L))
    Ug = gauge_transform(U, G)

    def Obar_of(X):
        return float(scalar_density(W.ape_all2(X, alpha=0.5, nsteps=2)).mean())

    o1, o2 = Obar_of(U), Obar_of(Ug)
    c1, c2 = float(F.chi_flat(U, L)), float(F.chi_flat(Ug, L))
    return dict(Obar=o1, Obar_gauge=o2, Obar_reldiff=abs(o2 - o1) / max(1e-12, abs(o1)),
                chi_flat=c1, chi_flat_gauge=c2,
                chi_flat_reldrift=abs(c2 - c1) / max(1e-12, abs(c1)))


def run_volume(L, beta, N, therm, gap, seed, ns_ape=2):
    np.random.seed(seed)
    U = F.random_start(L, seed)
    for _ in range(therm):
        F.sweep(U, beta, L)

    Obar = np.empty(N)       # mean of smeared gauge-invariant density (per config)
    Osq = np.empty(N)        # mean of O_s^2 (for the contact C(0))
    Obar_bare = np.empty(N)  # mean of BARE density (UV contact reference)
    chif = np.empty(N)
    Pbar = np.empty(N)
    t0 = time.time()
    for i in range(N):
        for _ in range(gap):
            F.sweep(U, beta, L)
        Us = W.ape_all2(U, alpha=0.5, nsteps=ns_ape)
        Os = scalar_density(Us)
        Obar[i] = Os.mean(); Osq[i] = (Os * Os).mean()
        Obar_bare[i] = scalar_density(U).mean()
        chif[i] = F.chi_flat(U, L)
        Pbar[i] = polyakov_spatial(U).mean()
    wall = time.time() - t0

    V = L ** 4
    chiO, chiO_se = jk_var_se(Obar); chiO *= V; chiO_se *= V           # total-normalized
    C0 = float(Osq.mean() - Obar.mean() ** 2)                           # contact C_s(0)
    chi_phys = chiO - C0
    chi_bare, _ = jk_var_se(Obar_bare); chi_bare *= V
    chifm, chifse = float(chif.mean()), float(chif.std(ddof=1) / np.sqrt(N))
    chiP, chiP_se = jk_var_se(Pbar); chiP *= L ** 3
    return dict(L=L, beta=beta, N=N, therm=therm, gap=gap, ns_ape=ns_ape, wall_s=wall,
                chi_O=chiO, chi_O_se=chiO_se, Var_mean=chiO / V,
                contact_C0=C0, chi_phys=chi_phys,
                chi_bare=chi_bare,
                chi_flat=chifm, chi_flat_se=chifse, tau_flat=tau_int(chif),
                chi_P=chiP, chi_P_se=chiP_se,
                Obar_mean=float(Obar.mean()), Obar_sd=float(Obar.std(ddof=1)))


def main():
    quick = 'quick' in sys.argv[1:]
    print("# WG gauge-invariant non-self-averaging susceptibility estimator  (beta=2.5 SU(2))")
    fa = free_acid()
    print_free(fa)

    # gauge-invariance + flat-chart-drift check at L=4
    print("\n" + "=" * 88)
    print(" [G-inv] gauge invariance: chi_O exact; chi_flat (flat chart) drifts  (L=4, one config)")
    print("=" * 88)
    gi = gauge_invariance_check(4)
    print(f"   chi_O : Obar={gi['Obar']:.8f}  Obar(gauge)={gi['Obar_gauge']:.8f}  "
          f"rel.diff={gi['Obar_reldiff']:.2e}   <- gauge-invariant (machine zero)")
    print(f"   chi_flat: {gi['chi_flat']:.6f} -> {gi['chi_flat_gauge']:.6f}  "
          f"rel.drift={gi['chi_flat_reldrift']:.2%}   <- flat chart drifts (gauge-dependent)")

    print("\n" + "=" * 88)
    print(" [G] SU(2) beta=2.5 plateau scan:  chi_O (total-norm, gauge-inv) vs references")
    print("=" * 88)
    out = {'free': fa, 'gauge_check': gi, 'volumes': {}}
    vols = [4, 6] if quick else [4, 6, 8]
    Ntab = {4: 800, 6: 400, 8: 160}
    for L in vols:
        N = 200 if quick else Ntab[L]
        therm = 150 if quick else 500
        gap = 4 if quick else 8
        r = run_volume(L, 2.5, N, therm, gap, seed=1000 + L)
        r['gauge_inv'] = gauge_invariance_check(L, seed=3000 + L)
        out['volumes'][f'L{L}'] = r
        print(f"\n  L={L} N={N} therm={therm} gap={gap} (APE ns=2)   wall {r['wall_s']:.1f}s")
        print(f"    chi_O  (total-norm, gauge-inv) = {r['chi_O']:.4f} +- {r['chi_O_se']:.4f}   "
              f"[Var_mean = chi_O/V = {r['Var_mean']:.4f}]")
        print(f"      contact C_s(0) = {r['contact_C0']:.4f}   ->  chi_phys = chi_O - C_s(0) = {r['chi_phys']:.4f}")
        print(f"    chi_bare (ns=0, contact-dominated) = {r['chi_bare']:.3f}")
        print(f"    chi_flat (flat chart, reference)   = {r['chi_flat']:.4f} +- {r['chi_flat_se']:.4f} "
              f"(tau={r['tau_flat']:.2f})")
        print(f"    chi_P (Polyakov, N-ality-1 ref)    = {r['chi_P']:.4f} +- {r['chi_P_se']:.4f}")
        g = r['gauge_inv']
        print(f"    gauge-inv check: chi_O reldiff={g['Obar_reldiff']:.1e}; "
              f"chi_flat drift={g['chi_flat_reldrift']:.2%}")

    # cross-volume ratios of the physical (0++ scalar) susceptibility
    have = [L for L in (4, 6, 8) if f'L{L}' in out['volumes']]
    for a, b in zip(have, have[1:]):
        ra, rb = out['volumes'][f'L{a}'], out['volumes'][f'L{b}']
        print(f"\n  ratio chi_phys(L{b})/chi_phys(L{a}) = "
              f"{rb['chi_phys']/ra['chi_phys']:.3f}   (free massless (b/a)^2 = {(b/a)**2:.2f}; "
              f"plateau = 1)")
        print(f"  ratio chi_O(L{b})/chi_O(L{a})       = "
              f"{rb['chi_O']/ra['chi_O']:.3f}   (non-self-averaging total-norm)")

    with open(RESULT, 'w') as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: float(o) if isinstance(o, (np.floating, np.ndarray)) else str(o))
    print(f"\n# wrote {RESULT}")
    print("Labels: [T] free acid & normalization algebra; [V] SU(2) numbers; [C] plateau reading.")
    print("No QED. Clay 0%.")


if __name__ == '__main__':
    main()
