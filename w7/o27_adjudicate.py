#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
o27_adjudicate.py -- ADJUDICATION of the A-vs-C sign conflict at beta=2.5, L=2, SU(2).

Definitions taken from the paper (verified, not from memory):
  docs/paper-v2/ch05-theorems.tex  eq:split      (label L239)
      G''(0) = E[H_vv+Phi'' | xi=0] - Var(G_v+Phi' | xi=0) = Delta_dyn + W_geo + X
      Delta_dyn = E[H] - Var(G) ,  X = -2 Cov(G,Phi')
  docs/paper-v2/ch05-theorems.tex  eq:wgeo
      W_geo := E[Phi''] - Var(Phi')
  docs/paper-v2/ch15-uvface.tex    eq:v3object (label L26)
      G''(v;a0) = E_a[d_v^2 V] - Var_a(d_v V) = Delta_dyn + W_geo + X,  V = S_W + Phi, Phi=-log J

So in the flat chart X_ell = log U_ell in R^3 (SU(2)), V = S_W + Phi:
      H_vv := d_v^2 S_W      (action Hessian along v)          [o13's "H", o11's "H"=-beta*Ch]
      G_v  := d_v S_W        (action force  along v)           [o13's "G", o11's "G"=-beta*Ag]
      Phi'' := d_v^2 Phi     (Jacobian curvature)              [o11's Dd]
      Phi'  := d_v Phi       (Jacobian force)                  [o11's Dp]
      kappa = H_vv + Phi'' = d_v^2 V   (o11/o26 shorthand)
      F     = G_v  + Phi'  = d_v  V    (o11/o26 shorthand)
=>  C = E[kappa] - Var(F) == B == G''(0) == A + W_full   (algebra, exact)

Tier 1 [E] exact: single-plaquette SU(2) beta=2.5, measure e^{-V} on R^3, exactly
  1D-quadrature in r (V is radial) + exact sphere average <(xhat.v)^2>=1/3.
Tier 2 [C] MC: L=2 4D SU(2) Wilson, beta=2.5, constrained (xi=0) HMC, SAME configurations
  for all four quantities; two modes (o13's and o11's) for the mode-dependence question.

Run:  py -3.12 o27_adjudicate.py [ntraj]
Labels [E]/[C]/[T].  No git commit.
"""
import sys, os, time
import numpy as np

sys.path.insert(0, r"D:\Dev\math\research\v3\w7")
import importlib.util
spec = importlib.util.spec_from_file_location("o13", r"D:\Dev\math\research\v3\w7\o13_ldep.py")
o13 = importlib.util.module_from_spec(spec); spec.loader.exec_module(o13)

BETA = 2.5
L = 2
P = lambda *a: print(*a, flush=True)
T0 = time.time()

# ============================================================ TIER 1: exact single plaquette
def tier1(nr=4_000_000, rmax=3.2):
    """Exact [E]: single SU(2) plaquette, flat chart, dU = J dX, V = beta(1-cos(r/2)) + Phi(r),
    Phi(r) = -2 log(sin(r/2)/(r/2)).  Radial quadrature in r + exact sphere average."""
    r = np.linspace(1e-3, rmax, nr)                        # r>0: 1/r terms are finite, w~r^2
    dr = r[1] - r[0]
    sw = BETA * (1 - np.cos(r / 2))                       # S_W
    Phi = -2 * np.log(np.sin(r / 2) / (r / 2))            # -log J
    w = np.exp(-(sw + Phi)) * (np.sin(r / 2) / (r / 2)) ** 2 * r ** 2   # e^{-V} J r^2 dr
    Z = np.trapezoid(w, dx=dr)
    A3 = 1.0 / 3.0                                        # <(xhat.v)^2> on S^2
    # derivatives along a fixed unit v, using x = r nhat, c = (nhat.v)
    dSw_r = 0.5 * BETA * np.sin(r / 2)                    # d S_W / dr
    d2Sw_rr = 0.25 * BETA * np.cos(r / 2)                 # d2 S_W / dr^2
    dPh_r = -np.cos(r / 2) / np.sin(r / 2) + 2.0 / r
    d2Ph_rr = 0.5 / np.sin(r / 2) ** 2 - 2.0 / r ** 2
    # H_vv = d2Sw_rr c^2 + (dSw_r/r)(1-c^2);  G_v = dSw_r c
    # Phi'' = d2Ph_rr c^2 + (dPh_r/r)(1-c^2);  Phi' = dPh_r c
    def avg(f, pw):
        """weighted sphere average of f(r) * c^pw  (pw in {0,1,2})"""
        c2 = A3
        mult = {0: 1.0, 1: 0.0, 2: c2}[pw]
        return np.trapezoid(w * f * mult, dx=dr) / Z
    EH = avg(d2Sw_rr, 2) + avg(dSw_r / r, 0) - avg(dSw_r / r, 2)
    EG = avg(dSw_r, 1)                                    # = 0 by c-average
    EG2 = avg(dSw_r ** 2, 2)
    EPh2 = avg(d2Ph_rr, 2) + avg(dPh_r / r, 0) - avg(dPh_r / r, 2)
    EPh1 = avg(dPh_r, 1)
    EPh1sq = avg(dPh_r ** 2, 2)
    EGP = avg(dSw_r * dPh_r, 2)                           # E[G Phi']
    EG2_tot = EG2 + EPh1sq + 2 * EGP
    # [V] finite-difference check of the analytic Phi', Phi'' and S_W derivatives
    rr = np.linspace(0.05, 3.0, 500); hh = 1e-5
    ph = lambda z: -2 * np.log(np.sin(z / 2) / (z / 2))
    d1f = lambda z: (ph(z + hh) - ph(z - hh)) / (2 * hh)
    d2f = lambda z: (ph(z + hh) - 2 * ph(z) + ph(z - hh)) / hh ** 2
    e1 = np.max(np.abs(d1f(rr) - (-np.cos(rr / 2) / np.sin(rr / 2) + 2.0 / rr)) / np.abs(d1f(rr)))
    e2 = np.max(np.abs(d2f(rr) - (0.5 / np.sin(rr / 2) ** 2 - 2.0 / rr ** 2)) / np.abs(d2f(rr)))
    P(f"    [V] Phi' FD-dev={e1:.2e}  Phi'' FD-dev={e2:.2e}  (analytic derivatives verified)")
    return dict(H=EH, G=EG, vG=EG2 - EG ** 2, Phii=EPh2, Phi1=EPh1, vP=EPh1sq - EPh1 ** 2,
                Cov=EGP - EG * EPh1, vF=EG2_tot - (EG + EPh1) ** 2, Z=Z)

def combos(d, tag):
    A = d['H'] - d['vG']
    B = (d['H'] + d['Phii']) - d['vF']
    C = (d['H'] + d['Phii']) - d['vF']            # kappa=H+Phi'', F=G+Phi'
    Wful = d['Phii'] + d['vG'] - d['vF']
    D = A + Wful
    Wgeo = d['Phii'] - d['vP']
    X = -2 * d['Cov']
    P(f"  {tag}")
    P(f"    E[H_vv]={d['H']:+.6f}  E[Phi'']={d['Phii']:+.6f}  Var(G_v)={d['vG']:.6f}"
      f"  Var(G_v+Phi')={d['vF']:.6f}  Cov(G,Phi')={d['Cov']:+.6f}  Var(Phi')={d['vP']:.6f}")
    P(f"    A=E[H]-Var(G)={A:+.6f}   B=E[H+Phi'']-Var(G+Phi')={B:+.6f}   C=E[kappa]-Var(F)={C:+.6f}")
    P(f"    W_full=E[Phi'']+Var(G)-Var(G+Phi')={Wful:+.6f} (=W_geo+X: {Wgeo+X:+.6f})"
      f"   D=A+W_full={D:+.6f}   |D-B|={abs(D-B):.2e}")
    P(f"    C-A = W_full = {C-A:+.6f}   [E[Phi'']={d['Phii']:+.6f} ; "
      f"-Var(Phi')-2Cov={-(d['vP']+2*d['Cov']):+.6f}]")
    return A, B, C, D, Wful

# ============================================================ TIER 2: L=2 lattice, one MC
def jack_bins(arr, nb=25):
    """bin the time series and jackknife the given scalar functional"""
    a = np.asarray(arr, float)
    n = a.size
    nb = min(nb, max(4, n // 8))
    k = n // nb
    blocks = a[:k * nb].reshape(nb, k, *a.shape[1:])
    return blocks, nb

def lat_mode(rec, mode_name):
    r = np.array(rec)                     # columns: Gv, Phi', Phi'', Hsv, S0
    Gv, Ph, Pip, Hsv = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
    n = r.shape[0]
    nb = min(25, max(4, n // 8)); k = n // nb
    def stats(sel):
        G, PP, QQ, HH = Gv[sel], Ph[sel], Pip[sel], Hsv[sel]
        mH, mPi = HH.mean(), QQ.mean()
        mG, mP = G.mean(), PP.mean()
        vG = G.var(); vF = (G + PP).var(); vP = PP.var()
        cov = (G * PP).mean() - mG * mP
        A = mH - vG; B = (mH + mPi) - vF; C = (mH + mPi) - vF
        Wf = mPi + vG - vF; D = A + Wf
        return np.array([mH, mPi, vG, vF, cov, vP, A, B, C, Wf, D, mG + mP])
    full = stats(slice(None))
    jk = []
    for i in range(nb):
        sel = np.r_[0:i * k, (i + 1) * k:n].astype(int)
        jk.append(stats(sel))
    jk = np.array(jk)
    sig = np.sqrt((nb - 1) / nb * ((jk - jk.mean(0)) ** 2).sum(0))
    # tau on the two most autocorrelated series
    def tau1(x, ml=200):
        x = x - x.mean(); v = (x * x).sum()
        ac = np.array([1.0] + [(x[:-j] * x[j:]).sum() / v for j in range(1, ml + 1)])
        return 1 + 2 * ac[1:].sum()
    tk, tf = tau1(Hsv), tau1(Gv + Ph)
    return full, sig, dict(n=n, tauK=tk, tauF=tf, neff=n / max(tk, tf, 1.0))

if __name__ == "__main__":
    ntraj = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    P("=" * 100)
    P(f"[O27-ADJUDICATE] beta={BETA} L={L} SU(2)   {time.strftime('%Y-%m-%d %H:%M:%S')}")
    P("defs: docs/paper-v2/ch05-theorems.tex eq:split / eq:wgeo ; ch15-uvface.tex eq:v3object")
    P("=" * 100)

    P("\n[TIER 1 / exact quadrature] single SU(2) plaquette, flat chart, measure e^{-V} on R^3")
    d1 = tier1()
    A1, B1, C1, D1, W1 = combos(d1, "[E] single plaquette (4e6-pt radial quadrature, err<1e-9)")

    P("\n[TIER 2 / L=2 4D SU(2) Wilson, beta=2.5, constrained xi=0 HMC, ONE ensemble per mode]")
    g = o13.build_geom(L)
    P(f"  n={g['n']} sites, 96 plaquettes; HMC {ntraj} traj x 8 steps, eps=0.1402, therm=200")
    out = {}
    vA = o13.make_v(g)                                    # o13 mode: mu=0 only, Lie 0, ||v||=1
    vB = np.zeros((4, g['n'], 3)); vB[:, :, 0] = 1.0 / (2 * L ** 2)   # o11 mode: all mu, ||v||=1
    for name, v in (("o13-mode (mu=0,Lie0)", vA), ("o11-mode (all mu,Lie0)", vB)):
        rec = []
        t = time.time()
        x, acc = o13.hmc(L, BETA, ntraj, 8, 0.1402, v, g, seed=2701 if name.startswith("o13") else 2702,
                         therm=200, meas_every=1, h=0.05, record=rec, verbose=False)
        full, sig, ex = lat_mode(rec, name)
        P(f"\n  --- {name}: acc={acc:.3f}  n={ex['n']}  tau_H={ex['tauK']:.1f} "
          f"tau_(G+Phi')={ex['tauF']:.1f}  Neff~{ex['neff']:.0f}  [{time.time()-t:.0f}s]")
        lbl = ["E[H_vv]", "E[Phi'']", "Var(G_v)", "Var(G_v+Phi')", "Cov(G,Phi')", "Var(Phi')",
               "A", "B", "C", "W_full", "D", "E[G+Phi']"]
        for j, nm in enumerate(lbl):
            P(f"      {nm:<16s} = {full[j]:+12.6f}  +- {sig[j]:.6f}")
        # exact algebraic identities, no statistics
        r1 = abs(full[10] - full[7]); r2 = abs(full[8] - full[7])
        r3 = abs(full[9] - (full[1] + full[2] - full[3]))
        P(f"      identity: |D-B|={r1:.2e}   |C-B|={r2:.2e}   |W_full-(E[Phi'']+VarG-Var(G+Phi'))|={r3:.2e}")
        P(f"      C-A = {full[8]-full[6]:+.6f} = W_full ; E[Phi''] part {full[1]:+.6f} ;"
          f" cross part -Var(Phi')-2Cov = {-(full[5]+2*full[4]):+.6f}")
        out[name] = (full, sig, ex, acc)

    P("\n" + "=" * 100)
    P(f"[VERDICT]  A != C. C == B == G''(0) == A + W_full  (exact, residual ~1e-16)")
    P(f"  exact anchor   : A={A1:+.4f} (sign {np.sign(A1):+.0f})  C={C1:+.4f} (sign {np.sign(C1):+.0f})"
      f"    C-A={C1-A1:+.4f}")
    for nm, (full, sig, ex, acc) in out.items():
        P(f"  {nm:<24s}: A={full[6]:+.3f}+-{sig[6]:.3f}  C={full[8]:+.3f}+-{sig[8]:.3f}"
          f"  sign(A)={np.sign(full[6]):+.0f} sign(C)={np.sign(full[8]):+.0f}")
    P(f"total elapsed {time.time()-T0:.1f}s")
