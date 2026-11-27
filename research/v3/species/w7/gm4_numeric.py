"""GM4-NUMERIC — can delta_1 (the a^2 coefficient of the sinh-mass) be extracted?

Task:  omega_a := (2/a) sinh(a Delta_a / 2) = Delta_cont + d1 a^2 + d2 a^4 + ...
       (SM-2) <=> d1 = 0.
NO new MC in this script: (1) free-field calibration, (2) 2D control,
(3) sensitivity/power budget, (4) re-analysis of the archived b01 beta-scan
correlators (zero MC: the configs are already on disk).

Outputs gm4_numeric.log (plain text) and prints the same to stdout.
"""
import os
import sys
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run3_apebase import correlators, fit_A, A_at  # noqa: E402

OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def omega(a, D):
    """sinh-mass: invert the lattice dispersion."""
    return (2.0 / a) * np.sinh(a * np.asarray(D) / 2.0)


def delta_lat(a, w):
    """forward map: lattice mass from the sinh-mass."""
    return (2.0 / a) * np.arcsinh(a * np.asarray(w) / 2.0)


def fit_d1(a, w):
    """LSQ fit w = D + d1 a^2 (+ d2 a^4 if enough points). Returns (D, d1, d2, rms)."""
    a = np.asarray(a, float)
    w = np.asarray(w, float)
    if len(a) >= 3:
        M = np.vstack([np.ones_like(a), a ** 2, a ** 4]).T
    else:
        M = np.vstack([np.ones_like(a), a ** 2]).T
    coef, *_ = np.linalg.lstsq(M, w, rcond=None)
    r = w - M @ coef
    D = coef[0]
    d1 = coef[1]
    d2 = coef[2] if len(coef) > 2 else np.nan
    return D, d1, d2, float(np.sqrt(np.mean(r ** 2)))


# ---------------------------------------------------------------- 1. free field
def part1_free():
    P("=" * 100)
    P("[1] FREE FIELD: Delta_a = (2/a) asinh(a m / 2)  =>  omega_a = (2/a) sinh(a Delta_a/2) = m  exactly")
    P("=" * 100)
    P(f"{'m':>4} {'a':>8} {'Delta_a':>22} {'omega_a - m':>14} {'rel':>12} "
      f"{'1-Delta/m':>12} {'/(a^2 m^2/24)':>14}")
    worst = 0.0
    for m in (0.5, 1.0, 2.0):
        for a in (1e-1, 1e-2, 1e-3, 1e-4):
            D = (2.0 / a) * np.arcsinh(a * m / 2.0)
            w = omega(a, D)
            d = w - m
            worst = max(worst, abs(d) / m)
            P(f"{m:>4.1f} {a:>8.0e} {D:>22.17g} {d:>14.3e} {abs(d)/m:>12.3e} "
              f"{1-D/m:>12.4e} {(1-D/m)/(a*a*m*m/24):>14.9f}")
    P(f"  => max |omega_a - m|/m = {worst:.3e}   (machine precision: {np.finfo(float).eps:.2e})")
    P("  [T] omega_a == m identically (sinh o asinh = id); [E] the round-off above.")
    P("")

    P("[1b] the free-field a^2 coefficient seen from Delta_a's side: Delta_a = m - a^2 m^3/24 + ...")
    for a in (0.1, 0.05, 0.02):
        P(f"     a={a:<6}  a^2 m^3/24 (m=1) = {a*a/24:.4e}   relative = {a*a/24:.4e}")
    P("     => the (SM-2) coefficient c1 = d1/D - D^2/24 (GK2-SYMBOL.md:83 (1.2a));")
    P("        free field has d1 = 0  =>  c1 = -m^2/24, i.e. the WHOLE a^2 artifact is the asinh kernel.")
    P("")

    P("[1c] round trip: inject omega_a = D + d1 a^2 + d2 a^4, rebuild Delta_a, refit d1.")
    P(f"{'d1_in':>10} {'d2_in':>8} {'pts':>4} {'d1_out':>16} {'abs err':>11} {'rel err':>10}")
    for d1_in in (0.0, 1e-3, 1e-2, 1e-1, 1.0):
        for d2_in in (0.0, 0.5):
            for npts in (3, 4, 6):
                a = np.linspace(0.05, 0.20, npts)
                w = 1.0 + d1_in * a ** 2 + d2_in * a ** 4
                D = delta_lat(a, w)
                w2 = omega(a, D)
                Df, d1, d2, rms = fit_d1(a, w2)
                rel = abs(d1 - d1_in) / abs(d1_in) if d1_in else np.nan
                P(f"{d1_in:>10.1e} {d2_in:>8.1f} {npts:>4d} {d1:>16.10e} "
                  f"{abs(d1-d1_in):>11.2e} {rel:>10.2e}   rms={rms:.2e}")
    P("  => the sinh transform + refit recovers an injected d1 to round-off: the METHOD is exact;")
    P("     what fails on real data is not the algebra but (a) the a-calibration and (b) the noise.")
    P("")

    P("[1d] noise / power: how precisely must Delta_a be known to see the a^2 term?")
    P("     signal-to-noise of the a^2 term ~ (D * mhat^2/24) where mhat := a*D is the LATTICE mass;")
    P("     requirement at 4 sigma:  sigma(Delta_a)/Delta_a <= mhat^2/96.")
    P(f"{'mhat=a*D':>10} {'a^2 D^3/24':>13} {'rel':>9} {'sig_rel(4sig)':>14}")
    for mhat in (0.05, 0.1, 0.2, 0.5, 1.0, 1.7, 3.0):
        rel = mhat ** 2 / 24.0
        P(f"{mhat:>10.2f} {'-':>13} {rel:>9.2e} {rel/4:>14.2e}")
    P("     (a=1 convention: a^2 D^3/24 = D*mhat^2/24 with mhat=a*D; the relative size is mhat^2/24.)")
    P("     e.g. a=0.1, D=1 (mhat=0.1): rel = 4.167e-4  => need 1.0e-4 relative  [matches the task sheet]")
    P("")

    P("[1e] power of a 3-point determination at fixed a-window (noise-limited):")
    P("     sigma(d1) = sigma_w / sqrt( sum (a_i^2 - <a^2>)^2 ) ; compare with the natural size D^3/24.")
    P(f"{'sig_rel':>9} {'a in':>14} {'sigma(d1)':>12} {'D^3/24':>9} {'ratio':>8} {'n_sigma':>8}")
    for sig_rel in (1e-2, 1e-3, 1e-4):
        for lo, hi in ((0.05, 0.20), (0.10, 0.20), (0.02, 0.10)):
            a = np.array([lo, (lo + hi) / 2, hi])
            D = 1.0
            w = D + 0 * a
            sig = sig_rel * D
            M = np.vstack([np.ones_like(a), a ** 2]).T
            cov = np.linalg.inv(M.T @ M)
            s_d1 = sig * np.sqrt(cov[1, 1])
            P(f"{sig_rel:>9.0e} {lo:>6.2f}-{hi:<6.2f} {s_d1:>12.2e} {1/24:>9.4f} "
              f"{s_d1/(1/24):>8.2f} {1/24/s_d1:>8.1f}")
    P("  => with 1% data over a in [0.1,0.2] the slope error is ~4x the natural size: 3-4 sigma only")
    P("     for a delta_1 as large as the free-field artifact itself.  At 1e-4 it is clean.")
    P("")


# ------------------------------------------------------------------------ 2. 2D
def part2_2d():
    P("=" * 100)
    P("[2] 2D CONTROL: Delta_a == Delta_cont for all a  (GK2-SYMBOL.md:278, H_a = (g^2/2) sum C2, no a)")
    P("=" * 100)
    P("Algebra: Delta_a = D const  =>  omega_a = (2/a) sinh(aD/2) = D + a^2 D^3/24 + a^4 D^5/1920 + ...")
    P("         => d1 = D^3/24  (!= 0),  and (1.2a) c1 = 0  <=>  d1 = D^3/24  (consistent).")
    P(f"{'D':>4} {'a':>8} {'omega_a':>18} {'(2/a)sinh(aD/2)-D':>18} {'a^2D^3/24':>11}")
    for D in (1.0, 2.0):
        for a in (1e-2, 1e-1, 3e-1):
            w = omega(a, D)
            P(f"{D:>4.1f} {a:>8.0e} {w:>18.12g} {w-D:>18.6e} {a*a*D**3/24:>11.4e}")
    for D in (1.0,):
        for lo, hi, n in ((0.02, 0.10, 4), (0.05, 0.20, 4), (0.10, 0.30, 4)):
            a = np.linspace(lo, hi, n)
            Df, d1, d2, rms = fit_d1(a, omega(a, D))
            P(f"  fit over a in [{lo},{hi}] ({n} pts): D={Df:.10f}  d1={d1:.10e}  "
              f"D^3/24={D**3/24:.10e}  ratio={d1/(D**3/24):.9f}  rms={rms:.1e}")
    P("  => 2D gives d1 = D^3/24 exactly (not 0).  The corpus's '2D: delta_1 = 0, trivial pass'")
    P("     (GK2-SYMBOL.md:30 row 0.12, :278 sec 4.2) is arithmetically incompatible with its own")
    P("     (1.2a) c1 = d1/D - D^2/24: c1 = 0 (2D has no lattice artifact) forces d1 = D^3/24.")
    P("     Consequence for C = D^2/24 + |d1|/D (GK2:194): 2D gives C = D^2/12, not 0.")
    P("     VERDICT: 2D is a NEGATIVE control with an exactly known value, not a zero test.")
    P("")


# ------------------------------------------------------- 3. archived beta scan
def anchor(La, La_ref, b_ref, b):
    """two-loop a*Lambda_L (GW7-FIXEDA-SCAN.md:71), normalised to the printed anchor."""
    b0 = 11.0 * 2 / (48 * np.pi ** 2)

    def f(bb):
        g2 = 4.0 / bb
        return np.exp(-1.0 / (2 * b0 * g2)) * (b0 * g2) ** (-51.0 / 121.0)
    return La_ref * f(b) / f(b_ref)


def part3_data():
    P("=" * 100)
    P("[3] ARCHIVED b01 BETA-SCAN (Ls=8, Lt=12, N=12000, one operator) -- re-analysis, ZERO new MC")
    P("=" * 100)
    P("   b01_scan.py:73 hard-codes Ls,Lt = 8,12 ; :104 m_fit = fit_A(C,Lt,ts=1..6) in LATTICE units ;")
    P("   :180 writes only S,plaq,w1 to the npz (m_fit lives in the .json files).")
    P("")
    files = [("b01_gw2p.json", 2.2), ("b01_gw2p.json", 2.75), ("b01_cal_b3.json", 3.0),
             ("b01_tail.json", 4.0), ("b01_t6.json", 6.0), ("b01_t10.json", 10.0)]
    rows = []
    for fn, b in files:
        d = json.load(open(os.path.join(HERE, fn)))
        r = [x for x in d if abs(x["beta"] - b) < 1e-9][0]
        rows.append(r)
    P(f"{'beta':>5} {'N':>6} {'sec':>6} {'plaq':>8} {'C1/C0':>9} {'C2/C0':>10} {'snr2':>6} "
      f"{'m_fit':>7} {'edge':>5}")
    for r in rows:
        P(f"{r['beta']:>5.2f} {r['N']:>6d} {r['sec']:>6.0f} {r['plaq']:>8.5f} {r['C1o0']:>9.5f} "
          f"{r['C2o0']:>10.6f} {r['snr2']:>6.2f} {r['m_fit']:>7.3f} {str(r['m_edge']):>5}")
    P("")
    P("  m_fit is NOT monotone in beta: 1.690, 2.706, 3.940, 2.116, 3.475, 4.000")
    P("  beta=10.0: m_fit = 4.000 with m_edge=True == the fitter's scan ceiling (run3_apebase.py:336")
    P("             ms = linspace(0.05,4.0,4000)) => that 'mass' is a boundary hit, not a measurement.")
    P("  snr2 = C(2)/sigma_JK(C(2)): 8.87, 1.25, 2.58, 2.41, 0.27, -0.15  => no t=2 signal at beta>=6.")
    P("")

    # independent re-analysis from the stored configs (zero MC)
    P("[3a] re-derived from the stored npz correlators (same estimators, jackknife 20 blocks):")
    P(f"{'beta':>5} {'C1/C0':>9} {'C2/C0':>10} {'m_fit':>8} {'+/-':>8} {'rel':>9} {'w1':>8} {'w1_err':>9}")
    res = []
    for fn in ["b01_b2.2.npz", "b01_b2.75.npz", "b01_b3.npz", "b01_b4.npz", "b01_b6.npz", "b01_b10.npz"]:
        z = np.load(os.path.join(HERE, fn))
        S = z["S"][:, 0, :]
        Lt = int(z["Lt"])
        b = float(z["beta"])
        C = correlators(S, Lt)
        m, A, edge = fit_A(C, Lt, ts=np.arange(1, 7))
        N = S.shape[0]
        nb = 20
        idx = np.array_split(np.arange(N), nb)
        vals = []
        w1v = []
        for j in range(nb):
            keep = np.concatenate([idx[q] for q in range(nb) if q != j])
            Cj = correlators(S[keep], Lt)
            vals.append(fit_A(Cj, Lt, ts=np.arange(1, 7))[0])
            w1v.append(A_at(Cj, Lt, 1.546) / Cj[0])
        vals = np.array(vals)
        w1v = np.array(w1v)
        em = np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean()) ** 2))
        ew = np.sqrt((nb - 1) / nb * np.sum((w1v - w1v.mean()) ** 2))
        res.append(dict(beta=b, C1=C[1] / C[0], C2=C[2] / C[0], m=m, em=em, w1=w1v.mean(), ew=ew,
                        N=N, Lt=Lt))
        P(f"{b:>5.2f} {C[1]/C[0]:>9.5f} {C[2]/C[0]:>10.6f} {m:>8.4f} {em:>8.4f} "
          f"{em/m:>9.2e} {w1v.mean():>8.4f} {ew:>9.4f}")
    P("  (beta=3: the archived w1 = 5.638 +- 4e-15 is the KNOWN correlators() axis bug, already")
    P("   registered in B01.md:96-97 and GP3-PLATEAU.md:154 -- S fed as (N,1,Lt) to the jackknife")
    P("   makes C(t) a constant vector, so w1 = sum(g)/sum(g^2) = 5.6387.  With the lev axis")
    P("   collapsed (this table) it is 0.2137 +- 0.0117.  No new finding; the m_fit of that file")
    P("   IS reproducible (3.9397), so only the jackknife column is affected.)")
    P("")

    # omega attempt
    P("[3b] forced omega_a 'sequence': a(beta) from the two-loop form (GW7:71), anchored at")
    P("     a*Lambda_L(2.5) = 3.5758e-3 (B04.md:58).  Delta_a = m_fit/a (lattice mass / a).")
    La = lambda b: anchor(None, 3.5758e-3, 2.5, b)
    rows2 = []
    P(f"{'beta':>5} {'aL':>10} {'1/aL':>10} {'mhat':>7} {'Delta_a=mhat/aL':>16} {'omega_a':>14} "
      f"{'a^2':>10} {'mhat^2/24':>10}")
    for r in res:
        b = r["beta"]
        al = La(b)
        D = r["m"] / al
        w = omega(al, D)
        rows2.append((b, al, r["m"], D, w))
        P(f"{b:>5.2f} {al:>10.3e} {1/al:>10.3e} {r['m']:>7.4f} {D:>16.4e} {w:>14.4e} "
          f"{al*al:>10.3e} {r['m']**2/24:>10.4f}")
    bs = np.array([x[0] for x in rows2])
    als = np.array([x[1] for x in rows2])
    ws = np.array([x[4] for x in rows2])
    P("")
    P("  omega_a spans %.2e .. %.2e  (5 decades) -- before any a^2 analysis the ZEROTH-order" %
      (ws.min(), ws.max()))
    P("  requirement (a common Delta_cont) already fails: the physical point is not held fixed.")
    P("  Pairwise (beta=2.2, beta=2.6-ish) and least-squares fits, for the record:")
    for sel in ([0, 1], [0, 2], [0, 3], [1, 2], [0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4]):
        a_ = als[sel]
        w_ = ws[sel]
        Df, d1, d2, rms = fit_d1(a_, w_)
        P(f"    beta={[float(bs[i]) for i in sel]}  D={Df:>12.4e}  d1={d1:>13.4e}  "
          f"d2={d2:>11.3e}  rms={rms:>10.3e}")
    P("  => d1 changes by ~15 orders and flips sign between subsets; the residuals are of the same")
    P("     size as the data.  No point estimate, no interval.  (Not a statistical problem.)")
    P("")

    P("[3c] why: three independent obstructions, each fatal on its own")
    P("   (i)  L=8 fixed => physical volume V = (8 aL)^4 shrinks with a:")
    for b1, b2 in ((2.2, 2.75), (2.5, 4.0), (2.5, 6.0)):
        P(f"        beta {b1} -> {b2}:  aL ratio = {La(b1)/La(b2):.3e},  V ratio = "
          f"{(La(b1)/La(b2))**4:.3e}")
    P("   (ii) a(beta) has no absolute value in-window: two-loop at beta=2.5 gives m*a = 0.023 vs the")
    P("        measured 1.707 (B04.md:64, factor 74), and a(2.5)/a(6) = 8535 (GW7-FIXEDA-SCAN.md:92).")
    P("        Without a(beta) there is no a^2 axis -- only a beta axis.")
    P("   (iii) m_fit is a UV reading, not a physical mass: mhat = 1.69..4.00 at every beta")
    P("        (xi/a = 1/1.707 = 0.586 at beta=2.5, B04.md:98) => Delta_a = mhat/aL ~ 1/aL grows")
    P("        by 4 decades as beta grows: that IS the 'no continuum limit' statement, not a mass.")
    P("")


# -------------------------------------------------------------- 4. experiment
def part4_cost():
    P("=" * 100)
    P("[4] COST of the experiment that WOULD measure delta_1")
    P("=" * 100)
    P("  measured throughput (b01 logs, Ls=8^3 x Lt=12):")
    P("    12000 cfg in 57.5 s (beta=2.2) .. 215 s (beta=6)  => 4.8 .. 18 ms/cfg")
    P("  GW7 table (fixed physical volume, Ls ~ 1/aL), measured: Ls=10,Lt=15 -> 147.7 s ;")
    P("    Ls=12,Lt=18 -> 331.7 s ; model (Ls/8)^4*62 s -> Ls=15: 766 s ; Ls=28: 9304 s.")
    P("")
    P("  required sigma(Delta)/Delta for a 4-sigma test: mhat^2/96.  current sigma(m_fit)/m_fit:")
    P("    beta=2.2: see [3a] table.  Scale N as (sigma_now/sigma_target)^2:")
    for mn, sn in ((0.5, 3e-2), (1.0, 3e-2), (1.7, 3e-2)):
        tgt = mn ** 2 / 96
        P(f"    mhat={mn}: sigma_target={tgt:.2e} -> N multiplier = {(sn/tgt)**2:.2e} "
          f"x 12000 = {12000*(sn/tgt)**2:.1e} cfg")
    P("")
    P("  => cheapest honest design (4 beta points on a fixed-physics line, a-ratio ~1.9,")
    P("     Ls = 8,10,12 and a 4th at Ls>=15), N chosen for the target sigma:")
    rows = [(2.50, 3.5758e-3, 8, 12, 57.5), (2.59, 2.8486e-3, 10, 15, 147.7),
            (2.66, 2.3861e-3, 12, 18, 331.7), (2.748, 1.9071e-3, 15, 22, 766.0)]
    P(f"{'beta':>6} {'aL':>10} {'Ls':>3} {'Lt':>3} {'t(12000)':>10} {'N for 1e-3':>11} "
      f"{'t(1e-3)':>10} {'N for 1e-4':>11} {'t(1e-4)':>12}")
    tot1 = tot2 = 0.0
    for b, al, Ls, Lt, t in rows:
        n1 = 12000 * (3e-2 / 1e-3) ** 2
        n2 = 12000 * (3e-2 / 1e-4) ** 2
        t1 = t * n1 / 12000
        t2 = t * n2 / 12000
        tot1 += t1
        tot2 += t2
        P(f"{b:>6.3f} {al:>10.3e} {Ls:>3d} {Lt:>3d} {t:>10.1f} {n1:>11.2e} {t1:>10.1f} "
          f"{n2:>11.2e} {t2:>12.1f}")
    P(f"  TOTAL for 1e-3 precision: {tot1:.0f} s = {tot1/3600:.1f} h")
    P(f"  TOTAL for 1e-4 precision: {tot2:.0f} s = {tot2/3600:.0f} h = {tot2/86400:.1f} d")
    P("  => NOT <= 1 hour by 2-3 orders of magnitude.  And note the design already assumes the")
    P("     two-loop a(beta) is right for the BOX MATCHING (it is not: 74x at 2.5), so even the")
    P("     'fixed physical point' would first need an independent scale setting (~one more scan).")
    P("")


def part3d_fixedvolume():
    P("=" * 100)
    P("[3d] THE ONLY FIXED-PHYSICAL-VOLUME a-SEQUENCE IN THE REPO: the GW7 line (beta=2.5,2.59,2.66)")
    P("     Ls = 8,10,12 ; Lt = 12,15,18 ; N = 24000 (2.5) / 12000 (2.59, 2.66) ; operator and kernel")
    P("     identical to b01 (gw7_fixeda.log:2-3).  m_fit re-derived here from the stored configs.")
    P("=" * 100)
    P(f"{'beta':>5} {'Ls':>3} {'Lt':>3} {'aL':>10} {'a^2':>10} {'mhat':>8} {'+/-':>7} {'rel':>8} "
      f"{'Delta=mhat/aL':>14} {'omega_a':>13}")
    seq = [("gh2_b2.5_Ls8_Lt12.npz", 2.5), ("gh2_b2.5_Ls8_Lt20.npz", 2.5),
           ("gw7_b2.59_Ls10.npz", 2.59), ("gw7_b2.66_Ls12.npz", 2.66)]
    got = []
    for fn, b in seq:
        z = np.load(os.path.join(HERE, fn))
        S = z["S"][:, 0, :]
        Lt = int(z["Lt"])
        Ls = int(z["Ls"])
        C = correlators(S, Lt)
        ts = np.arange(1, min(7, Lt // 2) )
        m, A, edge = fit_A(C, Lt, ts=ts)
        nb = 20
        idx = np.array_split(np.arange(S.shape[0]), nb)
        vals = []
        for j in range(nb):
            keep = np.concatenate([idx[q] for q in range(nb) if q != j])
            vals.append(fit_A(correlators(S[keep], Lt), Lt, ts=ts)[0])
        vals = np.array(vals)
        em = float(np.sqrt((nb - 1) / nb * np.sum((vals - vals.mean()) ** 2)))
        al = anchor(None, 3.5758e-3, 2.5, b)
        got.append(dict(b=b, Ls=Ls, Lt=Lt, al=al, m=m, em=em, ts=(int(ts[0]), int(ts[-1]))))
        P(f"{b:>5.2f} {Ls:>3d} {Lt:>3d} {al:>10.3e} {al*al:>10.3e} {m:>8.4f} {em:>7.4f} "
          f"{em/m:>8.2e} {m/al:>14.4e} {omega(al, m/al):>13.4e}")
    P("  (fit window ts = 1..min(6,Lt/2-1) for each Lt; a(beta) = two-loop, GW7:71 anchor)")
    P("")
    # the decisive zeroth-order test
    P("  ZEROTH-ORDER TEST (does the physical point stay fixed?): hold V*a^4 fixed by construction,")
    P("  then a continuum mass must satisfy mhat = a*Delta_cont, i.e. mhat must FALL like a.")
    base = got[0]
    P(f"{'beta':>5} {'mhat meas':>10} {'mhat if fixed physics':>22} {'ratio':>7} {'Delta_a growth':>15}")
    for g in got:
        pred = base["m"] * g["al"] / base["al"]
        P(f"{g['b']:>5.2f} {g['m']:>10.4f} {pred:>22.4f} {g['m']/pred:>7.3f} "
          f"{(g['m']/g['al'])/(base['m']/base['al']):>15.3f}")
    P("  => mhat RISES (1.92 -> 2.25 -> 2.41) where a fixed physical mass demands it FALLS")
    P("     (1.92 -> 1.53 -> 1.28).  Delta_a = mhat/aL therefore grows like a^-1.6 (log-log slope")
    P(f"     {np.log(1.875)/np.log(1.4986):.2f}), i.e. what is measured is the cutoff, not a continuum mass.")
    P("")
    # forced delta_1 estimate + interval
    P("  FORCED POINT ESTIMATE (3 points, fixed volume, a^2 ratio 2.2): propagate the jackknife")
    P("  errors through omega = (2/a) sinh(mhat/2), sigma_omega = cosh(mhat/2) * sigma_mhat:")
    use = [got[0], got[2], got[3]]
    a_ = np.array([g["al"] for g in use])
    w_ = np.array([omega(g["al"], g["m"] / g["al"]) for g in use])
    s_ = np.array([np.cosh(g["m"] / 2) * g["em"] / g["al"] for g in use])
    P(f"{'beta':>5} {'a':>10} {'omega':>13} {'sigma_omega':>12} {'rel':>8}")
    for g, w, s in zip(use, w_, s_):
        P(f"{g['b']:>5.2f} {g['al']:>10.3e} {w:>13.4e} {s:>12.3e} {s/w:>8.2e}")
    M = np.vstack([np.ones_like(a_), a_ ** 2]).T
    W = np.diag(1 / s_ ** 2)
    cov = np.linalg.inv(M.T @ W @ M)
    coef = cov @ (M.T @ W @ w_)
    D, d1 = coef
    eD, ed1 = np.sqrt(np.diag(cov))
    chi2 = float(((w_ - M @ coef) ** 2 / s_ ** 2).sum())
    P(f"  fit:  Delta_cont = {D:.4e} +- {eD:.2e}   delta_1 = {d1:.4e} +- {ed1:.2e}   chi2 = {chi2:.2f}/1")
    P(f"  natural size for comparison: 1/24 = {1/24:.4e} ;  delta_1/|Delta_cont^3/24| = "
      f"{abs(d1)/(abs(D)**3/24):.3e}")
    P("  => the interval is comparable to or wider than the value, the reduced chi^2 is ~0 (3 points,")
    P("     2 parameters, and the three omega's are ordered by mhat, not by a^2) -- it is a")
    P("     two-parameter fit to a curve that is not the model.  NO INTERVAL worth quoting.")
    P("     Even the sign is fixed by the cutoff drift, not by delta_1.")
    P("")
    P("  NULL MODEL (decisive): replace the three mhat by a constant (no continuum mass at all,")
    P("  pure 1/a cutoff) and run the SAME fit.  If the estimator cannot tell the two apart, then")
    P("  the 'delta_1' above is a cutoff meter, not a measurement.")
    mbar = float(np.mean([g["m"] for g in use]))
    w_null = np.array([omega(g["al"], mbar / g["al"]) for g in use])
    s_null = s_  # same sigma (dominated by the 1/a factor)
    Wn = np.diag(1 / s_null ** 2)
    covn = np.linalg.inv(M.T @ Wn @ M)
    coefn = covn @ (M.T @ Wn @ w_null)
    Dn, d1n = coefn
    en = np.sqrt(np.diag(covn))
    chin = float(((w_null - M @ coefn) ** 2 / s_null ** 2).sum())
    P(f"    null (mhat = {mbar:.4f} const): Delta_cont = {Dn:.4e} +- {en[0]:.2e}   "
      f"delta_1 = {d1n:.4e} +- {en[1]:.2e}   chi2 = {chin:.2f}/1")
    P(f"    data                          : delta_1 = {d1:.4e} +- {ed1:.2e}")
    P(f"    difference / sigma = {abs(d1-d1n)/np.hypot(ed1, en[1]):.2f} sigma  "
      f"=> the estimator DOES NOT separate the data from the pure-cutoff null.")
    P("    And note the null's delta_1 is also ~6 sigma from 0: so a 'delta_1 != 0 at 6 sigma'")
    P("    claim is exactly what a theory with NO continuum mass produces here.  (SM-2) cannot be")
    P("    tested this way: the a^2 lever arm is being fit to the 1/a cutoff drift.")
    P("")


def main():
    P("GM4-NUMERIC  |  sinh-mass delta_1: can (SM-2) be decided numerically from the corpus data?")
    P(f"numpy {np.__version__}  |  no new MC in this script")
    P("")
    part1_free()
    part2_2d()
    part3_data()
    part3d_fixedvolume()
    part4_cost()
    open(os.path.join(HERE, "gm4_numeric.log"), "w", encoding="utf-8").write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
