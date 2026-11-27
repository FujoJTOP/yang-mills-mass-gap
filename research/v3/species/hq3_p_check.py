"""hq3_p_check.py -- HQ3-P: zero-MC checks for P (no O(1) cross-volume plateau).

Sections
  [1] two-point (complete-monotonicity / log-convexity) test on every existing
      zero-momentum F^2 ensemble:  R(t) := c(2t)/c(t)^2  with  c(t)=C(t)/C(0).
  [2] the P-relevant point pair only:  R(Lt/4) = c(Lt/2)/c(Lt/4)^2 .
  [3] same test on EXACT families whose P-truth is known analytically
      (HN1 plateau family = P FALSE; spread family = P FALSE; delta_{1/sqrt L} = P TRUE;
      free field mode sum = P TRUE).  => pass/fail power of the test.
  [4] fixed-tau blindness: fixed-tau correlators cannot separate P-true from P-false.
  [5] log-convexity is forced by w_n>=0 on BOTH geometries (exp on the cylinder,
      cosh on the torus) -- numerical check of the cosh statement.

Estimator for the lattice correlator is copied verbatim from
research/v3/species/w7/run3_apebase.py::correlators (periodic zero-momentum
autocorrelator; nothing in that file is modified or re-run).
NO new MC.  Reads existing *.npz only.  Deterministic.
"""
import os
import sys
import math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = os.path.join(HERE, 'w7')
OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def head(s):
    P("")
    P("=" * 100)
    P(s)
    P("=" * 100)


# ---------------------------------------------------------------- estimator
def correlators(Sm, Lt):
    """verbatim run3_apebase.correlators: periodic zero-momentum autocorrelator."""
    N = Sm.shape[0]
    x = Sm - Sm.mean(axis=0)
    C = np.zeros(Lt)
    for t in range(Lt):
        C[t] = np.mean(np.sum(x * np.roll(x, -t, axis=1), axis=1)) / Lt
    return C


def jack_corr(Sm, Lt, nb=20):
    N = Sm.shape[0]
    C = correlators(Sm, Lt)
    idx = np.array_split(np.arange(N), nb)
    Cs = np.empty((nb, Lt))
    for b in range(nb):
        keep = np.concatenate([idx[j] for j in range(nb) if j != b])
        Cs[b] = correlators(Sm[keep], Lt)
    sig = np.sqrt((nb - 1.0) / nb * np.nansum((Cs - C) ** 2, axis=0))
    return C, sig, Cs


def ratio_R(C, Cs, t):
    """R(t) = c(2t)/c(t)^2  (c = C/C0); jackknife error."""
    def R_of(c):
        return (c[2 * t] / c[0]) / (c[t] / c[0]) ** 2
    R = R_of(C)
    Rs = np.array([R_of(c) for c in Cs])
    Rb = Rs.mean()
    e = math.sqrt((Cs.shape[0] - 1.0) / Cs.shape[0] * np.sum((Rs - Rb) ** 2))
    return R, e


def one(FILE, tag=''):
    d = np.load(FILE)
    S = d['S']
    Lt = int(d['Lt'])
    if S.ndim == 3:
        S = S[:, 0, :]
    beta = float(d['beta'])
    Ls = int(d['Ls']) if 'Ls' in d.files else 8
    N = S.shape[0]
    C, sig, Cs = jack_corr(S, Lt)
    c = C / C[0]
    snr = np.where(sig > 0, np.abs(C) / sig, np.inf)
    name = os.path.basename(FILE)
    P(f"--- {name:<26s} beta={beta:<5g} Ls={Ls} Lt={Lt:<3d} N={N:<6d} C(0)={C[0]:.5e} "
      f"<plaq>={float(np.mean(d['plaq'])):.5f} {tag}")
    P(f"{'t':>3} {'c(t)=C(t)/C(0)':>16} {'+-':>9} {'S/N':>7}   "
      f"{'R(t)=c(2t)/c(t)^2':>18} {'+-':>9}  {'2t':>3}")
    rows = []
    tmax = max(1, Lt // 4)                     # 2t <= Lt/2 : log-convex window
    for t in range(1, tmax + 1):
        R, e = ratio_R(C, Cs, t)
        P(f"{t:>3} {c[t]:>16.5e} {sig[t]/abs(C[0]):>9.2e} {snr[t]:>7.2f}   "
          f"{R:>18.5e} {e:>9.2e}  {2*t:>3d}")
        rows.append((t, c[t], sig[t] / abs(C[0]), snr[t], R, e))
    # P-relevant pair: horizon (Lt/2) vs half-horizon (Lt/4)
    th = Lt // 4
    if Lt % 4 == 0 and th >= 1:
        R, e = ratio_R(C, Cs, th)
        P(f"  [P-pair] horizon Lt/2={Lt//2} : c={c[Lt//2]:.5e} +- {sig[Lt//2]/abs(C[0]):.2e} "
          f"(S/N {snr[Lt//2]:.2f}) ; R(Lt/4)={R:.5e} +- {e:.2e} ; "
          f"c(Lt/4)^2={c[th]**2:.5e}")
    return dict(name=name, beta=beta, Ls=Ls, Lt=Lt, N=N, C0=float(C[0]),
                c=c.tolist(), sig=(sig / abs(C[0])).tolist(), snr=snr.tolist(), rows=rows)


# ---------------------------------------------------------------- [3] exact families
def exact_families():
    head("[3] EXACT FAMILIES: does the two-point test have any discriminating power?")
    P("cylinder C_L(tau)=sum_n w_n e^{-m_n tau};  test R(t)=c(2t)/c(t)^2 >= 1 .")
    P("")
    P("(3a) HN1 plateau family  nu_L = 1/2 delta_{1/L^2} + 1/2 delta_4   (HN1-H1.md:228-232)")
    P("     P is FALSE here (G_c -> 1/2 e^{-1} = 0.1839 , HZ1 0.6 / log [4]).")
    P(f"{'L':>6} {'c(L/2)':>12} {'c(L)':>12} {'R(L/2)=c(L)/c(L/2)^2':>22} {'verdict':>10}")
    bad = 0
    for L in (8, 64, 512, 4096, 65536):
        def c(t):
            return 0.5 * math.exp(-(1.0 / L) * t) + 0.5 * math.exp(-2.0 * t)
        R = c(L) / c(L / 2.0) ** 2
        bad += (R < 1.0)
        P(f"{L:>6} {c(L/2.0):>12.5e} {c(L):>12.5e} {R:>22.5f} {'PASS' if R >= 1 else 'FAIL':>10}")
    P(f"  => violations in 5 windows: {bad}   (P is FALSE, test passes => test cannot see P)")
    P("")
    P("(3b) spread family  nu_L = sum_{j<=L^3} (1/L^3) delta_{(c/L^2)(1+j/L^3)},  c=1")
    P("     P is FALSE here (G_c -> 1 , HZ1 0.2 / log [3]).")
    P(f"{'L':>6} {'c(L/2)':>12} {'c(L)':>12} {'R(L/2)':>14} {'verdict':>10}")
    bad = 0
    for L in (8, 16, 32, 64):
        N = L ** 3
        def c(t, L=L, N=N):
            s = 0.0
            for j in range(1, N + 1):
                mm = math.sqrt((1.0 / L ** 2) * (1 + j / float(N)))
                s += (1.0 / N) * math.exp(-mm * t)
            return s
        R = c(L) / c(L / 2.0) ** 2
        bad += (R < 1.0)
        P(f"{L:>6} {c(L/2.0):>12.5e} {c(L):>12.5e} {R:>14.5f} {'PASS' if R >= 1 else 'FAIL':>10}")
    P(f"  => violations: {bad}   (P is FALSE, test passes)")
    P("")
    P("(3c) P-TRUE family  nu_L = delta_{1/sqrt L}   (m_L L = sqrt L -> infinity)")
    for L in (8, 64, 512):
        def c(t):
            return math.exp(-(1.0 / math.sqrt(L)) * t)
        R = c(L) / c(L / 2.0) ** 2
        P(f"  L={L:<5d} m_L={1/math.sqrt(L):.4f}  c(L/2)={c(L/2.0):.6f}  c(L)={c(L):.6e}  "
          f"R(L/2)={R:.5f}  {'PASS' if R >= 1 else 'FAIL'}")
    P("")
    P("(3d) free massless 4D field on T^3_L x R (mode sum, sqrt c = 2 pi) : P TRUE (HZ1 log [1])")
    for L in (8, 16):
        def c(t, L=L):
            s = 0.0
            for n1 in range(-6, 7):
                for n2 in range(-6, 7):
                    for n3 in range(-6, 7):
                        if n1 == n2 == n3 == 0:
                            continue
                        k = 2 * math.pi * math.sqrt(n1 * n1 + n2 * n2 + n3 * n3) / L
                        s += math.exp(-k * t) / (2 * k)
            return s / L ** 3
        c0 = c(1e-9)
        R = (c(L) / c0) / (c(L / 2.0) / c0) ** 2
        P(f"  L={L:<4d} c(L/2)/c(0)={c(L/2.0)/c0:.5e}  c(L)/c(0)={c(L)/c0:.5e}  R(L/2)={R:.5f}")
    P("")
    P("  ==> the test is satisfied by BOTH P-true and P-false families: it is a consequence")
    P("      of w_n >= 0 only.  Per CLAUDE.md 1.1 (a tool must fail on the free field) it is")
    P("      VACUOUS as a discriminator: free field / 2D / M6 / HN1 / spread all pass.")


# ---------------------------------------------------------------- [4] fixed-tau blindness
def fixed_tau():
    head("[4] FIXED-tau BLINDNESS: no fixed-tau correlator decides P  [T]")
    P("(4a) P FALSE but fixed-tau correlator -> 0:   nu_L = 1/2 delta_{1/L^2} + 1/2 delta_{100}")
    P(f"{'L':>6} {'C(1)':>12} {'C(L)':>12} {'P?':>6}")
    for L in (8, 64, 512, 4096):
        C1 = 0.5 * math.exp(-1.0 / L) + 0.5 * math.exp(-100.0)
        CL = 0.5 * math.exp(-1.0) + 0.5 * math.exp(-100.0 * L)
        P(f"{L:>6} {C1:>12.5e} {CL:>12.5e} {'FALSE':>6}")
    P("")
    P("(4b) P TRUE but fixed-tau correlator -> 1:   nu_L = delta_{1/sqrt L}")
    for L in (8, 64, 512, 4096):
        C1 = math.exp(-1.0 / math.sqrt(L))
        CL = math.exp(-math.sqrt(L))
        P(f"{L:>6} {C1:>12.5e} {CL:>12.5e} {'TRUE':>6}")
    P("")
    P("  => the two classes share the fixed-tau asymptotics (C(1) -> 0 resp. -> 1 in both).")
    P("     The ONLY discriminating scale is tau = L.  Any test built from fixed-tau (or")
    P("     scale-covariant) correlator data is blind to P by construction.")


# ---------------------------------------------------------------- [5] cosh geometry
def cosh_geometry():
    head("[5] the SAME inequality on the TORUS (cosh geometry) is forced too  [T]")
    P("torus: C(tau) = sum_n w_n cosh(m_n(Lt/2-tau))/cosh(m_n Lt/2),  w_n>=0.")
    P("each cosh mode has log C convex in tau on [0,Lt/2]; a nonneg sum of log-convex")
    P("functions is log-convex (Hoelder) => c(2t)>=c(t)^2/c(0) holds for ALL w_n>=0.")
    P(f"{'Lt':>4} {'m':>6} {'R(1)':>10} {'R(2)':>10} {'R(3)':>10} {'R(4)':>10} {'R(5)':>10}")
    for Lt in (12, 16, 18, 20):
        for m in (0.05, 0.3, 1.0):
            def c(t):
                return math.cosh(m * (Lt / 2.0 - t)) / math.cosh(m * Lt / 2.0)
            vals = []
            for t in range(1, Lt // 4 + 1):
                vals.append(c(2 * t) / c(t) ** 2)
            P(f"{Lt:>4} {m:>6.2f} " + " ".join(f"{v:>10.5f}" for v in vals))
    P("  => all entries >= 1 : the torus data CANNOT violate the two-point test as long as")
    P("     the spectral weights are nonnegative.  A measured violation = estimator/positivity")
    P("     defect, not a statement about P.")


def main():
    head("[0] SETUP")
    P("operator: zero-momentum F^2 density, i.e. the per-time-slice spatial plaquette sum")
    P("          S(t) (run3/b01 'lev0 raw', normalised /3) -- same carrier as GP3-PLATEAU.")
    P("estimator: periodic autocorrelator of run3_apebase.correlators (verbatim),")
    P("           delete-one-block jackknife, nb=20.")
    P("dictionary: the finite-volume counterpart of the cylinder horizon L is the maximal")
    P("           temporal separation Lt/2 (torus); c is arbitrary in P (quantifier AAd>=0),")
    P("           so the O(1) factor 2 in 'L <-> Lt/2' only relabels c.")
    P("scale invariance: R(t)=c(2t)/c(t)^2 is invariant under C -> lambda C, so the")
    P("           normalisation of the operator cancels exactly.")

    head("[1][2] THE TWO-POINT TEST ON EVERY EXISTING ENSEMBLE (zero MC)")
    data = []
    D1 = ['b01_b2.2.npz', 'b01_b2.75.npz', 'b01_b3.npz', 'b01_b4.npz',
          'b01_b6.npz', 'b01_b10.npz']
    for f in D1:
        p = os.path.join(W7, f)
        if os.path.exists(p):
            data.append(one(p))
    for f in ['gh2_b2.5_Ls6_Lt12.npz', 'gh2_b2.5_Ls8_Lt12.npz',
              'gh2_b2.5_Ls6_Lt20.npz', 'gh2_b2.5_Ls8_Lt20.npz',
              'gw7_b2.59_Ls10.npz', 'gw7_b2.66_Ls12.npz',
              'gp3_b25_Ls8_Lt16_N24k.npz', 'gp3_b3_Ls8_Lt16.npz',
              'gp3_b3_Ls8_Lt16_N40k.npz', 'gp3_b3_Ls8_Lt16_N60k.npz',
              'gp3_b3_Ls8_Lt16_N100k.npz']:
        p = os.path.join(W7, f)
        if os.path.exists(p):
            data.append(one(p))

    head("[2b] SAME-TEMPORAL-EXTENT VOLUME PAIR (Ls=6 vs 8 at beta=2.5): does the level")
    P("     drop as the spatial volume grows?  (this is the actual cross-volume question)")
    for Lt in (12, 20):
        pa = os.path.join(W7, f'gh2_b2.5_Ls6_Lt{Lt}.npz')
        pb = os.path.join(W7, f'gh2_b2.5_Ls8_Lt{Lt}.npz')
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        da, db = np.load(pa), np.load(pb)
        Sa, Sb = da['S'][:, 0, :], db['S'][:, 0, :]
        Ca, sa, _ = jack_corr(Sa, Lt)
        Cb, sb, _ = jack_corr(Sb, Lt)
        P(f"  Lt={Lt}:  C0(Ls6)={Ca[0]:.5e}  C0(Ls8)={Cb[0]:.5e}  ratio={Cb[0]/Ca[0]:.4f} "
          f"(spatial-sum scaling would give (8/6)^3={ (8/6)**3:.3f})")
        P(f"{'t':>3} {'c_Ls6':>12} {'c_Ls8':>12} {'ratio 8/6':>12}")
        for t in range(0, Lt // 2 + 1):
            P(f"{t:>3} {Ca[t]/Ca[0]:>12.5e} {Cb[t]/Cb[0]:>12.5e} "
              f"{(Cb[t]/Cb[0])/(Ca[t]/Ca[0]):>12.4f}")

    exact_families()
    fixed_tau()
    cosh_geometry()

    head("[6] SUMMARY")
    viol = 0
    tested = 0
    for r in data:
        for (t, c2t, e2, snr2, R, eR) in r['rows']:
            tested += 1
            if R < 1.0 - 2.0 * eR:
                viol += 1
                P(f"  VIOLATION (beyond 2 sigma): {r['name']} t={t} R={R:.4e} +- {eR:.2e} "
                  f"S/N(2t)={snr2:.2f}")
    P(f"  tested point pairs: {tested} ; beyond-2-sigma violations: {viol}")
    P("  (a violation can only be an estimator/positivity defect -- see [5])")
    txt = "\n".join(OUT) + "\n"
    open(os.path.join(HERE, 'hq3_p_check.log'), 'w', encoding='utf-8').write(txt)


if __name__ == '__main__':
    main()
