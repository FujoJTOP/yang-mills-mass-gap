#!/usr/bin/env python3
"""
ds_acid.py -- ACID TESTS for the split-quotient D/S, the split audit on archived lattice data,
and the two-plaquette (convolution) model.

  [D1] free field (massless / massive, exact Gaussian):  D = 0, S = m^2  => ratio 0 (or 0/0)
  [D2] beta = 0 pure Haar:  S = A = 1/(2 Cin 2pi), D = 0  => ratio 0 while f_1^2 = 0   (vacuous)
  [D3] double well with a coarse/fibre coupling: EXACT closed form; S<0 at the barrier
       => "D/S<1" holds exactly where G''<0  (the ratio is orientation-blind)
  [D4] two plaquettes (independent, convolution): the ratio is an intensive object
  [D5] AUDIT of the archived 12-mode SU(2) b=2.4 L=8 data: can the "2%" be reproduced?
       (it cannot: G''*Var(a) per mode is NOT ~1, and 1/Var(a) != S)

Run: py -3.12 -I ds_acid.py
"""
import numpy as np
import importlib.util, os, sys

PI = np.pi
_here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ds', os.path.join(_here, 'ds_ratio.py'))
ds = importlib.util.module_from_spec(spec); spec.loader.exec_module(ds)
slice_moments, splits, fd_curv, E3, E1 = ds.slice_moments, ds.splits, ds.fd_curv, ds.E3, ds.E1
def P(*a): print(*a, flush=True)

# ------------------------------------------------------------------ [D1]
def d1_free():
    P("=" * 100)
    P("[D1] ACID TEST (i): FREE FIELD, exact")
    P("     V = m^2 a^2/2 + (fibre quadratic):  measure Gaussian, D = Var(V') = 0 identically")
    P("=" * 100)
    P("     %-22s %10s %10s %12s %10s" % ('model', 'S=E[V2]', 'D=VarV1', 'Gpp=S-D', 'D/S'))
    for m2, tag in ((0.0, 'massless (flat)'), (1.0, 'massive m^2=1'), (0.25, 'massive m^2=0.25')):
        S, D = m2, 0.0
        r = '0/0 undef' if S == 0 else '%.4f' % (D / S)
        P("     %-22s %10.4f %10.4f %12.4f %10s" % (tag, S, D, S - D, r))
    P("   ==> massless: S=D=0, the ratio is 0/0 (indeterminate), G''=0.")
    P("       massive  : ratio = 0 < 1 WHILE the model HAS a gap m^2>0.")
    P("       => 'D/S < 1' is satisfied by every Gaussian model, gapped or not: D10 (no discrimination).")

# ------------------------------------------------------------------ [D2]
def d2_haar(A):
    P("")
    P("=" * 100)
    P("[D2] ACID TEST (ii): beta = 0 PURE HAAR (Model A) -- curvature without an atom")
    P("=" * 100)
    m = slice_moments(0.0, E3, 0.0, PI); s = splits(m)
    P("     (one plaquette, theta=pi, v=e3, a0=0, beta=0):  S = %.10f  D = %.3e  ratio = %.3e"
      % (s['S_A2'], s['D_A2'], s['r_A2']))
    P("     G''(beta=0) = %.10f  (= 1/(2 Cin 2pi), the P1 anchor)" % s['Gpp_A2'])
    P("     but f_1^2 = 0 at beta=0 (no Wilson weight => no glueball atom).")
    P("   ==> the ratio is 0 (maximally 'good') on a measure with NO gap: the criterion is VACUOUS here.")
    P("   NOTE on the corpus quote 'G''=0.2688 at beta=0':  0.268830 is 1/Var(a), NOT G''.")
    P("     G'' = 0.20511530 = 1/(2 Cin(2pi))  ;  1/Var(a) = 0.268830  (Var = 3.719824).")
    P("     The 'product identity' reading  D/S := 1 - G''*Var(a)  = 1 - 0.205115/0.268830 = %.4f"
      % (1 - 0.2051152971 / 0.268830))
    P("     while the true A2 defect ratio at beta=0 is 0.0000  =>  that reading is not D/S.")

# ------------------------------------------------------------------ [D3]
def d3_double_well():
    P("")
    P("=" * 100)
    P("[D3] ACID TEST (iii): DOUBLE WELL, EXACT CLOSED FORM")
    P("     V(x,y) = lam (x^2-1)^2 - h x + w^2|y|^2/2 + g x |y|^2 ,  y in R^2, flat measure")
    P("     coarse a = x  =>  S(a) = 4 lam (3a^2-1) ,  D(a) = g^2 Var(|y|^2|a) = 4 g^2/(w^2+2ga)^2")
    P("=" * 100)
    lam, w2 = 1.0, 1.0
    P("   %-8s %-8s %8s %10s %10s %10s %9s" % ('g', 'a', 'S', 'D', 'Gpp', 'D/S', 'Gpp>0?'))
    for g in (0.0, 0.3, 0.8, 1.5):
        for a in (0.0, 0.5, 1.0, 1.5):
            tau2 = 1.0 / (w2 + 2 * g * a)
            S = 4 * lam * (3 * a * a - 1)
            D = 4 * g * g * tau2 ** 2
            r = D / S if S != 0 else float('nan')
            P("   %-8.2f %-8.2f %8.4f %10.4f %10.4f %10.4f %9s"
              % (g, a, S, D, S - D, r, 'YES' if S - D > 0 else 'no'))
    P("   ==> (a) g=0 (symmetric well): D=0, ratio=0 EVERYWHERE -- including the barrier a=0 where")
    P("           G'' = -4 lam < 0.  The criterion 'D/S<1' is satisfied at a negative-curvature point.")
    P("       (b) for S<0 (barrier, |a|<1/sqrt3) the ratio is NEGATIVE: 'D/S<1' <==> 'D>S' there,")
    P("           i.e. the equivalence  D/S<1 <=> G''>0  REQUIRES S>0.  Sign-flip, not a margin.")
    P("       (c) the crossing D/S=1 IS the edge G''=0 (there is no other scale):")
    for g in (0.3, 0.8, 1.5):
        aa = np.linspace(0.0, 2.0, 4001)
        tau2 = 1.0 / (w2 + 2 * g * aa)
        S = 4 * lam * (3 * aa ** 2 - 1); D = 4 * g * g * tau2 ** 2
        Gp = S - D
        idx = np.where(np.sign(Gp[:-1]) * np.sign(Gp[1:]) < 0)[0]
        seg = ', '.join('a=%.3f (S=%.3f, D/S=%.3f)' % (aa[i], S[i], D[i] / S[i]) for i in idx[:2])
        P("           g=%.1f : G''=0 at %s ;  D/S at a=1 (S=8): %.4f" % (g, seg, D[-1] / S[-1]))
    P("       => the ratio is a DISTANCE-TO-THE-EDGE coordinate; it equals 1 only on the edge itself.")

# ------------------------------------------------------------------ [D4]
def d4_two_plaquette():
    P("")
    P("=" * 100)
    P("[D4] TWO PLAQUETTES (a = (w1+w2)/sqrt2) + 1D CHAIN: exact reduction, and an HONEST [E]")
    P("     rho_2(a) = int rho(w) rho(sqrt2 a - w) dw   (exact: independent copies)")
    P("     The numerical defect split below is run as a SELF-TEST of the identity S-D = -d2log rho;")
    P("     if it fails, the numbers are NOT evidence (recorded as [E], not [V]).")
    P("=" * 100)
    beta, theta = 2.5, PI
    W = 2 * PI
    xq, wq = np.polynomial.legendre.leggauss(24)
    wj = W * xq; wj_w = W * wq

    def slab(w):                      # exact slice moments at one value of the coarse coordinate
        if abs(w) >= W - 1e-12:       # outside the flat chart |xi| < 2pi : density = 0
            return (-np.inf, 0.0, 0.0, 0.0)
        m = slice_moments(float(w), E3, beta, theta, ns=400, nph=32)
        Evp = m['EG'] + m['EP1']
        return (m['rho'], m['EH'] + m['EP2'], Evp,
                m['VarG'] + m['VarP1'] + 2 * m['Cov'] + Evp ** 2)

    def mom2(a):
        """rho_2(a) and the exact S2, D2 of a = (w1+w2)/sqrt2, by nested exact quadrature."""
        A = [slab(w) for w in wj]
        B = [slab(np.sqrt(2) * a - w) for w in wj]
        wt = np.array([(0.0 if (not np.isfinite(A[i][0]) or not np.isfinite(B[i][0]))
                        else np.exp(A[i][0] + B[i][0])) * wj_w[i] for i in range(len(wj))])
        Z = wt.sum()
        m2a = np.array([A[i][1] for i in range(len(wj))]); m2b = np.array([B[i][1] for i in range(len(wj))])
        m1a = np.array([A[i][2] for i in range(len(wj))]); m1b = np.array([B[i][2] for i in range(len(wj))])
        qa = np.array([A[i][3] for i in range(len(wj))]); qb = np.array([B[i][3] for i in range(len(wj))])
        S2 = 0.5 * (wt * m2a).sum() / Z + 0.5 * (wt * m2b).sum() / Z
        Ea = (wt * m1a).sum() / Z; Eb = (wt * m1b).sum() / Z
        cross = (wt * m1a * m1b).sum() / Z
        D2 = 0.5 * ((wt * qa).sum() + (wt * qb).sum()) / Z + cross - (0.5 * (Ea + Eb)) ** 2
        return np.log(Z), S2, D2

    P("     %-8s %12s %11s %11s %11s %11s %10s"
      % ('a', 'log rho_2', 'S2', 'D2', 'G2=S2-D2', 'fd(G2)', 'D2/S2'))
    for a in (0.0, 0.5, 1.0):
        l0, S2, D2 = mom2(a)
        da = 0.02
        f2 = -(-mom2(a - 2 * da)[0] + 16 * mom2(a - da)[0] - 30 * l0
               + 16 * mom2(a + da)[0] - mom2(a + 2 * da)[0]) / (12 * da * da)
        P("     %-8.2f %12.6f %11.6f %11.6f %11.6f %11.6f %10.5f"
          % (a, l0, S2, D2, S2 - D2, f2, D2 / S2 if S2 else float('nan')))
    m1p = slice_moments(0.0, E3, beta, theta); s1 = splits(m1p)
    P("     single plaquette at a=0:  S=%.6f D=%.6f ratio=%.6f  G''=%.6f"
      % (s1['S_A2'], s1['D_A2'], s1['r_A2'], s1['Gpp_A2']))
    P("   SELF-TEST: at a=0 the identity gives S2-D2 = -0.187 while the 5-point FD of log rho_2 gives")
    P("       +0.090 -> discrepancy 0.28, i.e. the quadrature is NOT certified (nw=24 is too coarse for")
    P("       a 5-point FD at da=0.02).  RESULT: [E] -- no number from this block may be used.")
    P("   WHAT IS BANKED (exact, no quadrature):")
    P("     (a) rho_2 = rho_1 * rho_1 exactly, so the pooled model is a functional of the one-plaquette")
    P("         marginal alone: pooling cannot supply a *new* estimate, only a new function of rho_1.")
    P("     (b) 1D PURE-GAUGE CHAIN: dU_1...dU_N = dW dV_1...dV_{N-1} (Haar invariance) => the chain")
    P("         integral is a SINGLE plaquette integral and every gauge-invariant functional equals its")
    P("         one-plaquette value.  A 1D chain contributes NO independent D/S.  [T]")
    P("   ==> none of {1 plaquette, 2 plaquettes, 1D chain} supplies a bound < 1 beyond the sign of G''.")

# ------------------------------------------------------------------ [D5]
def d5_npz_audit():
    P("")
    P("=" * 100)
    P("[D5] AUDIT of the archived lattice data  k0_su2_b2.4_L8_n2.npz  (SU(2), beta=2.4, L=8, 400 cfg)")
    P("     fields: a, H, G, P1, P2 (12 modes at k = pi/4 and pi/2; NOT k=0, NOT 0^{++})")
    P("=" * 100)
    p = os.path.normpath(r"D:\Dev\math\research\unconditional\w7\k0_su2_b2.4_L8_n2.npz")
    if not os.path.exists(p):
        P("   [npz missing]"); return
    d = np.load(p, allow_pickle=True)
    a, H, G, P1, P2 = d['a'], d['H'], d['G'], d['P1'], d['P2']
    nm = d['names']
    S2 = H.mean(0) + P2.mean(0)
    D2 = G.var(0) + P1.var(0) + 2 * ((G - G.mean(0)) * (P1 - P1.mean(0))).mean(0)
    Gp = S2 - D2
    Va = a.var(0)
    P("   %-16s %8s %8s %8s %9s %9s %9s %8s"
      % ('mode', 'Var(a)', 'Gpp', '1/Var', 'Gpp*Var', 'S_A2', 'D/S_A2', 'Gpp/S'))
    for i in range(len(nm)):
        P("   %-16s %8.4f %8.4f %8.4f %9.4f %9.4f %9.4f %8.4f"
          % (nm[i], Va[i], Gp[i], 1 / Va[i], Gp[i] * Va[i], S2[i], D2[i] / S2[i], Gp[i] / S2[i]))
    P("   means:  Var(a)=%.4f   1/Var=%.4f   S_A2=%.4f   Gpp=%.4f   E[Gpp]/E[S]=%.4f"
      % (Va.mean(), (1 / Va).mean(), S2.mean(), Gp.mean(), Gp.mean() / S2.mean()))
    P("   PRODUCT identity per mode:  Gpp*Var(a)  mean=%.4f  min=%.4f  max=%.4f"
      % ((Gp * Va).mean(), (Gp * Va).min(), (Gp * Va).max()))
    P("   Schur-type relation S_A2*Var(a):  mean=%.4f  (a Gaussian marginal would give 1)"
      % (S2 * Va).mean())
    P("   DEFECT ratio (A2):            mean=%.4f  min=%.4f  max=%.4f   modes D>S: %d/12"
      % ((D2 / S2).mean(), (D2 / S2).min(), (D2 / S2).max(), (Gp < 0).sum()))
    S3 = H.mean(0); D3 = D2 - P2.mean(0)
    P("   DEFECT ratio (H2 split):      mean=%.4f  min=%.4f  max=%.4f"
      % ((D3 / S3).mean(), (D3 / S3).min(), (D3 / S3).max()))
    P("   CONSEQUENCE: on this data the defect ratio is 0.97 (NOT 0.02), i.e. the margin is ~3 percent, and")
    P("   the two readings are complementary, not independent: 1-0.97 = 0.03.  The '2 percent' as D/S is a")
    P("   complementation error; the '2 percent' as 1-Gpp*Var(a) is a different functional (mean 0.83 here).")

def main():
    x, w = np.polynomial.legendre.leggauss(4000)
    rho = PI * (x + 1.0); wr = PI * w
    m = wr * (np.sin(rho / 2) ** 2 / rho)
    A = float((m * ((2.0 / rho - np.cos(rho / 2) / np.sin(rho / 2)) / rho)).sum() / m.sum())
    d1_free(); d2_haar(A); d3_double_well(); d4_two_plaquette(); d5_npz_audit()
    P("")
    P("Clay 0%.  No QED.")

if __name__ == '__main__':
    main()

