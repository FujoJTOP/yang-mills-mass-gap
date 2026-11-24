#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
transport_gradient_flow.py  —  Free-field acid tests for the Luescher gradient flow
as a candidate "transport machine" for the 4D Yang-Mills mass gap.

Exact continuum free-field integrals (no lattice job; light, 1-D quads).

Objects (Euclidean, d = 4, flow time t > 0, implicit flow  d_t phi = (d^2 - M^2) phi,
i.e. the L^2-gradient flow of the free massive action  S = 1/2 int (d phi)^2 + M^2 phi^2):

  (1) Flowed energy density (coincident  <(d phi)^2>):
        E_M(t) = int d^4p/(2pi)^4  p^2 e^{-2t(p^2+M^2)} / (p^2+M^2)
               = e^{-2tM^2} int_0^inf ds  e^{-s M^2} / (8 pi^2 (2t+s)^3)
      (the s-representation follows from 1/(p^2+M^2) = int_0^inf e^{-s(p^2+M^2)} ds
       and int d^4p p^2 e^{-alpha p^2}/(2pi)^4 = 1/(8 pi^2 alpha^3)).
      Closed massless form: E_0(t) = 1/(64 pi^2 t^2).

  (2) Flow-time effective mass (log-derivative; the "tilted spectral mean"):
        m_flow^2(t) := -1/2 d/dt ln E_M(t)   (= E[u] under u = spectral variable)
      Exact free value: m_flow^2(t) -> M^2 (M>0),  -> 0 (M=0).

  (3) Candidate scale-invariant monotones:  t*E,  t^2*E,  t^2*dE/dt.

  (4) Flowed two-point function at separation r (exact, 1-D):
        C_M(t,r) = int_0^inf ds e^{-s M^2}/(16 pi^2 (2t+s)^2) * exp(-r^2/(4(2t+s)))
      Effective mass from the r-decay:  m_eff(t,r) := -d/dr ln C_M(t,r).
      Exact asymptotics: m_eff -> M (M>0) for r>>sqrt(t);  m_eff ~ 2/r -> 0 (M=0),
      INDEPENDENT of t at leading order (the flow moves overlaps, not rates).

  (5) k->0 susceptibility / curvature: chi(t) = C(t,k=0) = e^{-2tM^2}/M^2 (M>0),
      so  G''(t) := 1/chi = M^2 e^{2tM^2}; at t=0  G''(0) = M^2 (M>0),  0 (M=0).
      This is the project's H1 object  G''(0)=1/chi.

Run:   py -3.12 transport_gradient_flow.py
"""

import numpy as np
from scipy.integrate import quad

PI2 = np.pi * np.pi


# --------------------------------------------------------------------------
# (1) flowed energy density
# --------------------------------------------------------------------------
def E_flow(t: float, M: float) -> float:
    if M == 0.0:
        return 1.0 / (64.0 * PI2 * t * t)
    f = lambda s: np.exp(-s * M * M) / (8.0 * PI2 * (2.0 * t + s) ** 3)
    val, _ = quad(f, 0.0, np.inf, limit=400)
    return np.exp(-2.0 * t * M * M) * val


# --------------------------------------------------------------------------
# (2) flow-time effective mass  m_flow^2(t) = -1/2 d/dt ln E
# --------------------------------------------------------------------------
def m_flow_sq(t: float, M: float, h: float = 1e-4) -> float:
    lnE_p = np.log(E_flow(t + h, M))
    lnE_m = np.log(E_flow(t - h, M)) if (t - h) > 1e-7 else np.log(E_flow(1e-7, M))
    return -0.5 * (lnE_p - lnE_m) / ((t + h) - (t - h if (t - h) > 1e-7 else 1e-7))


def dE_dt(t: float, M: float, h: float = 1e-4) -> float:
    return (E_flow(t + h, M) - E_flow(t - h, M)) / (2.0 * h)


# --------------------------------------------------------------------------
# (4) flowed two-point function and its r-decay effective mass
# --------------------------------------------------------------------------
def C_flow(t: float, r: float, M: float) -> float:
    f = lambda s: (np.exp(-s * M * M) / (16.0 * PI2 * (2.0 * t + s) ** 2)
                   * np.exp(-r * r / (4.0 * (2.0 * t + s))))
    # split large-s tail so quad on [0,inf) is stable; exponential in s dominates
    val, _ = quad(f, 0.0, np.inf, limit=400)
    return val


def m_eff_r(t: float, r: float, M: float, h: float = 1e-3) -> float:
    return -(np.log(C_flow(t, r + h, M)) - np.log(C_flow(t, r - h, M))) / (2.0 * h)


# --------------------------------------------------------------------------
# (5) k->0 susceptibility / curvature
# --------------------------------------------------------------------------
def Gpp0_of_t(t: float, M: float) -> float:
    # G''(t) = 1/chi(t), chi(t)=C(t,k=0)=e^{-2tM^2}/M^2  (M>0); 0 for M=0
    if M == 0.0:
        return 0.0
    return M * M * np.exp(2.0 * t * M * M)


def banner(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


def main():
    banner("(A) flowed energy density and the candidate scale-invariant monotones")
    print("    E_M(t) and  tE,  t^2 E,  t^2 dE/dt,   for M=0 (gapless) and M=1")
    print(f"    {'t':>6} | {'E0(t)':>12} {'t^2 E0':>10} {'tE0':>10} {'t^2 dE0/dt':>12}"
          f" | {'E1(t)':>12} {'t^2 E1':>10}")
    for t in [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]:
        E0 = E_flow(t, 0.0)
        E1 = E_flow(t, 1.0)
        print(f"    {t:6.2f} | {E0:12.4e} {t*t*E0:10.4e} {t*E0:10.4e} "
              f"{t*t*dE_dt(t,0.0):12.4e} | {E1:12.4e} {t*t*E1:10.4e}")
    print("\n    KEY: for the FREE MASSLESS field  t^2 E0(t) = 1/(64 pi^2) = "
          f"{1/(64*PI2):.6e}  at EVERY t  (a positive constant).")
    print("         Any lower bound built from a positive energy functional is therefore")
    print("         satisfied by the GAPLESS free field  =>  acid test FAILS (vacuous).")


    banner("(B) flow-time effective mass  m_flow^2(t) = -1/2 d/dt ln E   (DIRECTION)")
    print(f"    {'t':>6} | {'M=0':>12} {'M=1':>12}    (should be >= M^2 and decreasing to M^2)")
    prev0, prev1 = None, None
    mono0, mono1 = True, True
    for t in [0.1, 0.3, 1.0, 2.0, 4.0, 8.0, 16.0]:
        m0 = m_flow_sq(t, 0.0)
        m1 = m_flow_sq(t, 1.0)
        if prev0 is not None and m0 > prev0 + 1e-9:
            mono0 = False
        if prev1 is not None and m1 > prev1 + 1e-9:
            mono1 = False
        prev0, prev1 = m0, m1
        print(f"    {t:6.2f} | {m0:12.5f} {m1:12.5f}")
    print(f"\n    decreasing? M=0: {mono0}   M=1: {mono1}")
    print("    M=1: m_flow^2(t) >= 1 = M^2 for all t  (UPPER bound on the true mass^2),")
    print("         decreasing to M^2  ->  the finite-t effective mass OVER-estimates M^2.")


    banner("(C) flowed correlator r-decay effective mass m_eff(t,r) = -d/dr ln C")
    print("    massless (M=0): should track 2/r (power, NOT exponential) and be")
    print("                    t-INDEPENDENT -> true rate 0 (tool dies: no gap).")
    print("    massive  (M=1): should plateau at M=1  for r >> sqrt(t), t-INDEPENDENT.")
    ts = [0.1, 1.0, 3.0]
    rs = [2.0, 3.0, 4.0, 6.0, 8.0]
    for M in [0.0, 1.0]:
        print(f"\n    --- M = {M} ---")
        hdr = "    r  | " + " ".join(f"t={t:<5}" for t in ts) + "   |  2/r (M=0 ref)"
        print(hdr)
        for r in rs:
            cols = " ".join(f"{m_eff_r(t, r, M):6.3f}" for t in ts)
            print(f"    {r:4.1f}| {cols}   |  {2.0/r:6.3f}")
    print("\n    -> massive rate = 1.0 regardless of t;  massless rate = 2/r -> 0 regardless of t.")
    print("       The RATE is t-independent (flow moves overlaps c_n(t), never masses m_n).")


    banner("(D) k->0 curvature G''(t)=1/chi(t)  (the project's H1 object)")
    print(f"    {'t':>6} | {'G''_M0':>10} {'G''_M1':>12}")
    for t in [0.0, 0.5, 1.0, 2.0]:
        print(f"    {t:6.2f} | {Gpp0_of_t(t,0.0):10.4f} {Gpp0_of_t(t,1.0):12.5f}")
    print("\n    At t=0:  G''(0) = M^2 (massive) , 0 (massless)  -- the exact H1 free-field value.")
    print("    The flow multiplies it by e^{+2tM^2}>0 but the k->0 multiplier is 1, so the")
    print("    COARSE curvature G''(k->0) is flow-invariant;  G''=1/chi only UPPER-bounds m1.")


    banner("SUMMARY (machine-readable)")
    t = 1.0
    print(f"  free MASSIVE (M=1), t={t}:  m_flow^2 = {m_flow_sq(t,1.0):.5f}  (>= M^2=1, upper)")
    print(f"  free MASSLESS      , t={t}:  m_flow^2 = {m_flow_sq(t,0.0):.5f}  (finite >0, then ->0)")
    print(f"  free MASSLESS t^2 E(t) = {t*t*E_flow(t,0.0):.6e} = 1/(64 pi^2) = {1/(64*PI2):.6e}")
    print("  => positivity monotones are positive on the gapless field (acid test fails);")
    print("  => the rate-level 'effective mass' is an UPPER bound (m1^2 <= m_eff^2), t-independent.")


if __name__ == "__main__":
    main()
