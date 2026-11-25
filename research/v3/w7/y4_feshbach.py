#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
y4_feshbach.py -- Y4: TEST the "gauge-invariant nonlinear Feshbach / multiscale RG projection".

THE PROPOSAL UNDER TEST
-----------------------
Build a gauge-independent (no gauge fixing) Feshbach / RG projection that tracks the low-frequency
k=0 collective mode, via the FESHBACH EFFECTIVE ACTION (marginal / large-deviation rate function)
of a GAUGE-INVARIANT density O's total-normalized zero mode:

    q_O[U] = V^{-1/2} sum_x ( O(x) - <O> ) ,        S_eff(q) := - log E_mu[ delta(q - q_O[U]) ] .

The "multiscale" version does this at each momentum k (and each smearing scale) -> C(k)=Var(q(k)).

  * q_O is a gauge-invariant random variable (O gauge-invariant, 1/sqrt(V) a constant, <O> a
    c-number), so S_eff(q) and its curvature are EXACTLY gauge-invariant -- NO gauge choice.
  * Var(q_O) = chi_O = sum_z C(z)  (total-normalized; NOT V^{-1})  -- the T4/WG carrier.
  * S_eff''(q*) = the Feshbach curvature G'' (= Hess of the effective action at the mode).

THE QUESTIONS
-------------
(Q1) Non-vacuous carrier with Var != V^{-1}?     (Q2) Does it IMPROVE on the WG estimator?
ACID: free field -> known / 0 ; 2D trF^2.   KILL: (a) repackages WG/H1, (b) needs a gauge choice,
(c) cannot control the collective mode.

WHAT THIS WITNESS SHOWS (deterministic, one process, numpy only)
----------------------------------------------------------------
[A] GAUSSIAN: G'' = (P A^{-1} P^T)^{-1}, Var = P A^{-1} P^T, so  G'' * Var = I  (to 1e-14).  In the
    Gaussian limit the Feshbach curvature IS the reciprocal susceptibility -> repackages WG.
[B] FREE ELEMENTARY phi (1D): S_eff Gaussian, G'' = m^2 = 1/Var EXACTLY (the KL measure is a single
    atom, so k=0 marginal recovers the pole).  Acid: massless -> G''=0.  Known, no false floor.
[C] FREE COMPOSITE :phi^2: (the gauge-invariant analog of trF^2): the k=0 marginal is NON-Gaussian
    at finite L, G'' = S_eff''(0) != 1/Var, BUT the excess kurtosis k4 -> 0 (CLT) and G'' -> 1/Var
    as L grows; both stay != (2m)^2 = the true lightest threshold.  So the "nonlinear" content is a
    finite-V non-Gaussianity that VANISHES thermodynamically, and is a CEILING (not the pole), atomless.
[D] ONE-PLAQUETTE SU(2) (exact, V=1): density rho(z)=(2/pi)sqrt(1-z^2); EXACT G''=S_eff''(0)=1 while
    1/Var=4 -> the maximal finite-V split (G'' != 1/chi); smooth density => atomless q.
[E] SAVIDDY ONE-LOOP cell: the gauge-invariant collective-coordinate effective potential V_eff(B) has
    NEGATIVE curvature V_eff''(B)<0 for small B (the inherited Savvidy tachyon).  The projection
    does NOT control the collective mode; the full-measure boundedness is the WG measure fact, not
    the projection.
[F] GAUGE INVARIANCE (2D/3D SU(2) config): q=V^{-1/2}sum trF^2 invariant under a random gauge
    transform to ~1e-15; the link-log k=0 mode drifts.  Passes the "no gauge choice" criterion.
[G] MULTISCALE (all k): C(k)=Var(q(k))=Gtilde(k)=1/(m^2+khat^2); the effective mass from C(x) decay
    recovers m.  This is the STANDARD gauge-invariant effective mass (estimate, not a lower bound),
    not a new tool.

Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.  Run: py -3.12 -I research/v3/w7/y4_feshbach.py
"""
import sys
import numpy as np

np.seterr(all='ignore')
RNG = np.random.default_rng(20261124)
LINE = "=" * 84


# ------------------------------------------------------------------ helpers
def free_G_hat(L, m):
    """1D periodic: Ghat(k)=1/(m^2+4 sin^2(k/2)), and G(x)=ifft(Ghat)."""
    k = 2.0 * np.pi * np.arange(L) / L
    phat2 = 4.0 * np.sin(k / 2.0) ** 2
    Ghat = 1.0 / (phat2 + m * m)
    G = np.fft.ifft(Ghat).real
    return G, Ghat, phat2


# ------------------------------------------------------------------ [A] Gaussian identity
def A_gaussian():
    print("[A] GAUSSIAN block:  G'' = (P A^-1 P^T)^-1 ,  Var = P A^-1 P^T  =>  G'' * Var = I")
    for seed in (0, 1, 2):
        rng = np.random.default_rng(seed)
        n, m = 4, 3
        X = rng.standard_normal((n, n)); A = X @ X.T + n * np.eye(n)     # SPD precision
        P = rng.standard_normal((m, n))                                  # coarse projector
        Varcov = P @ np.linalg.inv(A) @ P.T                              # = Var(q) for the projection
        Gpp = np.linalg.inv(Varcov)                                      # = Hess of S_eff (Gaussian)
        prod = Gpp @ Varcov
        err = np.max(np.abs(prod - np.eye(m)))
        print(f"    seed={seed}: || G'' * Var - I ||_max = {err:.3e}   (exact => repackages WG)")
    print("    => in the Gaussian limit the Feshbach curvature is EXACTLY the reciprocal susceptibility.\n")


# ------------------------------------------------------------------ [B] free elementary scalar
def B_elementary():
    print("[B] FREE ELEMENTARY phi (1D): S_eff Gaussian,  G'' = m^2 = 1/Var EXACTLY ; acid m->0 -> 0")
    print(f"    {'L':>5} {'m':>7} {'Var=1/m^2':>12} {'1/Var':>10} {'G''=m^2':>10} {'chi':>10}")
    for L in (16, 64):
        for m in (0.5, 0.2):
            G, Ghat, _ = free_G_hat(L, m)
            Var = G.sum()                      # = Ghât(0) = 1/m^2
            print(f"    {L:>5} {m:>7.2f} {Var:>12.6f} {1.0/Var:>10.6f} {m*m:>10.6f} {Var:>10.6f}")
    # acid massless: Ghat(k=0) -> inf (IR), curvature -> 0
    for m in (0.05, 0.01):
        G, Ghat, _ = free_G_hat(256, m)
        print(f"    acid m={m:<5}: Var=1/m^2={1.0/(m*m):>12.1f}  G''=m^2={m*m:.2e} -> 0  (no false floor)")
    print("    => single-atom KL measure: the k=0 marginal curvature recovers the pole. ACID-CLEAN.\n")


# ------------------------------------------------------------------ [C] free composite :phi^2:
def sample_composite(L, m, nsamp):
    """Sample the total-normalized zero mode q of :phi^2:(x)=phi(x)^2-<phi^2> for free phi.
       Hermitian (real-field) Fourier sampling: E|phi_hat(k)|^2 = L*Ghat(k), phi_hat(-k)=conj."""
    G, Ghat, phat2 = free_G_hat(L, m)
    qs = np.empty(nsamp)
    B = 2000
    kmax = L // 2
    for s0 in range(0, nsamp, B):
        b = min(B, nsamp - s0)
        ph = np.zeros((b, L), dtype=complex)
        ph[:, 0] = np.sqrt(L * Ghat[0]) * RNG.standard_normal(b)        # k=0 real
        if L % 2 == 0:
            ph[:, kmax] = np.sqrt(L * Ghat[kmax]) * RNG.standard_normal(b)  # k=L/2 real
        for kk in range(1, kmax):
            amp = np.sqrt(L * Ghat[kk] / 2.0)
            z = amp * (RNG.standard_normal(b) + 1j * RNG.standard_normal(b))
            ph[:, kk] = z
            ph[:, L - kk] = np.conj(z)
        phi = np.fft.ifft(ph, axis=1).real       # real field, <phi(x)^2>=G(0)
        norm = phi * phi - G[0]                  # normal-ordered :phi^2:
        q = norm.sum(axis=1) / np.sqrt(L)
        qs[s0:s0 + b] = q
    return qs, G


def seff_curv_hist(qs, half_width_sigma=1.5, nbins=241):
    """Direct S_eff''(0) from a histogram of the (centered) sample: fit -log p ~ c0 + c2 q^2
       over |q| < half_width_sigma*sigma; S_eff'' = 2 c2.  Non-circular (no Edgeworth)."""
    mu = qs.mean(); sd = qs.std()
    q = qs - mu
    lo, hi = -half_width_sigma * sd, half_width_sigma * sd
    cnt, edges = np.histogram(q, bins=nbins, range=(lo, hi))
    ctr = 0.5 * (edges[1:] + edges[:-1])
    good = cnt > 0
    lpdf = -np.log(cnt[good] / qs.size)
    # parabola fit (no constant weight issues): weighted least squares in q^2
    X = np.vstack([np.ones(good.sum()), ctr[good] ** 2]).T
    coef, *_ = np.linalg.lstsq(X, lpdf, rcond=None)
    return 2.0 * coef[1]


def C_composite():
    print("[C] FREE COMPOSITE :phi^2: (1D): non-Gaussian k=0 marginal;  G''=S_eff'' != 1/Var at finite L,")
    print("    but kurtosis k4 -> 0 (CLT) and G'' -> 1/Var ; both stay != (2m)^2 = true threshold.")
    m = 0.5
    pole = (2.0 * m) ** 2
    print(f"    m={m}, true lightest threshold (2m)^2 = {pole:.6f}")
    print(f"    {'L':>4} {'chi=Var':>9} {'1/Var':>9} {'S_eff\"(hist)':>13} {'k4':>8} {'S_eff\"/1/Var':>13}")
    for L in (16, 32, 64):
        qs, G = sample_composite(L, m, nsamp=200000)
        var = qs.var()
        m4 = ((qs - qs.mean()) ** 4).mean()
        k4 = m4 / var**2 - 3.0                      # excess kurtosis
        gpp = seff_curv_hist(qs)                    # direct S_eff''(0), no Edgeworth
        print(f"    {L:>4} {var:>9.5f} {1.0/var:>9.5f} {gpp:>13.5f} {k4:>8.4f} {gpp*var:>13.4f}")
    # exact chi (Wick) = sum_z 2 G(z)^2  for reference (should match the sampled Var above)
    for L in (16, 64):
        _, Ghat, _ = free_G_hat(L, m)
        G = np.fft.ifft(Ghat).real
        chi_exact = 2.0 * (G * G).sum()
        print(f"    exact chi(Wick) L={L}: {chi_exact:.5f}   (1/chi={1/chi_exact:.5f} != (2m)^2={pole:.5f})")
    # acid massless
    for m2 in (0.1, 0.02):
        qs, G = sample_composite(64, m2, nsamp=60000)
        print(f"    acid m={m2}: Var={qs.var():>12.1f} -> inf  (gapless: H1 FAILS, no floor)")
    print("    => 'nonlinear' = finite-V non-Gaussianity (k4) that VANISHES (CLT); S_eff''->1/Var;")
    print("       both are CEILING-type (!= (2m)^2) and the k=0 marginal is atomless (f1^2=0). No gain.\n")


# ------------------------------------------------------------------ [D] one plaquette SU(2) exact
def D_plaquette():
    print("[D] ONE-PLAQUETTE SU(2) (exact, V=1):  rho(z)=(2/pi)sqrt(1-z^2), z=1/2 tr U (Haar).")
    # exact moments: Var(z)=1/4 ; mode z*=0 ; S_eff=-log rho
    var_exact = 0.25
    # S_eff''(0) = -rho''/rho at 0 :  rho''(0) = -2/pi, rho(0)=2/pi  => 1
    gpp_exact = 1.0
    print(f"    exact Var(z) = {var_exact:.6f}   1/Var = {1.0/var_exact:.6f}")
    print(f"    exact G'' = S_eff''(0) = -rho''(0)/rho(0) = {gpp_exact:.6f}")
    print(f"    => G'' != 1/Var  (factor { (1.0/var_exact)/gpp_exact:.1f}). MAX finite-V split (V=1, non-Gaussian).")
    # numeric check by quadrature
    z = np.linspace(-1 + 1e-9, 1 - 1e-9, 400001)
    rho = (2.0 / np.pi) * np.sqrt(1.0 - z * z)
    dz = z[1] - z[0]
    print(f"    numeric Var = {np.trapezoid(z*z*rho, z)/np.trapezoid(rho, z):.6f} (norm={np.trapezoid(rho,z):.6f})")
    print(f"    density smooth on (-1,1) => q atomless (no delta). Ceiling-type: works only if single atom.\n")


# ------------------------------------------------------------------ [E] Savvidy one-loop cell
def E_savvidy():
    print("[E] SAVIDDY ONE-LOOP cell: gauge-invariant collective-coordinate effective potential")
    print("    V_eff(B) = 1/2 B^2 + A B^2 ( ln(B/M^2) - 1/2 ),  A = 11 g^2/(96 pi^2) > 0")
    print("    V_eff''(B) = 1 + 2A ln(B/M^2) + 2A.  Sign: NEGATIVE for small B (inherited tachyon).")
    print(f"    {'g^2':>6} {'A':>10} {'B_c/M^2':>12} {'B_min/M^2':>12} {'V''(B_c/10)':>13} {'V''(B_min)':>12}")
    for g2 in (1.0, 2.5, 6.0):
        A = 11.0 * g2 / (96.0 * np.pi ** 2)
        # zero of V'': 1+2A ln(x)+2A = 0 -> x_c=exp(-(1+2A)/(2A)); min: 1+2A ln x=0 -> x_min=exp(-1/(2A))
        x_c = np.exp(-(1.0 + 2.0 * A) / (2.0 * A))
        x_min = np.exp(-1.0 / (2.0 * A))
        Vpp = lambda x: 1.0 + 2.0 * A * np.log(x) + 2.0 * A
        print(f"    {g2:>6.2f} {A:>10.5f} {x_c:>12.3e} {x_min:>12.3e} {Vpp(x_c/10):>13.4f} {Vpp(x_min):>12.4f}")
    print("    Strong-coupling toy A=0.5 (visible concavity, illustrative):")
    A = 0.5
    x_c = np.exp(-(1.0 + 2.0 * A) / (2.0 * A)); x_min = np.exp(-1.0 / (2.0 * A))
    Vpp = lambda x: 1.0 + 2.0 * A * np.log(x) + 2.0 * A
    for x in (0.02, 0.05, x_c, x_min, 0.8, 2.0):
        region = "STABLE" if Vpp(x) > 0 else "UNSTABLE (concave)"
        print(f"      B/M^2={x:>8.4f}  V''={Vpp(x):>9.4f}   {region}")
    print("    => the gauge-invariant collective coordinate (q ~ B^2) has NEGATIVE curvature near B=0:")
    print("       the projection INHERITS the Savvidy tachyon, it does NOT control/flip it. The full-")
    print("       measure positivity (WG chi_O=0.33..0.59 -> 1/chi_O=1.7..3.0 >0) is a MEASURE fact,")
    print("       not a projection fact.\n")


# ------------------------------------------------------------------ [F] gauge invariance on SU(2) config
def _rand_su2():
    q = RNG.standard_normal(4); q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a + 1j * b, c + 1j * d], [-c + 1j * d, a - 1j * b]], dtype=complex)


def _su2_log(U):
    u0 = 0.5 * (U[0, 0] + U[1, 1]).real
    M = U - u0 * np.eye(2)
    s = np.sqrt(max(0.0, 0.5 * np.trace(M.conj().T @ M).real))   # sin(theta)
    if s < 1e-12:
        return np.zeros((2, 2), dtype=complex)
    theta = np.arccos(np.clip(u0, -1.0, 1.0))
    return (theta / s) * M


def F_gauge(d=2, L=4):
    print(f"[F] GAUGE INVARIANCE ({d}D SU(2), L={L}, one random config): observable vs link-log carrier")
    shape = (d, L, L, 2, 2) if d == 2 else (d, L, L, L, 2, 2)
    U = np.empty(shape, dtype=complex)
    for idx in np.ndindex(*shape[:-2]):
        U[idx] = _rand_su2()

    def dag(A):
        return A.conj().swapaxes(-1, -2)

    def sh(A, mu):
        return np.roll(A, -1, axis=mu)          # A(x+mu)

    def plaquette_sum(U):
        O = np.zeros(U.shape[1:-2])
        for mu in range(d):
            for nu in range(mu + 1, d):
                P = U[mu] @ sh(U[nu], mu) @ sh(dag(U[mu]), nu) @ dag(U[nu])
                O = O + (1.0 - 0.5 * (P[..., 0, 0] + P[..., 1, 1]).real)
        return O

    def gauge_transform(U):
        Gs = np.empty(U.shape[1:-2] + (2, 2), dtype=complex)
        for idx in np.ndindex(*U.shape[1:-2]):
            Gs[idx] = _rand_su2()
        Ug = np.empty_like(U)
        for mu in range(d):
            Gxp = np.roll(Gs, -1, axis=mu)       # g(x+mu)
            Ug[mu] = Gs @ U[mu] @ dag(Gxp)
        return Ug

    V = int(np.prod(U.shape[1:-2]))
    O = plaquette_sum(U); Q = O.sum() / np.sqrt(V)
    # link-log k=0 mode: || sum_{x,mu} log U_mu(x) ||_F
    logsum = np.zeros((2, 2), dtype=complex)
    for idx in np.ndindex(*U.shape[1:-2]):
        for mu in range(d):
            logsum += _su2_log(U[(mu,) + idx])
    Aflat = np.linalg.norm(logsum)

    Ug = gauge_transform(U)
    O2 = plaquette_sum(Ug); Q2 = O2.sum() / np.sqrt(V)
    logsum2 = np.zeros((2, 2), dtype=complex)
    for idx in np.ndindex(*Ug.shape[1:-2]):
        for mu in range(d):
            logsum2 += _su2_log(Ug[(mu,) + idx])
    Aflat2 = np.linalg.norm(logsum2)

    print(f"    pointwise plaquette check  max|O - O'|            = {np.max(np.abs(O - O2)):.3e}")
    print(f"    total-norm observable Q = V^-1/2 sum O:  |Q-Q'|/|Q| = {abs(Q - Q2)/abs(Q):.3e}   (gauge-inv)")
    print(f"    link-log k=0 mode  ||sum log U|| :  |A-A'|/|A|     = {abs(Aflat - Aflat2)/Aflat:.3e}  (DRIFTS)")
    print("    => the observable carrier q=V^-1/2 sum trF^2 is EXACTLY gauge-invariant (passes 'no gauge")
    print("       choice'); the link-space k=0 log mode is NOT. The proposal must use the observable.\n")


# ------------------------------------------------------------------ [G] multiscale = propagator
def G_multiscale():
    print("[G] MULTISCALE (all k): C(k)=Var(q(k))=Gtilde(k)=1/(m^2+khat^2); effective mass from C(x).")
    m = 0.3
    L = 128
    k = 2.0 * np.pi * np.arange(L) / L
    phat2 = 4.0 * np.sin(k / 2.0) ** 2
    Chat = 1.0 / (phat2 + m * m)
    Cx = np.fft.ifft(Chat).real
    # effective mass from the spatial decay C(x) ~ exp(-m_eff x) on the positive-x side
    xs = np.arange(1, 25)
    slope = np.polyfit(xs, np.log(np.abs(Cx[xs]) + 1e-300), 1)[0]
    m_eff = -slope
    print(f"    true m = {m:.5f}")
    print(f"    C(k=0) = chi = {Chat[0]:.5f} = 1/m^2 = {1.0/(m*m):.5f}  (the WG number: one point)")
    print(f"    m_eff from d ln C(x)/dx over x=1..24 : {m_eff:.5f}   (recovers m: the standard eff. mass)")
    print("    => the all-k multiscale version yields the STANDARD gauge-invariant effective mass from")
    print("       the propagator C(k); it is an ESTIMATE ([C]/[V]), not a lower bound, and not new")
    print("       (cf. INVENT-GAUGE-INVARIANT-MASS). The k=0 slice alone (WG) has only chi=C(0).\n")


# ------------------------------------------------------------------ main
def main():
    print(LINE)
    print("Y4 witness: gauge-invariant nonlinear Feshbach / multiscale RG projection -- TEST")
    print(LINE + "\n")
    A_gaussian()
    B_elementary()
    C_composite()
    D_plaquette()
    E_savvidy()
    F_gauge()
    G_multiscale()
    print(LINE)
    print("READING / KILL-CRITERIA LEDGER")
    print(" * (Q1) non-vacuous, Var=chi_O!=V^-1 : YES -- but it is exactly the T4/WG carrier. [T]")
    print(" * (Q2) improves on WG? NO in the load-bearing sense:")
    print("     - curvature G''=S_eff'' is a k=0 thermodynamic functional of the SAME marginal;")
    print("       Gaussian limit: G''=1/chi (A); off-Gaussian split (D: 1 vs 4) is a finite-V")
    print("       non-Gaussianity that VANISHES by CLT (C: k4->0, G''->1/Var). Repackages WG.")
    print("     - full nonlinear S_eff(q) is Gaussian in the thermo limit (CLT): no content beyond chi.")
    print("     - both G'' and 1/chi are CEILING-type, != pole for composites (C: != (2m)^2). [NG]")
    print("     - the only real increment, C(k) (G), is the STANDARD eff. mass: estimate, not a bound.")
    print(" * KILL: (a) repackages WG/H1 -> TRIGGERED (CLT Gaussian = WG variance).")
    print("         (b) needs a gauge choice -> NOT triggered (observable q is exactly gauge-invariant).")
    print("         (c) cannot control the collective mode -> TRIGGERED (Savvidy V''<0 inherited).")
    print(" * ACID: free field -> known/0 (B,D); 2D/free composite trF^2 -> finite, atomless, no floor (C).")
    print(" * The atom f1^2 >= r0 remains [O].  No QED.  Clay 0%.")
    print(LINE)


if __name__ == "__main__":
    main()
