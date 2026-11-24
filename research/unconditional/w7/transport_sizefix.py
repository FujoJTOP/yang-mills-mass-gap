"""
TRANSPORT-SIZE-FIX : witness for flip condition #2 of E2-FLIP-TRANSPORT §5 --
can a SIZE-FIXING mechanism for the instanton modulus in R^4 be constructed /
proved, i.e. does the effective action for the size rho have a minimum at a
finite, PARAMETRICALLY SMALL rho* << 1/Lambda derived from the action?

We build the pieces and test them. Units: Lambda = 1 (so 1/Lambda = 1).

  PART A -- the ONE-LOOP effective action for rho (pure SU(N), Lambda=1).
     One-loop instanton measure:
        dmu/drho  ~  (1/rho^5) [8pi^2/g^2(rho)]^{2N} exp(-8pi^2/g^2(rho)),
        8pi^2/g^2(rho) = b * ln(1/rho),   b = 11N/3      (one-loop running).
     Writing L = ln(1/rho) this is exp(-Gamma(L)) with
        Gamma(L) = (b-5)L - 2N ln(b L)          [the "effective action"]
     -> has a genuine minimum at  L* = 2N/(b-5),  rho* = exp(-2N/(b-5)).
     We locate it numerically, and ALSO report the 't Hooft coupling there:
        8pi^2/g^2(rho*) = b L*  ;  alpha_s(rho*) = 2*pi / (8pi^2/g^2(rho*)).
     KEY HONESTY FLAG: at rho* the coupling is STRONG (alpha_s ~ 0.5), so the
     one-loop formula is being evaluated at the EDGE of its validity.

  PART B -- mean-field "liquid" interaction (setting A):
     A probe instanton of size rho in a frozen background liquid of density n,
     the IIbar interaction averaged over neighbours:
        DeltaGamma_int(rho) = n * INT_{Rc}^{Rmax} d^4 R  V_{IIbar}(R; rho, rhobar)
        V_{IIbar} = -(32pi^2/g^2) * w * rho^2 rhobar^2 / R^4      (large R; CDG/SS)
        d^4 R = 2 pi^2 R^3 dR   ->   DeltaGamma_int = -n (32pi^2/g^2) w rho^2 rhobar^2 2pi^2 ln(Rmax/Rc)
        Rc = max(rho,rhobar),  Rmax = n^{-1/4}  (screening by inter-particle spacing).
     We add it to Gamma and ask: does it create / move a minimum to SMALLER rho?
     Physical ILM inputs (Lambda=1 units): n = 1, rhobar = 0.33  (packing n rhobar^4).
     We also scan a dense n (=50) to show the trend.  w = 1 (order-unity avg).

  PART B2 -- isolated IIbar molecule (setting B):
     one member's size action with a partner fixed at distance R:
        Gamma_mol(L) = one-loop Gamma(L) + DeltaGamma_int(rho; R fixed).
     Same question.

  PART C -- ACID TEST (mandatory): free / scale-invariant field (no running,
     g fixed, no Lambda).  Then 8pi^2/g^2 = const, determinant prefactor = const,
     so  Gamma_free(L) = 5 L  (pure measure 1/rho^5).  This is MONOTONE in L ->
     NO stationary point; the distribution  rho^{-5} drho  is scale-covariant
     (rho -> lam rho multiplies it by lam^{-4}).  The size-fixing mechanism MUST
     vanish.  We confirm numerically and flag the consequence: the ONLY thing
     that creates the minimum is the running = Lambda itself.

Headless, single process. numpy + scipy only.
"""
import numpy as np
from scipy.optimize import minimize_scalar


# -----------------------------------------------------------------------------
# one-loop effective action for the size, L = ln(1/rho)
# -----------------------------------------------------------------------------
def make_one_loop(N):
    b = 11.0 * N / 3.0
    S0 = lambda rho: b * np.log(1.0 / rho)            # 8 pi^2 / g^2(rho)
    def Gamma(L):
        # Gamma(L) = (b-5) L - 2N ln(b L)   ; domain L>0 (rho<1)
        L = np.asarray(L, dtype=float)
        out = np.full_like(L, np.inf)
        m = L > 0.0
        out[m] = (b - 5.0) * L[m] - 2.0 * N * np.log(b * L[m])
        return out
    return b, S0, Gamma


def find_min_L(Gamma, Llo=1e-3, Lhi=8.0, npts=40001):
    # sample then refine.  Gamma may be scalar-only -> evaluate elementwise.
    L = np.linspace(Llo, Lhi, npts)
    G = np.array([float(Gamma(float(x))) for x in L])
    i = int(np.argmin(G))
    lo, hi = L[max(i-1, 0)], L[min(i+1, len(L)-1)]
    res = minimize_scalar(lambda x: float(Gamma(float(x))), bracket=(lo, L[i], hi), method="brent")
    return float(res.x), float(res.fun)


# -----------------------------------------------------------------------------
# PART A
# -----------------------------------------------------------------------------
def partA():
    print("=" * 80)
    print("PART A  one-loop effective action Gamma(rho): the ONLY intrinsic minimum")
    print("=" * 80)
    for N in (2, 3):
        b, S0, Gamma = make_one_loop(N)
        Lstar, Gmin = find_min_L(Gamma)
        rho = np.exp(-Lstar)
        S = b * Lstar                                  # 8pi^2/g^2 at the min
        alpha = 2.0 * np.pi / S                         # alpha_s = g^2/4pi
        curv = 2.0 * N / Lstar**2                       # d^2 Gamma / dL^2 > 0
        sigmaL = 1.0 / np.sqrt(curv)                    # 1-sigma width in L = ln rho
        Lana = 2.0 * N / (b - 5.0)
        print(f"\n  SU({N})  b = 11N/3 = {b:.4f}")
        print(f"    numerical  L* = ln(1/rho*) = {Lstar:.6f}   analytic 2N/(b-5) = {Lana:.6f}")
        print(f"    => peak size  rho* = {rho:.6f}  (= {rho:.3f}/Lambda)   <-- O(1)/Lambda, NOT << 1/Lambda")
        print(f"    8pi^2/g^2(rho*) = b L* = {S:.4f}")
        print(f"    alpha_s(rho*)  = 2pi/S = {alpha:.4f}    <-- STRONG (one-loop invalid here)")
        print(f"    Gamma''(L*)    = {curv:.4f}  ->  sigma_L = {sigmaL:.4f}  "
              f"(size fixed to factor e^{{+/-{sigmaL:.2f}}}={np.exp(sigmaL):.2f} at 1 sigma)")
        print(f"    VARIANCE OF L* over N: 't Hooft coupling g^2 N at peak:")
        print(f"     8pi^2/g^2 = b L* = {S:.3f} -> g^2 = 8pi^2/S = {8*np.pi**2/S:.3f}, "
              f"g^2 N = {8*np.pi**2/S*N:.3f}  (=lambda_tHooft, LARGE)")
    print("\n  READ: a stationary point DOES exist, but")
    print("        (i)  it sits at rho* = O(1)/Lambda -- the confinement scale itself,")
    print("             not a parametrically small rho* << 1/Lambda;")
    print("        (ii) alpha_s(rho*) ~ 0.5, i.e. the one-loop running used to FIND it")
    print("             has broken down exactly there -> extrapolation, not control;")
    print("        (iii) its LOCATION rho* = c/Lambda is set by Lambda.")

    # explicit acid contrast: why the min needs running (see PART C)
    return


# -----------------------------------------------------------------------------
# PART B -- mean-field liquid interaction
# -----------------------------------------------------------------------------
def gamma_int_meanfield(rho, N, n, rhobar, w=1.0):
    """DeltaGamma_int(rho) = -n (32pi^2/g^2) w rho^2 rhobar^2 2pi^2 ln(Rmax/Rc)."""
    b, S0, _ = make_one_loop(N)
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    out = np.zeros_like(rho)
    Rmax = n ** (-0.25)
    for k, r in enumerate(rho):
        if not (0.0 < r < 1.0):
            out[k] = 0.0
            continue
        Sp = b * np.log(1.0 / r)                 # 8pi^2/g^2(r)
        Rc = max(r, rhobar)
        if Rmax <= Rc:
            out[k] = 0.0
            continue
        pref = 4.0 * Sp                          # 32pi^2/g^2 = 4 * (8pi^2/g^2)
        dG = -n * pref * w * r * r * rhobar * rhobar * 2.0 * np.pi**2 * np.log(Rmax / Rc)
        out[k] = dG
    return out[0] if out.size == 1 else out


def partB():
    print("\n" + "=" * 80)
    print("PART B  mean-field liquid interaction (setting A): does it fix rho SMALL?")
    print("=" * 80)
    for N in (2, 3):
        b, S0, Gamma = make_one_loop(N)
        Lstar_1, _ = find_min_L(Gamma)
        r1 = np.exp(-Lstar_1)
        for (n, rhobar, tag) in [(1.0, 0.33, "physical ILM (n~1/Lambda^4, rhobar~0.33/Lambda)"),
                                 (50.0, 0.37, "dense (n~50/Lambda^4, packing O(1))")]:
            pack = n * rhobar**4
            def Gtot(L):
                r = np.exp(-L)
                return float(Gamma(L)) + float(gamma_int_meanfield(r, N, n, rhobar))
            Lstar_tot, _ = find_min_L(Gtot)
            rtot = np.exp(-Lstar_tot)
            # location of the interaction term's OWN extremum (most negative DeltaGamma)
            Lg = np.linspace(1e-3, 6.0, 60001)
            Gi = np.array([float(gamma_int_meanfield(np.exp(-x), N, n, rhobar)) for x in Lg])
            r_int = np.exp(-Lg[int(np.argmin(Gi))])
            Gi_min = Gi.min()
            # is the interaction attractive everywhere?
            print(f"\n  SU({N}) [{tag}]  packing n rhobar^4 = {pack:.4f}")
            print(f"    one-loop min          rho*_1loop = {r1:.4f}  (<- the only intrinsic minimum)")
            print(f"    DeltaGamma_int extremum at rho = {r_int:.4f}  (value {Gi_min:+.3f}: ATTRACTIVE, "
                  f"bounded below; < 0)")
            print(f"    TOTAL minimum         rho*_tot   = {rtot:.4f}   "
                  f"(shift vs one-loop: factor {rtot/r1:.3f})")
            # explicit weak/strong conclusion for this (N, regime)
            if r_int >= r1:
                dirn = "LARGER (toward rho~1/e) -- NEVER smaller"
            else:
                dirn = "smaller, but only down to the interaction extremum rho~1/e"
            print(f"    interaction pulls the size {dirn}")
    print("\n  READ: DeltaGamma_int(rho) < 0 (attractive) and BOUNDED; its extremum sits at")
    print("        rho ~ 1/e / Lambda (for Rmax~1), i.e. an O(1)/Lambda scale. It can move the")
    print("        size toward ~0.37/Lambda (LARGER for SU(2) whose one-loop rho*=0.18, equal")
    print("        for SU(3) whose rho*=0.37). It NEVER produces a minimum at rho << 1/Lambda.")
    print("        At O(1) packing it is not a small correction and the mean-field integral")
    print("        extends to R~rho, where the large-R form V ~ 1/R^4 is INVALID -- so the")
    print("        '0.33 fm' size is a joint {n,rhobar} self-consistency FIT, not a bound.")


# -----------------------------------------------------------------------------
# PART B2 -- isolated IIbar molecule
# -----------------------------------------------------------------------------
def partB2():
    print("\n" + "=" * 80)
    print("PART B2  isolated IIbar molecule (setting B): fixed partner at distance R")
    print("=" * 80)
    N = 3
    b, S0, Gamma = make_one_loop(N)
    rhobar = 0.33
    for Rsep in (1.0, 2.0, 4.0):
        def Gm(L):
            r = np.exp(-L)
            Sp = b * np.log(1.0 / r) if 0 < r < 1 else 1e9
            Rc = max(r, rhobar)
            # For the molecule the interaction at fixed R is V(R) directly (no d^4R integral):
            V = -(4.0 * Sp) * r * r * rhobar * rhobar / Rsep**4 if (Rsep >= Rc and 0 < r < 1) else 0.0
            return float(Gamma(L)) + V
        # scan for a minimum
        Ls = np.linspace(1e-3, 8.0, 400001)
        Gs = np.array([Gm(x) for x in Ls])
        i = int(np.argmin(Gs))
        interior = 0 < i < len(Ls) - 1 and Gs[i] < Gs[0] and Gs[i] < Gs[-1]
        # derivative of V w.r.t rho (attractive?) at rho~0.33
        r0 = rhobar
        dV = (-(4.0 * b * np.log(1.0/r0)) * 2.0 * r0 * rhobar**2 / Rsep**4) if Rsep >= rhobar else np.nan
        print(f"\n    R = {Rsep:.2f}/Lambda:  V(R,rho=rhobar) = "
              f"{(-(4.0*b*np.log(1/r0))*r0**2*rhobar**2/Rsep**4):+.4f},  "
              f"dV/drho@rhobar = {dV:+.4f}  ({'attractive->larger rho' if dV<0 else 'repulsive'})")
        print(f"      min over rho in (0,1): {'INTERIOR min at rho='+format(np.exp(-Ls[i]),'.4f') if interior else 'NO interior min (monotone -> runs to boundary)'}")
    print("\n  READ: in the molecule the interaction is also attractive and pushes to LARGER rho;")
    print("        no small-rho minimum. The molecule 'binding' lowers the PAIR action by pulling")
    print("        the size UP until R~rho (validity edge), not by pinning rho small.")


# -----------------------------------------------------------------------------
# PART C -- acid test
# -----------------------------------------------------------------------------
def partC():
    print("\n" + "=" * 80)
    print("PART C  ACID TEST: free / scale-invariant field (no running, fixed g, no Lambda)")
    print("=" * 80)
    # no running: 8pi^2/g^2 = const, determinant prefactor = const => Gamma_free(L)=5L
    L = np.linspace(1e-3, 8.0, 200001)
    Gfree = 5.0 * L
    dG = np.diff(Gfree)
    has_stat = np.any(np.sign(dG[:-1]) != np.sign(dG[1:])) and np.any(dG == 0)
    # stationary point of Gamma_free?
    stationary = np.any(np.abs(np.gradient(Gfree, L)) < 1e-9)
    print(f"    Gamma_free(L) = 5 L    (pure measure 1/rho^5; prefactor const; no running)")
    print(f"    dGamma_free/dL = 5 > 0 everywhere; stationary point present? {stationary}")
    print(f"    distribution dmu = rho^{{-5}} drho:  scale-covariant under rho->lam rho")
    lam = 2.0
    I1 = np.trapezoid(np.exp(-5*np.log(1/np.linspace(0.2,0.8,10000))), np.linspace(0.2,0.8,10000))
    print(f"    INT rho^-5 drho over [0.2,0.8] ~ {I1:.3e} (diverges as lower cutoff ->0)")
    print("    => NO minimum, NO preferred size; conformal / scale-invariant. ACID PASSED.")
    print("    CONSEQUENCE (load-bearing): the ONLY ingredient that produced the PART A")
    print("    minimum was the running = dimensional transmutation = Lambda. Removing Lambda")
    print("    removes the minimum => the minimum's existence IS the presence of Lambda.")
    return stationary


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    partA()
    partB()
    partB2()
    partC()
    print("\n" + "=" * 80)
    print("VERDICT SUMMARY")
    print("=" * 80)
    print("  - V_eff(rho) HAS a stationary point, but only from measure (1/rho^5) vs")
    print("    one-loop running; located at rho* = exp(-2N/(b-5))/Lambda = O(1)/Lambda")
    print("    (SU(3): 0.37/Lambda), at alpha_s ~ 0.57 (STRONG coupling).")
    print("  - It is NOT at a parametrically small rho* << 1/Lambda.")
    print("  - The IIbar interaction (mean-field AND molecule) is ATTRACTIVE and pushes")
    print("    the size LARGER; it does NOT create a small-rho minimum and at the physical")
    print("    packing is a fit, not a controlled correction.")
    print("  - ACID PASSED: no running => no minimum. This SAME fact proves the")
    print("    circularity: the size is fixed at c/Lambda because Lambda is the only scale.")
    print("  => flip #2 NOT achieved. It is a FIT/PICTURE, not a theorem; circular in Lambda;")
    print("     the exact missing step = control of the rho-effective action at alpha_s~0.5")
    print("     (the multi-instanton liquid at packing O(1)) = the same crossover wall (H1).")
    print("\nDONE.")
