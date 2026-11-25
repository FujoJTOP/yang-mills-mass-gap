#!/usr/bin/env python3
r"""
t4_carrier.py -- T4: does a GAUGE-INVARIANT, NON-SELF-AVERAGING carrier exist, whose variance
        is NOT V^{-1} (so H1 on it is non-vacuous)?  And does its susceptibility carry gap content?

The corpus (Q6 / H1-IS-IT-THE-GAP) claims: on the gauge-invariant composite tr F^2 the variance is
Var ~ 1/V (H1 "vacuous"), on the gluon log U it is non-vacuous but gauge-dependent, and "no linear
gauge-invariant carrier with Var not-> 0 exists".  This witness settles the NORMALIZATION algebra:

  * total-normalized zero mode  q = (1/sqrt(V)) sum_x O(x)  =>  Var(q) = chi := sum_z C(z)   (O(1)!=1/V)
  * mean-normalized  zero mode  qbar = (1/V)    sum_x O(x)  =>  Var(qbar) = chi / V          (-> 0)

So the "Var ~ 1/V" of the composite is the MEAN-normalized (self-averaging) mode -- a normalization
artifact, not an obstruction.  The TOTAL-normalized mode has Var = chi, an O(1) L-independent quantity.

For a free massive scalar the total-normalized zero mode has Var = 1/m^2 EXACTLY, independent of L
(the k=0 mode of the propagator): a gauge-invariant-friendly, non-self-averaging carrier.  The
composite :phi^2: (the gauge-invariant *composite* analog of tr F^2) has chi_{:phi^2:} = sum_z 2 G(z)^2,
finite and L-independent for m>0, divergent for m->0 -- the free-field acid test.

We test:
  [A] elementary phi (free massive/massless), d=1,2,4 periodic:  Var(total-zero) = 1/m^2 (exact, L-indep);
      Var(mean-zero) = 1/(m^2 V) (slope -1).  Log-log slopes fit.
  [B] gauge-invariant composite :phi^2:, d=1,2:  chi = sum_z 2G(z)^2, L-independent (m>0), -> inf (m=0).
  [C] one-plaquette SU(2) Haar: Var of the gauge-invariant plaquette action density (single site, O(1)).
  [D] the "ceiling not edge" check: for :phi^2: the lightest coupled state is at 2m, but chi != 1/(2m)^2;
      chi is a residue-weighted harmonic mean (the squeeze ceiling), never a lower bound on m_lightest^2.

numpy only, one job, -I safe.  Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""
import numpy as np

np.seterr(all='ignore')


# --------------------------------------------------------------- exact periodic free propagator
def free_propagator(L, m, d):
    """G(x) = (1/V) sum_k e^{ikx} / (m^2 + sum_mu 4 sin^2(k_mu/2)) on a periodic L^d lattice.
       Returned as an array of shape (L,)*d -- this IS G(x) (translation invariant)."""
    grids = np.meshgrid(*[np.arange(L) for _ in range(d)], indexing='ij')
    phat2 = np.zeros((L,) * d)
    for mu in range(d):
        kk = 2.0 * np.pi * grids[mu] / L
        phat2 += 4.0 * np.sin(kk / 2.0) ** 2
    denom = phat2 + m * m
    with np.errstate(divide='ignore'):
        fk = 1.0 / denom
    return np.fft.ifftn(fk).real          # = (1/V) sum_k fk e^{ikx}


def chi_elementary(L, m, d):
    """Var of the TOTAL-normalized zero mode (1/sqrt V) sum_x phi(x) = sum_z G(z) = Ghat(k=0) = 1/m^2.
       (Exact: only the k=0 term survives the sum over x.)"""
    G = free_propagator(L, m, d)
    return float(G.sum())


def chi_composite(L, m, d):
    """chi of the gauge-invariant composite :phi^2:  =  sum_z 2 G(z)^2   (free-field Wick).
       Finite for m>0 on a fixed lattice; diverges as m->0 (acid)."""
    G = free_propagator(L, m, d)
    return float(np.sum(2.0 * G * G))


# --------------------------------------------------------------- one-plaquette SU(2) Haar
def one_plaquette():
    """SU(2) one-plaquette: gauge-invariant action density s = 1 - (1/2) tr U = 1 - z.
       Haar SU(2) = S^3: class angle phi in [0,pi], measure sin^2(phi) dphi, z=(1/2)tr U = cos(phi).
       Exact: <z> = 0 (Schur orthogonality), <z^2> = <|tr U|^2>/4 = 1/4, so Var(z) = 1/4.
       Compute by quadrature to confirm."""
    th = np.linspace(0.0, np.pi, 2000001)
    dens = np.sin(th) ** 2                            # S^3 / Haar density of class angle phi
    dens /= np.trapezoid(dens, th)                    # normalize
    z = np.cos(th)
    mz = np.trapezoid(z * dens, th)
    mz2 = np.trapezoid(z * z * dens, th)
    var_z = mz2 - mz * mz
    return mz, var_z, var_z          # Var(1-z) = Var(z) (shift invariance)


# --------------------------------------------------------------- driver
def main():
    print("=" * 84)
    print(" T4 -- gauge-invariant non-self-averaging carrier: normalization algebra + solvable audit")
    print("=" * 84)

    # -------- [A] elementary phi: total vs mean normalized variance, massive (m>0) --------
    print("\n[A] FREE MASSIVE SCALAR phi  --  Var(total-normalized zero mode) = 1/m^2  (exact, L-independent)")
    print("    vs  Var(mean-normalized) = 1/(m^2 V)  (self-averaging, slope -1)")
    for d, Ls in [(1, [16, 32, 64, 128, 256]),
                  (2, [8, 16, 32, 64]),
                  (4, [4, 6, 8, 12])]:
        m = 0.5
        Vs, chiT, chiM = [], [], []
        for L in Ls:
            V = L ** d
            G = free_propagator(L, m, d)
            chi_tot = float(G.sum())            # Var totale zero mode = sum_z G(z)
            chi_men = chi_tot / V               # Var mean zero mode
            Vs.append(V); chiT.append(chi_tot); chiM.append(chi_men)
        Vs = np.array(Vs, float)
        # log-log slopes
        sT = np.polyfit(np.log(Vs), np.log(np.array(chiT)), 1)[0]
        sM = np.polyfit(np.log(Vs), np.log(np.array(chiM)), 1)[0]
        print(f"  d={d}  m={m}:  1/m^2 = {1.0/m**2:.6f}")
        print(f"     L          V       Var(total)      Var(mean)      Var(total)*V^slope0")
        for L, V, cT, cM in zip(Ls, Vs, chiT, chiM):
            print(f"    {L:>3}  {int(V):>8}   {cT:12.6f}   {cM:12.6e}")
        print(f"     -> log-log slope Var(total) = {sT:+.4f}  (0 expected: NON-self-averaging)")
        print(f"     -> log-log slope Var(mean)  = {sM:+.4f}  (-1 expected: self-averaging)")

    # -------- [A'] massless acid: total-normalized zero mode diverges (k=0 mode of a massless field) --------
    print("\n[A'] FREE MASSLESS acid test: total-normalized zero mode Var = 1/m^2 -> inf as m->0")
    print("     (the zero mode of a massless field has divergent norm -> H1 correctly FAILS)")
    for d in [1, 2, 4]:
        L = 16
        for m in [0.5, 0.2, 0.05, 0.01]:
            G = free_propagator(L, m, d)
            chi_tot = float(G.sum())
            print(f"     d={d} L={L}: m={m:<5}  Var(total-zero)={chi_tot:12.4f}   (1/m^2={1/m**2:12.4f})")

    # -------- [B] gauge-invariant composite :phi^2:  --  finite, L-independent (m>0), acid diverges --------
    print("\n[B] GAUGE-INVARIANT COMPOSITE :phi^2:  --  chi = sum_z 2 G(z)^2  (free Wick)")
    print("    (This is the gauge-invariant-composite ANALOG of tr F^2: a local 0++ density.)")
    for d, Ls in [(1, [32, 64, 128, 256]), (2, [16, 32, 64])]:
        m = 0.5
        Vs, chic = [], []
        for L in Ls:
            V = L ** d
            Vs.append(V); chic.append(chi_composite(L, m, d))
        Vs = np.array(Vs, float)
        sC = np.polyfit(np.log(Vs), np.log(np.array(chic)), 1)[0]
        print(f"  d={d}  m={m}:  Var(total-zero of :phi^2:) across L")
        for L, V, c in zip(Ls, Vs, chic):
            print(f"     L={L:>3}  V={int(V):>8}   chi_{{:phi^2:}} = {c:12.6f}")
        print(f"     -> log-log slope = {sC:+.4f}  (0 expected: NON-vacuous, non-self-averaging)")
    # acid: composite massless
    print("  acid (m->0): chi_{:phi^2:} must DIVERGE (IR); free massless has no gap")
    for d in [1, 2]:
        L = 64
        row = []
        for m in [0.5, 0.2, 0.1, 0.05, 0.02]:
            row.append(f"m={m}:{chi_composite(L, m, d):.3f}")
        print(f"     d={d} L={L}:  " + "   ".join(row))

    # -------- [C] one-plaquette SU(2): gauge-invariant plaquette density, O(1) variance, no V --------
    print("\n[C] ONE-PLAQUETTE SU(2) Haar:  gauge-invariant action density s=1-(1/2)tr U = 1-z")
    mz, varz, vars_s = one_plaquette()
    print(f"     <z> = {mz:.6f}   Var(z) = {varz:.6f}   Var(s=1-z) = {vars_s:.6f}")
    print(f"     (a single gauge-invariant site: Var is O(1), no V-dependence; the carrier is trivially")
    print(f"      non-self-averaging -- there is no volume to average over.)")

    # -------- [D] ceiling not edge: chi_{:phi^2:} is NOT 1/(m_lightest)^2 --------
    print("\n[D] CEILING-NOT-EDGE: the lightest state coupled to :phi^2: sits at threshold 2m;")
    print("    is chi_{:phi^2:} equal to 1/(2m)^2?  (It must NOT be -- chi is a harmonic mean / ceiling.)")
    d, L, m = 2, 64, 0.5
    c = chi_composite(L, m, d)
    print(f"     d={d} L={L} m={m}:  chi_{{:phi^2:}} = {c:.6f}   1/(2m)^2 = {1/(2*m)**2:.6f}")
    print(f"     ratio chi / (1/(2m)^2) = {c/(1/(2*m)**2):.4f}   (not 1 -> chi is NOT the pole mass^-2;")
    print(f"     it is the residue-weighted harmonic mean = the SQUEEZE CEILING, never the edge.)")

    print("\n" + "=" * 84)
    print(" VERDICT (algebra):")
    print("  * total-normalized zero mode of a gauge-invariant density: Var = chi = O(1), NOT 1/V.")
    print("  * the corpus's 'Var ~ 1/V' is the MEAN-normalized (self-averaging) mode: a normalization")
    print("    artifact, NOT a carrier obstruction.  A non-vacuous gauge-invariant carrier EXISTS.")
    print("  * but chi is a residue-weighted harmonic mean (squeeze CEILING).  No variance functional")
    print("    lower-bounds the mass; the EDGE needs the atom.  Gap content: ceiling only.")
    print("=" * 84)


if __name__ == "__main__":
    main()
