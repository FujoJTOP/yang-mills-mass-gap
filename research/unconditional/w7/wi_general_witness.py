"""
wi_general_witness.py -- WI-GENERAL witness (fork). Lightweight, numpy only.
Run: py -3.12 -I research/unconditional/w7/wi_general_witness.py

Model. A finite system gives Z(tau) = w0 + sum_i a_i exp(-E_i tau),
  E_i = spectrum of H (operator-INDEPENDENT -> the SUPPORT of the spectral measure),
  a_i = |<i|O|0>|^2, w0 = |<0|O|0>|^2 (operator-DEPENDENT -> the WEIGHTS/residues).
"Reweighting" = change the residues {a_i}, w0 at FIXED poles {E_i} (mix the operator O -> O+tO').

We test the claim: a functional A(Z) that depends only on the SUPPORT (pole positions),
not the weights, is reweighting-invariant IFF the edge is an ISOLATED atom.

Tests:
  (T1) atom (1 isolated pole [+V EV]):  A = pi/|Im tau_zero| should equal E0, weight-independent.
  (T2) isolated atom + separate continuum (gap delta>0): A stays ~E0 (exp-small correction).
  (T3) continuum at the edge (dense poles at E0): A MOVES under reweighting.
  (T4) zero ladder: atom -> odd multiples pi/E0*(1,3,5,...); cut -> irregular.
  (T5) real-tau effective mass m_eff(tau) = -d/dtau ln C: atom flat (=E0), cut drifts (lattice-visible).
"""
import numpy as np
np.seterr(all='ignore')
rng = np.random.default_rng(0)


def Zsum(tau, poles, weights, w0):
    poles = np.asarray(poles, float); weights = np.asarray(weights, float)
    return w0 + np.sum(weights * np.exp(-poles * tau))


def Zprime(tau, poles, weights):
    poles = np.asarray(poles, float); weights = np.asarray(weights, float)
    return np.sum(-poles * weights * np.exp(-poles * tau))


def find_zeros(poles, weights, w0, re_lo=-5.0, re_hi=5.0, nre=44,
               ims=(0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 25.0)):
    roots = []
    def f(t): return Zsum(t, poles, weights, w0)
    def fp(t): return Zprime(t, poles, weights)
    for re in np.linspace(re_lo, re_hi, nre):
        for im in ims:
            for sgn in (1.0, -1.0):
                t = complex(re, sgn * im)
                for _ in range(80):
                    z = f(t); d = fp(t)
                    if abs(d) < 1e-300:
                        break
                    step = z / d
                    t = t - step
                    if abs(step) < 1e-13:
                        break
                if abs(f(t)) < 1e-9 and 1e-5 < abs(t.imag) < 60.0:
                    if not any(abs(t - r) < 1e-4 for r in roots):
                        roots.append(t)
    roots.sort(key=lambda t: abs(t.imag))
    return roots


def A_nearest(poles, weights, w0):
    roots = find_zeros(poles, weights, w0)
    if not roots:
        return None, None, 0
    r0 = roots[0]
    return np.pi / abs(r0.imag), roots, len(roots)


def reweight(weights, w0, scale=0.6):
    """Positive multiplicative noise on residues at FIXED poles = operator reweighting."""
    noise = np.exp(scale * rng.standard_normal(len(weights)))
    return np.asarray(weights) * noise, w0 * float(np.exp(scale * rng.standard_normal()))


def train():
    print("=" * 78)
    print("(T1) SINGLE ISOLATED ATOM  E0=1.0 :  A = pi/|Im tau_zero|  should equal E0")
    print("=" * 78)
    poles = [1.0]
    for w0, a in [(1.0, 0.5), (1.0, 0.2), (1.0, 2.0), (0.3, 5.0), (7.0, 0.05)]:
        A, roots, n = A_nearest(poles, [a], w0)
        print(f"  w0={w0:5.2f} a={a:6.3f}  ->  A={A:.8f}   (zeros found: {n}, |Im|min={abs(roots[0].imag):.6f})")
    print("  => A is E0 regardless of (w0,a): the atom mass rides on Im(tau), the weight on Re(tau).\n")

    print("=" * 78)
    print("(T2) ISOLATED ATOM (E0=1.0) + CONTINUUM (from 2.0 up, gap=1.0): nearest-zero is UNRELIABLE;")
    print("     the ROBUST WI signature is the effective-mass PLATEAU m_eff(tau) -> E0.")
    print("=" * 78)
    poles = np.concatenate(([1.0], 2.0 + 0.1 * np.arange(0, 20)))
    base = np.concatenate(([0.5], np.exp(-1.0 * np.arange(0, 20) * 0.1)))
    def meffw(tau, p, w):
        e = np.asarray(w, float) * np.exp(-np.asarray(p, float) * tau)
        return np.sum(np.asarray(p, float) * e) / np.sum(e)
    print("  m_eff(tau): " + "  ".join(f"({t})={meffw(t,poles,base):.4f}" for t in (0.3, 1, 2, 4, 8)))
    print("  -> plateau at E0=1.0 for large tau (atom dominates; continuum e^{-2tau} dies).")
    # nearest-Im zero is polluted by the continuum (pi/2 < pi/1): show it is NOT atom-valued
    vals = [A_nearest(poles, *reweight(base, 1.0))[0] for _ in range(5)]
    print(f"  nearest-zero A reweighted: {[round(v,4) for v in vals]} -> NOT E0 (grabs continuum zero).")
    # the plateau IS reweighting-invariant:
    plat = [meffw(8.0, poles, reweight(base, 1.0)[0]) for _ in range(6)]
    print(f"  plateau m_eff(8) reweighted: {[round(v,6) for v in plat]}  spread={max(plat)-min(plat):.1e}")
    print("  => WI atom signature = effective-mass plateau (operator-independent mass), NOT nearest-zero.\n")

    print("=" * 78)
    print("(T3) CONTINUUM AT THE EDGE (dense poles from E0=1.0, no isolation)")
    print("=" * 78)
    poles = 1.0 + 0.05 * np.arange(0, 60)
    wts = np.exp(-0.8 * np.arange(0, 60) * 0.05)
    vals = []
    for _ in range(6):
        w, w0 = reweight(wts, 1.0)
        A, _, n = A_nearest(poles, w, w0)
        vals.append(A)
    print(f"  nearest-edge A under reweighting: {[round(v,5) for v in vals]}")
    print(f"  spread = {max(vals)-min(vals):.3f}  => NOT invariant (continuum: density sets the zeros).\n")

    print("=" * 78)
    print("(T4) ZERO LADDER: atom = odd multiples pi/E0*(1,3,5,..); cut = irregular")
    print("=" * 78)
    Aa, roots_a, _ = A_nearest([1.0], [0.5], 1.0)
    ims_a = sorted(abs(r.imag) for r in roots_a)[:6]
    print(f"  ATOM |Im| ladder (E0=1): {[round(v,4) for v in ims_a]}")
    print(f"    ratios to pi/E0=pi: {[round(v/np.pi,3) for v in ims_a]}  (expect 1,3,5,...)")
    Ac, roots_c, _ = A_nearest(poles, wts, 1.0)
    ims_c = sorted(abs(r.imag) for r in roots_c)[:6]
    print(f"  CUT  |Im| ladder:        {[round(v,4) for v in ims_c]}")
    print(f"    ratios:                {[round(v/np.pi,3) for v in ims_c]}  (irregular)\n")

    print("=" * 78)
    print("(T5) REAL-tau EFFECTIVE MASS  m_eff(tau) = sum a_i E_i e^{-E_i tau} / sum a_i e^{-E_i tau}")
    print("=" * 78)
    def meff(tau, poles, wts):
        poles = np.asarray(poles, float); wts = np.asarray(wts, float)
        e = wts * np.exp(-poles * tau)
        return np.sum(poles * e) / np.sum(e)
    print("  ATOM E0=1.0:  " + "  ".join(f"m_eff({t})={meff(t,[1.0],[0.5]):.4f}" for t in (0.5, 1, 2, 4)))
    print("  CUT  (dense): " + "  ".join(f"m_eff({t})={meff(t,poles,wts):.4f}" for t in (0.5, 1, 2, 4)))
    print("  => atom: flat (operator-independent mass, WI); cut: drifts (density-contaminated).\n")

    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print("  Atom (isolated):  A invariant, ladder periodic, m_eff flat  -> WI atom signature.")
    print("  Continuum@edge:   A moves,   ladder irregular, m_eff drifts -> NOT WI.")
    print("  --- WI (operator-reweighting invariance) is an atom signature iff the edge is isolated;")
    print("      uniform isolation across (L,a) at beta=2.5 is H1. no QED.")


train()
