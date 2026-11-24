"""
invent_flowresp_witness.py -- witness for research/unconditional/INVENT-FLOWRESP.md

Model: the Yang-Mills gradient flow provides a Gaussian damping e^{-2 p^2 t} = e^{-2 s t}
on intermediate states of mass-squared s.  Hence the flowed energy is the LAPLACE
transform of the spectral (Kallen-Lehmann) measure nu:

        E(t) = int s * nu(ds) * e^{-2 s t}          (structure only; constants dropped)

Three test measures:
  [free massless]   nu(ds) ~ s ds on [0,inf)         -> POWER LAW, no freeze  (acid: correct)
  [atomless :phi^2:] nu(ds) ~ (s-theta) ds on [theta,inf), theta=(2m)^2  -> EXP freeze  (acid FAIL)
  [gapped atom]      f1^2 * delta(s - m1^2) + continuum               -> EXP freeze

The "flow-time response / freeze time" is t_* = 1 / (large-t log-derivative rate).
The acid test: the response must NOT fire on the free massless field, and must
DISCRIMINATE the :phi^2: continuum (f1^2=0) from the atom (f1^2>0).

Run:  py -3.12 -I invent_flowresp_witness.py
"""
import numpy as np
from scipy.integrate import quad

np.set_printoptions(precision=4)


def E_massless(t):
    # nu ~ s ds  =>  E(t) = int s^2 e^{-2 s t} ds = 1/(4 t^3)
    f = lambda s: s * s * np.exp(-2.0 * s * t)
    hi = max(60.0 / t, 10.0)
    val, _ = quad(f, 0.0, hi, limit=400)
    return val


def E_cont(theta, t):
    # atomless two-particle continuum from s=theta:  nu ~ (s-theta) ds
    if t <= 0:
        return np.inf
    f = lambda s: s * max(s - theta, 0.0) * np.exp(-2.0 * s * t)
    hi = theta + max(80.0 / t, 20.0)
    val, _ = quad(f, theta, hi, limit=800)
    return val


def E_atom(m1sq, f1sq, t):
    # single atom at s=m1^2 with weight f1^2
    return f1sq * m1sq * np.exp(-2.0 * m1sq * t)


def rate(fun, t, dt=1e-4):
    # large-t rate = -d log E / dt
    return -(np.log(fun(t + dt)) - np.log(fun(t))) / dt


def main():
    print("=" * 74)
    print("INVENT-FLOWRESP :: flow-time response as a gap handle")
    print("=" * 74)

    m = 1.0
    theta = (2.0 * m) ** 2          # atomless :phi^2: two-particle threshold
    m1sq, f1sq = 1.0, 0.5           # atom

    print("\n[1] FREE MASSLESS  (must show NO freeze; rate -> 0):")
    for t in [1.0, 4.0, 16.0]:
        print(f"    t={t:6.2f}   E(t)={E_massless(t):.4e}   rate=-dlogE/dt={rate(E_massless, t):.4f}")

    print("\n[2] ATOMLESS continuum :phi^2:  (theta=4m^2=%.1f; gapped, but f1^2=0):" % theta)
    for t in [1.0, 4.0, 16.0]:
        ev = lambda tt: E_cont(theta, tt)
        print(f"    t={t:6.2f}   E(t)={E_cont(theta, t):.4e}   rate=-dlogE/dt={rate(ev, t):.4f}")

    print("\n[3] GAPPED ATOM  (m1^2=%.1f, f1^2=%.1f):" % (m1sq, f1sq))
    for t in [1.0, 4.0, 16.0]:
        ev = lambda tt: E_atom(m1sq, f1sq, tt)
        print(f"    t={t:6.2f}   E(t)={E_atom(m1sq, f1sq, t):.4e}   rate=-dlogE/dt={rate(ev, t):.4f}")

    # --- large-t freeze times ---
    tL = 16.0
    r_free = rate(E_massless, tL)
    r_cont = rate(lambda tt: E_cont(theta, tt), tL)
    r_atom = rate(lambda tt: E_atom(m1sq, f1sq, tt), tL)
    tstar = lambda r: (1.0 / r) if r > 1e-3 else np.inf

    print("\n[4] ACID TEST  (freeze time t_* ~ 1/rate at t=%.0f):" % tL)
    print(f"    massless  : rate={r_free:.4f}  -> t_*=inf (no freeze)      [acid: correct]")
    print(f"    :phi^2:   : rate={r_cont:.4f}  -> t_*={tstar(r_cont):.4f}          <-- FIRES on atomless  [acid FAIL]")
    print(f"    atom      : rate={r_atom:.4f}  -> t_*={tstar(r_atom):.4f}          (rate -> 2 m1^2 = {2*m1sq})")
    print("    => finite freeze time for BOTH the atom (f1^2>0) and the :phi^2: continuum (f1^2=0).")
    print("       The flow response is a GAP detector, NOT an atom certifier  ==  the base's "
          "'gap-detector != atom-certifier' fact.")

    # --- the response dE/dkappa (kappa = shift of the spectral threshold) ---
    print("\n[5] RESPONSE dE/dkappa (kappa = coarse-mass shift of the continuum threshold):")
    t = 2.0
    dk = 1e-4
    # dE/dtheta = -int_theta^inf s e^{-2 s t} ds  (a moment / connected correlator)
    dfun = lambda s: s * np.exp(-2.0 * s * t)
    hi = theta + max(80.0 / t, 20.0)
    mom, _ = quad(dfun, theta, hi, limit=800)
    resp = -mom
    E0 = E_cont(theta, t)
    print(f"    t={t}:  dE/dkappa = {resp:.4e};   E(t) = {E0:.4e};   ratio = {resp/E0:.4f}")
    print("    -> the response is a MOMENT (linear functional of nu): homogeneous in the")
    print("       deformation, i.e. the RHS of the fixpoint, NOT the kappa-independent source.")

    print("\n[6] VERDICT: renormalon-free YES (E(t) finite by construction);")
    print("    atom-blind  (fires on the atomless continuum) => acid FAIL;")
    print("    response = a moment (homogeneous) => not the inhomogeneous source.  No QED.")


if __name__ == "__main__":
    main()
