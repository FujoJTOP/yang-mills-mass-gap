"""
base_minimal_model.py -- witness for research/unconditional/BASE-MINIMAL-MODEL.md

DEMONSTRATION ONLY (small 3D Ising Metropolis, numpy, <2 s).
3D Z2 lattice gauge theory is dual to the 3D Ising model (Wegner 1971).
Under the duality dictionary:
    gauge (electric, local)  <->  Ising (magnetic, order parameter)
    Wilson loop (electric)   <->  dual domain-wall / high-T expansion
    't Hooft loop (magnetic) <->  Ising magnetization m           <-- the NON-LOCAL order parameter
    confined (gapped) phase  <->  Ising ORDERED phase (m != 0)
    gapless / deconfined     <->  Ising critical point (xi -> inf)

We measure, as a function of beta, BOTH:
    |m|  -- the (dual 't Hooft) ORDER PARAMETER  -> sharply turns on below beta_c
    E    -- the LOCAL energy density (plaquette analog) -> smooth through beta_c
to illustrate (NOT prove) that the discriminator is the non-local order parameter.
The rigorous "no local gauge-invariant order parameter" statement is Elitzur's
theorem [T]; the numerics here only illustrate, they do not establish it.

Run:  py -3.12 -I base_minimal_model.py
"""
import numpy as np

def nx(s, ax):
    return np.roll(s, -1, axis=ax)

def metropolis(L, beta, n_eq=600, n_meas=600, rng=None):
    s = rng.choice([-1, 1], size=(L, L, L)).astype(np.int8)
    idx = np.indices((L, L, L))
    parity = (idx[0] + idx[1] + idx[2]) & 1
    masks = [parity == 0, parity == 1]
    msum = 0.0; esum = 0.0; cnt = 0
    for sweep in range(n_eq + n_meas):
        for mk in masks:
            nb = (np.roll(s, 1, 0) + np.roll(s, -1, 0)
                  + np.roll(s, 1, 1) + np.roll(s, -1, 1)
                  + np.roll(s, 1, 2) + np.roll(s, -1, 2))
            dE = 2 * s * nb                       # energy change of flipping s
            r = rng.random(s.shape)
            flip = mk & ((dE <= 0) | (r < np.exp(-beta * dE)))
            s[flip] = -s[flip]
        if sweep >= n_eq:
            msum += abs(s.mean())
            # energy per spin: E = -(1/2) <s * n>  (each bond counted twice)
            nb = (np.roll(s, 1, 0) + np.roll(s, -1, 0)
                  + np.roll(s, 1, 1) + np.roll(s, -1, 1)
                  + np.roll(s, 1, 2) + np.roll(s, -1, 2))
            esum += -0.5 * (s * nb).mean()
            cnt += 1
    return msum / cnt, esum / cnt

def main():
    rng = np.random.default_rng(20261124)
    L = 8
    betas = [0.150, 0.190, 0.210, 0.2217, 0.235, 0.260, 0.300]
    print("=" * 74)
    print("3D Ising (= dual 3D Z2 gauge) -- local energy vs non-local order parameter")
    print("L=%d, beta_c(3D Ising, J=1) ~ 0.2217" % L)
    print("=" * 74)
    print("  beta     |m| (dual 't Hooft)      E (local / plaquette)   phase")
    print("  " + "-" * 70)
    for b in betas:
        m, e = metropolis(L, b, rng=rng)
        tag = "ORDERED (dual: confined/gapped)" if m > 0.5 else \
              ("critical" if 0.2 < m < 0.5 else "DISORDERED (dual: deconfined/gapless)")
        print("  %.4f   %.4f                   %+.4f                %s" % (b, m, e, tag))
    print()
    print("Reading: |m| (the NON-LOCAL 't Hooft/dual-magnetization order parameter)")
    print("turns ON below beta_c; the LOCAL energy E is smooth through beta_c.")
    print("In the gauge theory this is: the confinement order parameter is MAGNETIC")
    print("(non-local) -- local electric (plaquette) data cannot read it (Elitzur).")
    print()
    print("NOT a proof of anything about 4D SU(N). Illustration only. No QED.")

if __name__ == "__main__":
    main()
