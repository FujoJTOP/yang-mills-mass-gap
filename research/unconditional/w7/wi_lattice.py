"""
wi_lattice.py -- empirical calibration of the WI (weight-invariance) atom signature
on solvable finite systems.  LIGHTWEIGHT: exact linear algebra, no MC.

The WI signature (from CONSTRUCT-IR-ANALYTIC-NUM, PNT-IR-ANALYTIC):
  For Z_O(tau) = Sum_n |<n|O|Omega>|^2 exp(-(E_n-E_0) tau), the nearest complex-tau zero
  gives A_O = pi / Im(tau_zero).  CLAIM to test:
    * ISOLATED ATOM (a single dominant isolated level)  -> A_O is WEIGHT-INDEPENDENT
       (zero spread across operators O of the same quantum numbers)  [T]
    * CUT / CONTINUUM (dense / multi-level, no isolated atom)       -> A_O MOVES
       with the operator (nonzero spread)                          [V-num]

Sections:
  A. exact exponential sums: 2-term atom vs dense cut  (rigorous control)
  B. 1D transverse-field Ising chain, exact diagonalization, N=10 (real interacting model):
       even 0++ : O1 = sum sigma^z sigma^z (bond energy), O2 = sum sigma^x  -> CUT
       odd  ctrl: O  = sum sigma^z (single kink)                            -> ATOM

Run:  py -3.12 -I research/unconditional/w7/wi_lattice.py
"""
import numpy as np

# ----------------------------------------------------------------------------
# root finder for Z(tau)=sum_i w_i exp(-E_i tau), all E_i>=0, w_i>0 (real)
# ----------------------------------------------------------------------------
def _Z(tau, E, w):
    return np.sum(np.asarray(w) * np.exp(-np.asarray(E) * tau))

def _Zp(tau, E, w):
    return np.sum(-np.asarray(E) * np.asarray(w) * np.exp(-np.asarray(E) * tau))

def _newton(tau, E, w, iters=400):
    bad = complex(1e9, 1e9)
    for _ in range(iters):
        dz = _Zp(tau, E, w)
        z = _Z(tau, E, w)
        if not np.isfinite(z) or not np.isfinite(dz) or abs(dz) < 1e-300:
            return bad
        step = z / dz
        if not np.isfinite(step):
            return bad
        if abs(step) > 2.0:                      # damp large steps
            step = step / abs(step) * 2.0
        tau = tau - step
        if tau.real < -50 or tau.real > 50 or abs(tau.imag) > 90:
            return bad
        if abs(step) < 1e-13:
            break
    return tau

def nearest_zero_A(E, w, ntry=800, seed=0, re_rng=(-12, 12), im_rng=(0.03, 30)):
    """Return (A=pi/|Im tau_0|, tau_0) for the zero nearest the real axis; None if not found."""
    E = np.asarray(E, float); w = np.asarray(w, float)
    keep = w > 1e-13 * w.max()
    E, w = E[keep], w[keep]
    if len(E) == 1:
        return None, None
    rng = np.random.default_rng(seed)
    best = None
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for _ in range(ntry):
            t0 = complex(rng.uniform(*re_rng), rng.uniform(*im_rng))
            tau = _newton(t0, E, w)
            if not np.isfinite(tau.real):
                continue
            if abs(_Z(tau, E, w)) < 1e-7 * w.sum() and abs(tau.imag) > 1e-3:
                if best is None or abs(tau.imag) < abs(best.imag):
                    best = tau
    if best is None:
        return None, None
    return np.pi / abs(best.imag), best

def spread(A_list):
    A_list = [a for a in A_list if a is not None]
    if not A_list:
        return None
    return max(A_list) - min(A_list)

# ----------------------------------------------------------------------------
# Section A: exact exponential sums
# ----------------------------------------------------------------------------
def section_A():
    print("=" * 78)
    print("A.  EXACT EXPONENTIAL SUMS  (operational definition of the WI test)")
    print("=" * 78)
    m = 1.0
    print("\nA1. ATOM:  E=[0, m],  m=%.3f ; vary the weight ratio w0/w1" % m)
    E_atom = [0.0, m]
    print("    %-16s %-12s" % ("weights (w0,w1)", "A = pi/Im(tau0)"))
    A_atom = []
    for w0, w1 in [(1, 1), (1, 0.3), (1, 5), (0.2, 1), (7, 1), (1, 0.05)]:
        A, tau = nearest_zero_A(E_atom, [w0, w1])
        A_atom.append(A)
        print("    %-16s %-12.6f  (Re tau=%.3f)"
              % ("(%.2f,%.2f)" % (w0, w1), A, tau.real))
    print("    -> ATOM spread = %.2e   (expected ~0 : A=m=%.3f, weight-independent) [T]"
          % (spread(A_atom), m))
    print("    -> A_signed relative spread = %.2e"
          % (spread(A_atom) / np.mean(A_atom)))

    print("\nA2. CUT:  E=[0] + dense grid on [2m, 8m] (40 levels) ; vary the density")
    Ec = np.linspace(2 * m, 8 * m, 40)
    # three different spectral densities on the SAME energies (= three operators)
    def dens(kind):
        if kind == 0:
            return np.ones_like(Ec)               # flat
        if kind == 1:
            return np.exp(-(Ec - 2 * m) / 3.0)    # decaying
        return 1.0 + 0.5 * np.cos(3 * Ec)         # modulated
    print("    %-18s %-12s" % ("density on cut", "A = pi/Im(tau0)"))
    A_cut = []
    for k in range(3):
        w = np.concatenate([[1.0], dens(k)])
        E = np.concatenate([[0.0], Ec])
        A, tau = nearest_zero_A(E, w)
        A_cut.append(A)
        print("    %-18s %-12.6f" % (["flat", "decaying", "modulated"][k], A))
    print("    -> CUT spread = %.4f   (expected > 0 : A moves with the operator) [V-num]"
          % spread(A_cut))

    print("\nA3. NEAR-DEGENERATE pair (the atom is NOT isolated): E=[0, m, m+0.02m]")
    print("    the 'single level' m has a close neighbour -> does A move?")
    E3 = [0.0, m, m + 0.02]
    A3 = []
    for w2 in [0.01, 0.3, 1.0, 5.0]:
        A, tau = nearest_zero_A(E3, [1.0, 1.0, w2])
        A3.append(A)
        print("    weight on m+0.02 = %-6.2f  A = %.5f" % (w2, A))
    print("    -> NOT-isolated spread = %.4f   (moves: isolation is the KEY variable) [V-num]"
          % spread(A3))
    return spread(A_atom), spread(A_cut), spread(A3)

# ----------------------------------------------------------------------------
# Section B: 1D transverse-field Ising (TFIM), exact diagonalization, N=10
#   H = -J sum_i sz_i sz_{i+1} - h sum_i sx_i ,  periodic
#   Z2 symmetry P = prod sx ;  P-even (0++) : products of an even # of sz ; also sx
# ----------------------------------------------------------------------------
def _kron_list(mats):
    r = mats[0]
    for mm in mats[1:]:
        r = np.kron(r, mm)
    return r

def tfim_build(N, J, h):
    sx = np.array([[0., 1.], [1., 0.]])
    sz = np.array([[1., 0.], [0., -1.]])
    I2 = np.eye(2)
    dim = 2 ** N
    H = np.zeros((dim, dim))
    def op_at(op, i):
        return _kron_list([op if j == i else I2 for j in range(N)])
    for i in range(N):
        H -= J * (op_at(sz, i) @ op_at(sz, (i + 1) % N))
    for i in range(N):
        H -= h * op_at(sx, i)
    return H, sz, sx, I2

def sector_A_op(evals, evecs, Omat, E0, w0):
    """spectrum {E_n-E_0} with weights |<n|O|Omega>|^2 ; Omega = ground (index w0)."""
    v0 = evecs[:, w0]
    Ov = Omat @ v0
    c = evecs.T @ Ov                    # c_n = <n|O|Omega>
    w = c ** 2
    E = evals - E0
    keep = w > 1e-12 * w.max()
    return E[keep], w[keep]

def section_B():
    print("\n" + "=" * 78)
    print("B.  1D TRANSVERSE-FIELD ISING,  N=10, exact diagonalization")
    print("=" * 78)
    N = 10
    J, h = 0.5, 1.0           # DISORDERED phase: unique ground state, no symmetry-breaking doublet
    H, sz, sx, I2 = tfim_build(N, J, h)
    evals, evecs = np.linalg.eigh(H)
    E0 = evals[0]; w0 = 0
    print("    J=%.2f h=%.2f  E0=%.6f  E1-E0=%.6f  (disordered/gapped, unique ground state)"
          % (J, h, E0, evals[1] - E0))

    def op_at(op, i):
        return _kron_list([op if j == i else I2 for j in range(N)])
    # ---- even 0++ operators ----
    O_bond = sum(op_at(sz, i) @ op_at(sz, (i + 1) % N) for i in range(N))
    O_sx   = sum(op_at(sx, i) for i in range(N))
    # a second even operator "smeared" variant: bond with a different profile
    O_bond2 = sum((1.0 + 0.3 * ((-1) ** i)) * (op_at(sz, i) @ op_at(sz, (i + 1) % N))
                  for i in range(N))
    print("\n    EVEN sector (0++):  the two-kink CUT")
    for name, O in [("sum sz sz (bond)", O_bond), ("sum sx (field)", O_sx),
                    ("smeared bond", O_bond2)]:
        E, w = sector_A_op(evals, evecs, O, E0, w0)
        A, tau = nearest_zero_A(E, w)
        print("      %-20s A=%-10s  n_states=%d  lowest exc=%.4f"
              % (name, ("None" if A is None else "%.5f" % A), len(E), E.min() if len(E) else -1))
    # ---- odd control: single kink = ATOM ----
    O_sz = sum(op_at(sz, i) for i in range(N))
    O_sz2 = sum((1.0 + 0.25 * ((-1) ** i)) * op_at(sz, i) for i in range(N))
    print("\n    ODD sector (control, single kink):  the ATOM")
    A_odd = []
    for name, O in [("sum sz", O_sz), ("smeared sz", O_sz2)]:
        E, w = sector_A_op(evals, evecs, O, E0, w0)
        A, tau = nearest_zero_A(E, w)
        A_odd.append(A)
        print("      %-20s A=%-10s  n_states=%d  lowest exc=%.4f"
              % (name, ("None" if A is None else "%.5f" % A), len(E), E.min() if len(E) else -1))
    print("    -> ODD (atom) spread = %.4f ; lowest excitation (single-kink mass)=%.4f"
          % (spread(A_odd), evals[1] - E0))

    # even-sector spread
    even_As = []
    for O in [O_bond, O_sx, O_bond2]:
        E, w = sector_A_op(evals, evecs, O, E0, w0)
        A, _ = nearest_zero_A(E, w)
        even_As.append(A)
    print("    -> EVEN (0++, cut) spread = %.4f" % spread(even_As))
    return spread(even_As), spread(A_odd), evals[1] - E0

# ----------------------------------------------------------------------------
if __name__ == "__main__":
    sa_atom, sa_cut, sa_near = section_A()
    sb_even, sb_odd, m_kink = section_B()
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print("  exact-sum ISOLATED ATOM spread = %.2e  (-> 0)   [T]" % sa_atom)
    print("  exact-sum CUT          spread = %.4f  (-> >0)  [V-num]" % sa_cut)
    print("  exact-sum NOT-ISOLATED spread = %.4f  (moves)  [V-num]" % sa_near)
    print("  TFIM 0++ (even, CUT)   spread = %.4f   [V-num]" % sb_even)
    print("  TFIM odd (kink band)   spread = %.4f   lowest level = %.4f  [V-num]"
          % (sb_odd, m_kink))
    print("  VERDICT: zero spread <=> ISOLATED atom (gap to every other level).")
    print("           Not-isolated (band/near-degenerate/cut) -> nonzero spread.")
    print("           TFIM 0++ is a two-kink CUT (moves), as is tr F^2; a clean isolated")
    print("           0++ atom is absent in free-fermion TFIM (the YM bound-state problem).")
