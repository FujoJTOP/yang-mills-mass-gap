"""
wi_realtime.py -- WI flood: recover the IR-analytic atom signature

    A = pi / |Im tau_zero|   (nearest complex-tau zero of Z(tau))

and its reweighting-spread from REAL-tau multi-exponential data with noise.

Z(tau) = sum_j c_j exp(-E_j tau),  tau real in [0, beta].
  atom  : w0 (c0 at E0=0) + one atom f1^2 at m1  -> p=2 modes
  cut   : w0 + many exponentials approximating a continuum -> p large

Estimator: matrix-pencil / ESPRIT (Prony) on a Hankel matrix  ->  modes (c,E)
           -> Newton zeros of the recovered exponential sum in complex tau
           -> nearest |Im tau|  ->  A.
Weight-invariance: reweight the operator (change c_j but keep E_j); A is stable
for an atom, moves for a cut.

Run:  py -3.12 -I research/unconditional/w7/wi_realtime.py
"""
import numpy as np

rng = np.random.default_rng(7)


def gen(c, E, tau):
    c = np.asarray(c, complex)
    E = np.asarray(E, complex)
    return np.array([np.sum(c * np.exp(-E * t)) for t in tau])


def prony(tau, y, p, dt):
    """ESPRIT / matrix-pencil: recover (c, E) from real-tau samples y=Z(tau)."""
    N = len(tau)
    L = N // 2
    H = np.array([[y[i + j] for j in range(N - L)] for i in range(L)], dtype=complex)
    U, S, _ = np.linalg.svd(H, full_matrices=False)
    r = min(p, U.shape[1])            # model order actually used
    Us = U[:, :r]
    U1, U2 = Us[:-1, :], Us[1:, :]
    Psi, *_ = np.linalg.lstsq(U1, U2, rcond=None)
    z = np.linalg.eigvals(Psi)
    z = z[np.abs(z) > 1e-9]
    n = np.arange(N)
    V = z[None, :] ** n[:, None]
    c, *_ = np.linalg.lstsq(V, y, rcond=None)
    E = -np.log(z) / dt
    return c, E


def _f(c, E, tau):
    ex = np.exp(-E * tau)
    return np.sum(c * ex), np.sum(-E * c * ex)


def nearest_A(c, E, re_lo=-6.0, re_hi=6.0, im_hi=8.0, n=24, iters=80, ftol=1e-7):
    """Nearest complex-tau zero with Im>0; return A=pi/Im (or nan).

    Snap numerically-nonzero 'constant' modes (|E|<1e-9) to exactly 0 so the
    function does not underflow to a spurious zero at Re tau -> +inf; reject
    any zero found with |Re tau|>10 (underflow artifacts).
    """
    c = np.asarray(c, complex)
    E = np.where(np.abs(np.asarray(E, complex)) < 1e-9, 0.0, np.asarray(E, complex))
    scale = max(1.0, float(np.sum(np.abs(c))))
    best = np.inf
    for re0 in np.linspace(re_lo, re_hi, n):
        for im0 in np.linspace(0.05, im_hi, n):
            tau = complex(re0, im0)
            ok = True
            for _ in range(iters):
                f, df = _f(c, E, tau)
                if df == 0 or not np.isfinite(df):
                    ok = False
                    break
                step = f / df
                tau = tau - step
                if not np.isfinite(tau) or abs(tau.real) > 30 or abs(tau.imag) > 60:
                    ok = False
                    break
                if abs(step) < 1e-13:
                    break
            if not ok:
                continue
            f, _ = _f(c, E, tau)
            if abs(f) < ftol * scale and abs(tau.real) <= 10 and 1e-3 < tau.imag < best:
                best = tau.imag
    return np.pi / best if np.isfinite(best) else np.nan


def recover_A(y, tau, p, dt):
    c, E = prony(tau, y, p, dt)
    return nearest_A(c, E)


def sample(c, E, beta=12.0, N=120, noise=0.0):
    dt = beta / (N - 1)
    tau = np.arange(N) * dt
    y = gen(c, E, tau)
    y = y + noise * np.abs(y).mean() * rng.standard_normal(N)   # relative noise
    return y, tau, dt


def true_A(c, E):
    return nearest_A(np.asarray(c, complex), np.asarray(E, complex))


def main():
    print("=" * 72)
    print("TRUE values (noise-free, full model):")
    atom_c = [1.0, 0.5]; atom_E = [0.0, 1.0]
    print(f"  ATOM  w0=1 f1^2=0.5 m1=1  ->  A_true = {true_A(atom_c, atom_E):.5f}  (expect 1.0 = m1)")
    # cut: w0 + continuum on [1,4] approximated by p exp's
    cuts = {}
    for p in (4, 8, 16):
        E = np.linspace(1.0, 4.0, p)
        # weights ~ density * dE  (uniform density), positive
        w = np.full(p, 2.0) / p
        cuts[p] = (list(w), list(E))
        print(f"  CUT   p={p:2d} w0=1 continuum[1,4]  ->  A_true = {true_A([1.0]+list(w), [0.0]+list(E)):.5f}")

    print("\n" + "=" * 72)
    print("(A) ATOM recovery vs noise  (expected A=1.00000, spread~0)")
    for noise in (0.0, 1e-3, 1e-2, 5e-2):
        vals = []
        for _ in range(6):
            y, tau, dt = sample(atom_c, atom_E, noise=noise)
            vals.append(recover_A(y, tau, p=2, dt=dt))
        v = np.array([x for x in vals if np.isfinite(x)])
        if len(v):
            print(f"  noise={noise:5.3f}:  A = {v.mean():.4f} +- {v.std():.4f}   err={abs(v.mean()-1):.4f}   spread={v.max()-v.min():.4f}")
        else:
            print(f"  noise={noise:5.3f}:  FAILED (no finite A)")

    print("\n" + "=" * 72)
    print("(B) CUT recovery vs #exponentials p and noise  (expected: degrades, spread>0)")
    for p, (w, E) in cuts.items():
        for noise in (0.0, 1e-2):
            vals = []
            pmax = min(p + 1, 40)
            for _ in range(6):
                y, tau, dt = sample([1.0] + list(w), [0.0] + list(E), noise=noise)
                vals.append(recover_A(y, tau, p=pmax, dt=dt))
            v = np.array([x for x in vals if np.isfinite(x)])
            if len(v):
                print(f"  p={p:2d} noise={noise:5.3f}:  A = {v.mean():.4f} +- {v.std():.4f}   spread={v.max()-v.min():.4f}")
            else:
                print(f"  p={p:2d} noise={noise:5.3f}:  FAILED")

    print("\n" + "=" * 72)
    print("(C) WEIGHT-INVARIANCE signature (reweight operator, keep E):")
    # atom: reweight f1^2 (w0 fixed); A must stay = m1
    print("  ATOM reweight f1^2 in {0.2,0.5,1.5,5}:")
    a_vals = []
    for f2 in (0.2, 0.5, 1.5, 5.0):
        y, tau, dt = sample([1.0, f2], atom_E, noise=1e-3)
        a_vals.append(recover_A(y, tau, p=2, dt=dt))
    a_vals = np.array([x for x in a_vals if np.isfinite(x)])
    print(f"    A = {np.round(a_vals,4)}   SPREAD = {a_vals.max()-a_vals.min():.4f}  (atom: ~0)")
    # cut: reweight each f_i^2 randomly; A moves
    print("  CUT(p=8) reweight each mode's weight by random factor:")
    w = np.array(cuts[8][0]); E = cuts[8][1]
    c_vals = []
    for _ in range(4):
        mult = rng.uniform(0.3, 3.0, size=len(w))
        y, tau, dt = sample([1.0] + list(w * mult), [0.0] + list(E), noise=1e-3)
        c_vals.append(recover_A(y, tau, p=9, dt=dt))
    c_vals = np.array([x for x in c_vals if np.isfinite(x)])
    if len(c_vals):
        print(f"    A = {np.round(c_vals,4)}   SPREAD = {c_vals.max()-c_vals.min():.4f}  (cut: moves)")

    print("\n" + "=" * 72)
    print("(D) CONDITIONING: Hankel singular-value spread vs model order (noise-free)")
    print("    (why raw A is unreliable for a dense continuum: recovering many")
    print("     closely-spaced exponentials from finite real-tau data is ill-conditioned)")
    for tag, (w, E) in [("atom", ([0.5], [1.0]))] + [(f"cut{p}", cuts[p]) for p in cuts]:
        y, tau, dt = sample([1.0] + list(w), [0.0] + list(E), noise=0.0)
        N = len(tau); L = N // 2
        H = np.array([[y[i + j] for j in range(N - L)] for i in range(L)], dtype=complex)
        S = np.linalg.svd(H, compute_uv=False)
        keep = min(len(w) + 1, len(S))
        cond = S[0] / S[keep - 1]
        print(f"    {tag:6s}: sv_ratio S0/S_(p-1) = {cond:.2e}  (large => ill-conditioned)")

    print("\n" + "=" * 72)
    print("SUMMARY: A(atom) stable & weights don't move it; A(cut) raw value is unreliable")
    print("(ill-conditioned) BUT the reweighting-SPREAD is unambiguously larger than the atom's.")
    print("=> the robust finite-L atom signature is the SPREAD (A-invariance), not the raw A value.")


if __name__ == "__main__":
    main()
