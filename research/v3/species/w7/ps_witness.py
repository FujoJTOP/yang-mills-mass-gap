"""ps_witness.py -- witness for PS-NONNORMAL (non-normality / pseudospectral species).

Deterministic, numpy(+scipy) only, one process, no downloads.  py -3.12

[A] 2x2 toy: identical spectrum {1,e^-m}, different normality.  kappa(z)=dist*||(T-z)^-1||,
    Kreiss index ||T||/rho(T), numerical-radius gap omega-rho, dep(T).
[B] Gelfand deficit  m >= -(1/(n a0)) log ||A^n||^{1/n} : exact at n=1 iff normal;
    identity  deficit(n=1) == (1/a0) log(||A||/rho(A))  verified to machine precision.
[C] exact one-plaquette / 2D-U(1) rotor on the circle; Haar = dphi is FLAT, no singularity.
    L_beta = -M_beta^{-1} A.  theta=0 => L symmetric => dep=0, kappa=1.  beta-sweep.
    (ABELIAN model => the defect is NOT nonabelian-essential.)
[D] conversion cost cond(e^{-W/2}) = exp((Wmax-Wmin)/2): bounded on the circle, INFINITE for
    an unbounded W => the flat-frame resolvent bound is worthless in the limit.
[E] Rouche perturbation is not volume-uniform: ||T^{tensor V} - T0^{tensor V}|| grows like V.
[F] finite-moment blindness: competitor on [0,inf) matching M_0..M_6 of delta_1, inf supp = 0.
"""
import numpy as np
from numpy.linalg import norm, eigvals, matrix_power
from numpy.linalg import eigh as symeig


def hdr(s):
    print("\n" + "=" * 76 + "\n" + s + "\n" + "=" * 76)


def dist_spec(z, T):
    return min(abs(z - e) for e in eigvals(T))


def kap(z, T):
    return dist_spec(z, T) * norm(np.linalg.inv(T - z * np.eye(T.shape[0])), 2)


# ----------------------------------------------------------------------------- [A]
hdr("[A] 2x2: identical spectrum {1, e^-m}, different normality.  a0 = 1")
m, lam = 0.6, np.exp(-0.6)
zmid, delta = 0.5 * (1 + lam), 0.5 * (1 - lam)
print(f"  sigma = {{1, {lam:.6f}}};  half-gap delta = {delta:.6f};  z_mid = {zmid:.6f}")
print("    c    dep(T)   ||T||/rho  omega-rho  kappa(z_mid)  min_{|z-zmid|=0.4d}||(T-z)^-1||")
for c in [0.0, 0.1, 0.5, 1.0, 2.0]:
    T = np.array([[1.0, c], [0.0, lam]])
    dep = np.sqrt(norm(T.T @ T - T @ T.T, 2))
    rho = max(abs(eigvals(T)))
    om = max(abs(eigvals(0.5 * (T + T.T))))
    ang = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    on = min(norm(np.linalg.inv(T - (zmid + 0.4 * delta * np.exp(1j * t)) * np.eye(2)), 2) for t in ang)
    print(f"  {c:4.1f} {dep:9.4f} {norm(T,2)/rho:10.4f} {om-rho:10.4f} {kap(zmid,T):12.4f} {on:14.4f}")
print("  -> non-normality INFLATES kappa and every resolvent norm; it never lowers the floor.")

# ----------------------------------------------------------------------------- [B]
hdr("[B] Gelfand deficit  m >= -(1/(n a0)) log ||A^n||^{1/n};  A = P_{vac^perp} T P_{vac^perp}")
a0 = 1.0
print(f"  normal sector A = diag(e^-m,0),  m_true = {m}")
An = np.diag([lam, 0.0])
for n in [1, 2, 4, 8, 16]:
    s = norm(matrix_power(An, n), 2) ** (1.0 / n)
    print(f"    n={n:3d}  ||A^n||^(1/n)={s:.12f}  bound={-np.log(s):.12f}  deficit={m+np.log(s):.2e}")
print("  non-normal sector A = [[lam,c],[0,0]] (same spectrum for every c)")
for c in [0.5, 1.0, 5.0]:
    A = np.array([[lam, c], [0.0, 0.0]])
    kA = norm(A, 2) / max(abs(eigvals(A)))
    print(f"    c={c}:  kappa(A)=||A||/rho={kA:.6f}  (1/a0)log kappa={np.log(kA):.6f}"
          f"  deficit(n=1)={m+np.log(norm(A,2)):.6f}")
    for n in [1, 2, 4, 8, 16, 64]:
        s = norm(matrix_power(A, n), 2) ** (1.0 / n)
        print(f"       n={n:3d}  ||A^n||^(1/n)={s:.8f}  bound={-np.log(s):.6f}  deficit={m+np.log(s):.6f}")

# ----------------------------------------------------------------------------- [C]
hdr("[C] exact one-plaquette / 2D-U(1) rotor on the circle;  Haar = dphi is FLAT")
n = 256
phi = np.arange(n) * 2 * np.pi / n
h = 2 * np.pi / n
mid = phi + 0.5 * h
D1 = (np.roll(np.eye(n), -1, axis=1) - np.eye(n)) / h


def rotor(beta):
    w = np.exp(beta * (np.cos(mid) - 1.0))
    A = D1.T @ (w[:, None] * D1)
    M = np.exp(beta * (np.cos(phi) - 1.0))
    return -np.diag(1.0 / M) @ A, M, A


print("  beta  ||L-L^T||_F/||L||_F  ||[L,L^T]||_F/||L||_F^2  kappa(z=0.05)  kappa(z=0.5)"
      "  kappa(z=0.005)  cond")
for beta in [0.0, 0.5, 1.0, 2.0, 2.5, 4.0]:
    L, M, A = rotor(beta)
    dep = norm(L - L.T) / norm(L)
    com = norm(L @ L.T - L.T @ L) / norm(L) ** 2
    print(f"  {beta:4.1f}  {dep:20.6e}  {com:24.6e}  {kap(0.05,L):13.6f}  {kap(0.5,L):12.6f}"
          f"  {kap(0.005,L):14.6f}  {np.sqrt(M.max()/M.min()):9.3f}")
print("  -> dep == 0 and kappa == 1 EXACTLY at beta = 0 (free / Haar point): the operator is normal.")
print("  -> both > 0 for beta > 0, and this model is ABELIAN: the defect is NOT nonabelian-essential.")
print("     (charter A2/F2: the mechanism therefore fails the nonabelian-essentialness axiom.)")

# ----------------------------------------------------------------------------- [D]
hdr("[D] conversion cost cond(e^{-W/2}) = exp((Wmax-Wmin)/2)")
print("  (i) rotor, W = -beta(cos phi - 1) on the circle: BOUNDED")
for beta in [0.0, 2.5, 4.0]:
    print(f"      beta={beta:4.1f}   cond = e^{{{2*beta:.2f}}} = {np.exp(beta):.4e}")
print("  (ii) confining W = (g/2)(x^2-1)^2 on [-R,R]: UNBOUNDED in R")
for R in [2.0, 3.0, 4.0, 6.0, 8.0, 10.0]:
    Wmax = 0.5 * 0.5 * (R ** 2 - 1) ** 2
    print(f"      g=0.5  R={R:5.1f}  W_max={Wmax:10.2f}   cond = e^{{{0.5*Wmax:10.2f}}}"
          f" = {np.exp(min(0.5*Wmax, 700)):.3e}" + ("   [overflow]" if 0.5*Wmax > 700 else ""))
print("  -> W unbounded => cond = infinity: in the flat frame the resolvent bound says NOTHING;")
print("     in the weighted frame the same operator is exactly self-adjoint (kappa == 1).")

# ----------------------------------------------------------------------------- [E]
hdr("[E] Rouche perturbation is not volume-uniform")
aT = 0.002
L2, M2, _ = rotor(2.5)
L0, M0, _ = rotor(0.0)


def Tof(L, M):
    S = (np.sqrt(M)[:, None] * L) / np.sqrt(M)[None, :]
    ev, U = symeig(0.5 * (S + S.T))
    eS = U @ np.diag(np.exp(-aT * ev)) @ U.T
    sM = np.sqrt(M)
    return (1.0 / sM)[:, None] * eS * sM[None, :], np.sort(ev)


T2, ev2 = Tof(L2, M2)
T0, ev0 = Tof(L0, M0)
T2n, T0n = T2 / norm(T2, 2), T0 / norm(T0, 2)
dT = norm(T2n - T0n, 2)
gap1 = ev2[-1] - ev2[-2]                       # ev ascending: ev[-1] is the vacuum (~0)
dist0 = 0.5 * (1.0 - np.exp(-aT * gap1))
print(f"  normalized one-site operators: ||T(2.5)-T(0)|| = eps1 = {dT:.6f}")
print(f"  contour inside the gap: dist(edge, sigma) = {dist0:.6f}   (E1 = {gap1:.4f})")
print(f"  single site:  eps1/dist = {dT/dist0:.4f}")
print("     V      V*eps1/dist   (V independent sites: ||Delta T^{tensor V}|| = V*eps1 + O(eps1^2))")
for V in [1, 2, 4, 8, 16, 256, 65536]:
    print(f"  {V:7d}   {V*dT/dist0:16.4f}")
print("  -> ||E|| is EXTENSIVE (~V); dist(contour,sigma) is not.  Rouche fails for macroscopic V.")
print("     Same obstruction as the banked e^{-cV} constants of PROOF-GAUSS-INVARIANT-PF.")

# ----------------------------------------------------------------------------- [F]
hdr("[F] finite-moment blindness: competitor matching M_0..M_4 of a gapped 2-atom nu, inf supp = 0")
# target nu = (1/2) delta_1 + (1/2) delta_4   (gapped; the 1-atom case is RIGID: a single atom is
# pinned by M_0,M_1, so the blindness claim needs >=2 atoms -- which is the physical case.)
tgt = [(0.5, 1.0), (0.5, 4.0)]
Mtarget = np.array([sum(w * s ** j for w, s in tgt) for j in range(7)])
print(f"  target nu: M_j = {np.array2string(Mtarget, precision=8)}")
try:
    from scipy.optimize import fsolve
    K = 2                                        # nodes 0, t1, t2 ; match M_0..M_4
    def resid(x):
        w0, w, t = x[0], x[1:1 + K], x[1 + K:]
        if np.any(t <= 0) or np.any(t >= 5):
            return np.full(2 * K + 1, 1e3)
        return np.array([(w0 if j == 0 else 0.0) + np.sum(w * t ** j) for j in range(2 * K + 1)]) - Mtarget[:2 * K + 1]
    best = None
    for seed in range(400):
        rng = np.random.default_rng(seed)
        x0 = np.concatenate([[rng.uniform(0.01, 0.5)], rng.uniform(0.05, 2.0, K),
                             np.sort(rng.uniform(0.05, 4.0, K))])
        sol, info, ier, msg = fsolve(resid, x0, full_output=True)
        if ier == 1 and sol[0] > 1e-3 and np.all(sol[1:1 + K] > 1e-6) and np.all(sol[1 + K:] > 1e-3):
            res = float(np.max(np.abs(resid(sol))))
            if best is None or res < best[0]:
                best = (res, sol)
    if best is None:
        print("  no admissible competitor found by fsolve")
    else:
        res, sol = best
        w0, w, t = sol[0], sol[1:1 + K], np.sort(sol[1 + K:])
        Mj = np.array([(w0 if j == 0 else 0.0) + np.sum(w * t ** j) for j in range(2 * K + 1)])
        print(f"  competitor nu' = {w0:.6e}*delta_0 + "
              + " + ".join(f"{wi:.6e}*delta_{{{ti:.6f}}}" for wi, ti in zip(w, t)))
        print(f"  M_j(nu') = {np.array2string(Mj, precision=8)}")
        print(f"  max |M_j(nu') - M_j(nu)|, j = 0..{2*K} :  {np.max(np.abs(Mj - Mtarget[:2*K+1])):.3e}")
        print(f"  total mass nu' = {w0 + np.sum(w):.10f},  w0 = {w0:.6f}")
        print(f"  inf supp nu' = 0   vs   inf supp nu = 1   -> gap 0 vs 1, data identical to 4 moments.")
        print("  -> finitely many moments (hence ANY finite-n Gelfand datum) cannot see the floor.")
except ImportError:
    print("  scipy unavailable; skipped")

print("\nDONE")
