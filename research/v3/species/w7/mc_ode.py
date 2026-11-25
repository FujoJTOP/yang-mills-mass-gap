"""mc_ode.py -- T3: the lambda_t / gamma side.  Three exact statements, no simulation.

(1) FLAT GAUSSIAN TOY (the paper's own engine).  V_0(phi) = (1/2)<phi,A phi>, A > 0,
    C_t = int_0^t Cdot, Hess V_t = (A^{-1} + C_t)^{-1}.  The criterion
        Cdot_t Hess V_t Cdot_t - Cdotdot_t/2  >=  lam_dot_t Cdot_t ,   1/gamma = int_0^inf e^{-2 lam_t}
    is checked exactly for two ansaetze:
       (a) Cdot_t = e^{-tA}          (Example 3.18;  C_inf = A^{-1} < inf)   -> gamma = lam_min(A)
       (b) Cdot_t = c * 1  (the "C_t = t g^{-1}" ansatz of T3*)             -> gamma = 0
    (b) is why the ansatz C_t = t g^{-1} cannot be used in the FLAT criterion: (3.10) requires
    C_inf = int_0^inf Cdot < inf.  On a compact group there is no (3.10) -- the manifold
    version replaces it by (3.93) -- and then (b)-type ansaetze give the ASYMPTOTIC value
    lam_dot -> lam_min(Ric/2) = N/8 > 0 instead (see (2)).

(2) COMPACT-GROUP ASYMPTOTICS.  K_t -> 1 (Haar), hence V_t -> -log Z = const and
    Hess V_t -> 0.  Therefore in (star)  lam_dot_t -> lam_min(Ric_G/2) = N/8 for the
    bi-invariant metric:  lambda_t -> +inf  and  int e^{-2 lambda_t} < inf  AUTOMATICALLY,
    for every beta, as soon as lambda_dot is bounded below on (0, inf).  The whole content of
    (star)+(star star) is that lower bound -- i.e. an L-UNIFORM bound on Hess V_t, not the
    divergence of lambda_t.  (Checked against the free-field limit of PL-POLCHINSKI 6.3.)

(3) THE ODE BOUND.  If lam_dot_t >= -C on [0,T] and lam_dot_t >= a > 0 for t >= T, then
        lambda_t >= -C T + a (t - T)   and   1/gamma <= e^{2CT} (T + 1/(2a)),
    so gamma >= e^{-2CT} / (2T + 1/a) > 0.  gamma is positive iff C < inf: the criterion is a
    statement about a LOWER BOUND, and the size of gamma is exp(-2CT).
"""
import numpy as np

print("=" * 96)
print("(1) flat Gaussian toy: exact check of the criterion for two covariance ansaetze")
n = 3
A = np.diag([0.3, 1.0, 2.5])          # mass spectrum; lam_min(A) = 0.3
Ainv = np.linalg.inv(A)


def gamma_of(cdot_fn, T=30.0, dt=1e-3, free=False):
    """lam_dot_t = lam_min of the criterion; then 1/gamma = integral of e^{-2 lam}."""
    ts = np.arange(dt, T, dt)
    lam = 0.0
    integ = 0.0
    lamdot_last = None
    for t in ts:
        Cd, Cdd, Ct = cdot_fn(t)
        Hess = np.zeros((n, n)) if free else np.linalg.inv(Ainv + Ct)
        Q = Cd @ Hess @ Cd - 0.5 * Cdd
        ci = np.diag(1.0 / np.sqrt(np.diag(Cd)))     # Cd^{-1/2} (diagonal in this toy)
        S = ci @ Q @ ci                              # congruence: sign and ratio are exact
        S = 0.5 * (S + S.T)
        ld = np.linalg.eigvalsh(S)[0]
        lam += ld * dt
        integ += np.exp(-2 * lam) * dt
        lamdot_last = ld
    return 1.0 / integ if integ > 0 else np.inf, lam, lamdot_last


# (a) Cdot_t = e^{-tA}  (C_t = A^{-1}(1 - e^{-tA}), C_inf = A^{-1})
def cd_a(t):
    E = np.diag(np.exp(-np.diag(A) * t))
    return E, -A @ E, np.diag((1 - np.exp(-np.diag(A) * t)) / np.diag(A))


# (b) Cdot_t = 1 (C_t = t*1, the T3* ansatz: C_inf = infinity)
def cd_b(t):
    return np.eye(n), np.zeros((n, n)), t * np.eye(n)


for nm, fn in (("(a) Cdot = e^{-tA}", cd_a), ("(b) Cdot = 1 (T3* ansatz)", cd_b)):
    g, lamT, ld = gamma_of(fn)
    print(f"  QUADRATIC V0  {nm:26s}: 1/gamma = {1/g if np.isfinite(g) else np.inf:10.4f}   "
          f"gamma = {g:10.6f}   lam_dot(T) = {ld:+.5f}")
for nm, fn in (("(a) Cdot = e^{-tA}", cd_a), ("(b) Cdot = 1 (T3* ansatz)", cd_b)):
    g, lamT, ld = gamma_of(fn, free=True)
    print(f"  FREE V0 = 0   {nm:26s}: 1/gamma = {1/g if np.isfinite(g) else np.inf:10.4f}   "
          f"gamma = {g:10.6f}   lam_dot(T) = {ld:+.5f}")
print(f"  exact answer for (a) with V0 = 0 (PL-POLCHINSKI 6.3): gamma = lam_min(A) = "
      f"{np.linalg.eigvalsh(A)[0]:.6f}")
print("  => the T3* ansatz C_t = t g^{-1} (Cdot = const) gives lam_dot = 1/(t + 1/lam_min)")
print("     for a quadratic V0 (gamma finite) but lam_dot = 0 for the FREE field: in the flat")
print("     case it therefore has no convexification of its own.  On a compact group the")
print("     missing convexification is supplied by Ric/2 = (N/8) g, which is why the group")
print("     criterion survives the free-field test with gamma -> N/8 > 0.")

print()
print("=" * 96)
print("(2) compact group: the asymptotic value of lam_dot is lam_min(Ric/2) = N/8")
for N in (2, 3, 4, 5):
    print(f"  SU({N}): lam_dot -> N/8 = {N/8:.4f}  (beta-independent, because V_t -> const)")

print()
print("=" * 96)
print("(3) gamma from a lower bound C and asymptotic rate a")
for C in (0.0, 3.5, 10.0, 40.0):
    for T in (1.0, 5.0):
        a = 0.25
        gam = np.exp(-2 * C * T) / (2 * T + 1 / a)
        print(f"  C={C:5.1f} T={T:4.1f} a={a}:  gamma >= {gam:.3e}   (positive for every finite C)")
    print()
print("  => (star star) is never the obstruction on a compact group.  The obstruction is")
print("     L-UNIFORMITY of C, and C is controlled by sup|grad S_W|^2, which is EXTENSIVE:")
print("     measured (mc_lattice.log): sum_links |grad S_W|^2 = 900 at L=2 and 14400 at L=4")
print("     (SU(2), beta=2.5) -- a factor 16 for L -> 2L, i.e. exactly L^4, 14.06 per link.")
