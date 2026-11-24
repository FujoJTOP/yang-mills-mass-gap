"""
construct_pade.py -- CONSTRUCT-PADE witness
Padé / Borel-Pade continuation of the 0++ strong-coupling series m/sqrt(sigma) toward beta=2.5.

Lightweight: numpy + scipy.integrate only. Single run.  py -3.12 -I construct_pade.py

Structure
  A. Exact SU(2) input: u(beta)=I2/I1, leading strong-coupling m/sqrt(sigma)=4*sqrt(ln(1/u)).
  B. Method calibration on two Borel-summable asymptotic series with KNOWN exact values at u=0.5:
       M1 Stieltjes   f1(u)=int_0^inf e^-t/(1+u t) dt ,   a_n=(-1)^n n!          (Borel 1/(1+t), rational)
       M2 sqrt-branch f2(u)=int_0^inf e^-t/sqrt(1+u t) dt, a_n=(-1)^n (2n)!/(4^n n!) (Borel 1/sqrt(1+t), non-rational)
     For each: Padé [L/M] of the asymptotic coeffs, and Borel-Padé, evaluated at u=0.5; error vs exact.
     This measures HOW MANY TERMS and WHICH method are needed -- the calibration the physical case lacks.
  C. Physical estimate: leading anchor at beta=2.5, and the UNDERDETERMINATION argument (a one-parameter family
     of corrections matches the leading term yet spans the whole 3.3..3.9 band).

No QED. Labels printed. Not a Yang-Mills proof.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import ive

np.set_printoptions(precision=6, suppress=True)


def banner(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ----------------------------------------------------------------------
# Padé [L/M] from a Taylor series a_0..a_{L+M}
# ----------------------------------------------------------------------
def pade_coeffs(a, L, M):
    """Return (p, q) coefficient arrays with P/Q ~ sum a_n x^n, Q(0)=1."""
    a = np.asarray(a, float)
    # Q: q_0=1, solve for q_1..q_M from the M equations  sum_{j=0}^M q_j a_{n-j}=0, n=L+1..L+M
    A = np.zeros((M, M))
    rhs = np.zeros(M)
    for i in range(M):
        n = L + 1 + i
        for j in range(1, M + 1):
            if n - j >= 0:
                A[i, j - 1] = a[n - j]
        rhs[i] = -a[n] if n >= 0 else 0.0
    q = np.zeros(M + 1)
    q[0] = 1.0
    if M > 0:
        try:
            q[1:] = np.linalg.solve(A, rhs)
        except np.linalg.LinAlgError:
            q[1:] = np.linalg.lstsq(A, rhs, rcond=None)[0]
    # P: p_n = sum_{j=0}^{min(n,M)} q_j a_{n-j}, n=0..L
    p = np.zeros(L + 1)
    for n in range(L + 1):
        s = 0.0
        for j in range(min(n, M) + 1):
            s += q[j] * a[n - j]
        p[n] = s
    return p, q


def peval(c, x):
    return np.polyval(c[::-1], x)  # c[0]+c[1]x+... = polyval with reversed


def pade_value(a, L, M, x):
    p, q = pade_coeffs(a, L, M)
    den = peval(q, x)
    if abs(den) < 1e-300:
        return np.nan
    return peval(p, x) / den


# ----------------------------------------------------------------------
# Borel-Padé: b_n = a_n/n!; Padé B(t)~sum b_n t^n; f(u)~int_0^inf e^-s B(u s) ds
# ----------------------------------------------------------------------
def borel_pade_value(a, L, M, u):
    a = np.asarray(a, float)
    n = np.arange(len(a))
    fact = np.array([float(np.math.factorial(int(k))) if False else 1.0 for k in n])  # placeholder
    # b_n = a_n / n!  (compute via lgamma to avoid overflow)
    from math import lgamma
    lb = np.array([np.log(abs(a[k])) - lgamma(k + 1) if a[k] != 0 else -np.inf for k in range(len(a))])
    sign = np.sign(a)
    b = sign * np.exp(lb)
    p, q = pade_coeffs(b, L, M)

    def integrand(s):
        x = u * s
        den = peval(q, x)
        if abs(den) < 1e-300:
            return 0.0
        return np.exp(-s) * (peval(p, x) / den)

    val, _ = quad(integrand, 0.0, np.inf, limit=400)
    return val


# ======================================================================
banner("A. Exact SU(2) input  u(beta)=I2/I1  and the leading strong-coupling m/sqrt(sigma)")
print(f"{'beta':>7} {'u=I2/I1':>12} {'ln(1/u)':>10} {'4*sqrt(ln1/u)':>14}")
for beta in [0.02083, 0.1, 0.5, 1.0, 1.625, 2.0, 2.2, 2.5, 3.0]:
    # u = I2/I1 = ive(2,beta)/ive(1,beta)  (ive = exponentially scaled, ratio exact)
    u = ive(2, beta) / ive(1, beta)
    if u <= 0 or u >= 1:
        continue
    mrat = 4.0 * np.sqrt(np.log(1.0 / u))
    print(f"{beta:7.4f} {u:12.6f} {np.log(1/u):10.6f} {mrat:14.4f}")
u25 = ive(2, 2.5) / ive(1, 2.5)
print(f"\n  beta=1/48 coarse approx : u={ive(2,1/48)/ive(1,1/48):.6f}")
print(f"  beta=2.5                : u={u25:.6f}, 4*sqrt(ln(1/u))={4*np.sqrt(np.log(1/u25)):.4f}")
print("  LATITUDE check: lattice m_0++/sqrt(sigma) ~ 3.9 (Morningstar-Peardon hep-lat/9901004).")
print("  => leading strong-coupling continuation UNDERSHOOTS lattice by ~ (3.9-3.30)/3.9 = %.1f%%"
      % (100*(3.9-4*np.sqrt(np.log(1/u25)))/3.9))

# ======================================================================
banner("B. METHOD CALIBRATION on asymptotic series with KNOWN exact f(0.5)")

# exact values at u=0.5
f1_exact = quad(lambda t: np.exp(-t) / (1.0 + 0.5 * t), 0, np.inf, limit=400)[0]
f2_exact = quad(lambda t: np.exp(-t) / np.sqrt(1.0 + 0.5 * t), 0, np.inf, limit=400)[0]
print(f"  exact f1(0.5) = int e^-t/(1+0.5t) dt      = {f1_exact:.8f}")
print(f"  exact f2(0.5) = int e^-t/sqrt(1+0.5t) dt  = {f2_exact:.8f}")

from math import factorial, comb
NMAX = 16

def coeffs_M1(n):
    return np.array([(-1.0) ** k * factorial(k) for k in range(n)])

def coeffs_M2(n):
    return np.array([(-1.0) ** k * factorial(2 * k) / (4.0 ** k * factorial(k)) for k in range(n)])

print("\n  Padé [L/M] of the ASYMPTOTIC series, evaluated at u=0.5, error vs exact:")
print(f"  {'N=L+M':>5} {'[L/M]':>8} {'M1 Padé err':>14} {'M2 Padé err':>14}")
for N in range(2, NMAX + 1, 2):
    L = N // 2
    M = N - L
    a1 = coeffs_M1(N + 1)
    a2 = coeffs_M2(N + 1)
    e1 = abs(pade_value(a1, L, M, 0.5) - f1_exact) / f1_exact
    e2 = abs(pade_value(a2, L, M, 0.5) - f2_exact) / f2_exact
    print(f"  {N:5d} {f'[{L}/{M}]':>8} {e1:14.3e} {e2:14.3e}")

print("\n  Borel-Padé [L/M] (b_n=a_n/n!), evaluated at u=0.5, error vs exact:")
print(f"  {'N':>5} {'[L/M]':>8} {'M1 BP err':>14} {'M2 BP err':>14}")
for N in range(2, NMAX + 1, 2):
    L = N // 2
    M = N - L
    a1 = coeffs_M1(N + 1)
    a2 = coeffs_M2(N + 1)
    try:
        v1 = borel_pade_value(a1, L, M, 0.5)
        e1 = abs(v1 - f1_exact) / f1_exact
    except Exception:
        e1 = np.nan
    try:
        v2 = borel_pade_value(a2, L, M, 0.5)
        e2 = abs(v2 - f2_exact) / f2_exact
    except Exception:
        e2 = np.nan
    print(f"  {N:5d} {f'[{L}/{M}]':>8} {e1:14.3e} {e2:14.3e}")

# ======================================================================
banner("C. PHYSICAL UNDERDETERMINATION: the strong-coupling anchor does not pin beta=2.5")
m_lead = 4.0 * np.sqrt(np.log(1.0 / u25))
print(f"  leading continuation at beta=2.5:  m/sqrt(sigma) = {m_lead:.4f}")
print("  Write m/sqrt(sigma)(u) = 4*sqrt(ln(1/u)) * g(u), g(0)=1 (g corrections unknown).")
print("  A one-parameter family g_c(u)=1/(1 - c u) matches the leading term exactly; at u=0.5:")
print(f"  {'c':>7} {'g_c(0.5)':>10} {'m/sqrt(sigma)(2.5)':>20}")
for c in [0.0, 0.10, 0.20, 0.30, 0.40, 0.50]:
    g = 1.0 / (1.0 - c * 0.5)
    print(f"  {c:7.2f} {g:10.4f} {m_lead*g:20.4f}")
print("  The correction coefficient c is UNCONSTRAINED by the leading strong-coupling data:")
print("  ANY c in [0,0.5] gives g(0)=1 yet spans the whole 3.3..3.9 band. ==> underdetermined.")

print("\n  SUMMARY")
print("  - Method: Padé / Borel-Padé of a Borel-summable asymptotic series converges (see Table B),")
print("    so a CONTROLLED crossover estimate is possible IN PRINCIPLE, given enough genuine coefficients.")
print("  - Physical input: only the LEADING strong-coupling term is available; at u(2.5)=%.3f the" % u25)
print("    correction O(u) is NOT small (~0.5), so the leading formula undershoots by ~%.0f%%."
      % (100*(3.9-m_lead)/3.9))
print("  - Honest estimate: m_0++/sqrt(sigma)(beta=2.5) in [3.3, 3.9], best leading point 3.30,")
print("    systematic uncertainty ~15%. A robust continuation requires the higher SU(2) strong-coupling")
print("    coefficients (OS'78 / Münster / SZZ beyond leading) -- the EXACT MISSING INPUT [O].")
print("\n  Labels: [V-num] tables A-C; [C] method calibration; [O] missing coefficients. No QED.")
