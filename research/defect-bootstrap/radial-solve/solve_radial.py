"""Numerical spectral solve attempt for the defect radial Casimir block.

Two objects:
(A) The 1D hypergeometric block f(x)=x^D 2F1(D,D;2D;x), which every sub-agent agrees
    is the 'xi-part' / product factor.  Solve the hypergeometric ODE spectrally
    (Chebyshev) for g=2F1 and verify against scipy.hyp2f1.
(B) The H^2 Laplacian eigen-equation in the K3-invariant variable z=cosh d (and
    u=(z-1)/(z+1)), with the naive Casimir eigenvalue D(D-1).  Check whether it is a
    regular (hypergeometric-type) ODE or irregular -- i.e. whether 'solve the single
    ODE in u' is well-posed.
"""
import numpy as np
from scipy.special import hyp2f1
from scipy.integrate import solve_ivp


def cheb_D(N):
    if N == 0:
        return np.zeros((1, 1)), np.array([1.0])
    k = np.arange(N + 1)
    x = np.cos(np.pi * k / N)          # x[0]=1, x[N]=-1
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** k
    X = np.tile(x[:, None], (1, N + 1))
    dX = X - X.T
    D = (c[:, None] * (1.0 / c)[None, :]) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def solve_hypergeom(Dval, a=1e-9, b=0.85, N=160):
    """Solve x(1-x) g'' + [2D-(2D+1)x] g' - D^2 g = 0, g(a)=1, g'(a)=D/2."""
    Dd, xc = cheb_D(N)                 # xc in [-1,1], xc[N]=-1 -> x=a
    mid, half = (a + b) / 2, (b - a) / 2
    t = mid + half * xc
    Dt = Dd / half
    L = Dt @ Dt
    ab = Dval * Dval
    csz = 2 * Dval
    A = np.zeros((N + 1, N + 1)); rhs = np.zeros(N + 1)
    for i in range(N - 1):             # ODE on rows 0..N-2
        tt = t[i]
        A[i, :] = tt * (1 - tt) * L[i, :] + (csz - (ab + 1) * tt) * Dt[i, :] - ab * np.eye(N + 1)[i, :]
    A[N - 1, :] = np.eye(N + 1)[N, :]  # g at t=a
    rhs[N - 1] = 1.0
    A[N, :] = Dt[N, :]                 # g' at t=a
    rhs[N] = Dval / 2.0
    g = np.linalg.solve(A, rhs)
    return t, g


def check_h2_ode(Dval, sign=+1):
    """(z^2-1) f'' + z f' + sign*D(D-1) f = 0, shoot from z0=1.01 (f=1,f'=0)."""
    lam = Dval * (Dval - 1)

    def rhs(z, y):
        f, fp = y
        if abs(z * z - 1) < 1e-12:
            return [fp, 0.0]
        return [fp, (-z * fp - sign * lam * f) / (z * z - 1.0)]

    sol = solve_ivp(rhs, [1.01, 2.5], [1.0, 0.0], max_step=0.01, rtol=1e-8, atol=1e-10)
    z = sol.t; f = sol.y[0]
    return z, f


print("=== (A) hypergeometric block, spectral vs hyp2f1 ===")
for Dval in [1, 2, 3]:
    t, g = solve_hypergeom(Dval)
    m = t > 1e-8
    err = np.max(np.abs(g[m] - hyp2f1(Dval, Dval, 2 * Dval, t[m])))
    # tabulate f = t^D * g at sample x
    samples = np.array([0.1, 0.3, 0.5, 0.7])
    fg = np.interp(samples, t[::-1], g[::-1])   # t ascending? xc desc, t desc -> reverse
    fvals = samples ** Dval * hyp2f1(Dval, Dval, 2 * Dval, samples)
    print(f"Delta={Dval}: max|spectral-hyp2f1| = {err:.2e} on t in ({t.min():.1e},{t.max():.2f})")
    print(f"   f=x^D*2F1 at x={list(samples)}: {np.round(fvals,6)}")

print("\n=== (B) H^2 Laplacian ODE  (z^2-1)f'' + z f' +/- D(D-1) f = 0 ===")
for Dval in [2, 3]:
    for sign in [+1, -1]:
        z, f = check_h2_ode(Dval, sign)
        fin = np.isfinite(f).all()
        print(f"Delta={Dval}, sign={sign:+d}: f finite={fin}, f(z=1.01)={f[0]:.3e}, "
              f"f(z=1.5)={np.interp(1.5,z,f):.3e}, f(z=2.5)={f[-1]:.3e}")

# Indicial analysis at z=1: substitute z=1+w, look at leading balance
print("\n=== indicial check at z=1 (z=1+w, 2w f'' + f' + sign*lam f = 0) ===")
print("For f=w^alpha: alpha(2alpha-1) w^{alpha-1} + sign*lam w^alpha = 0")
print("=> alpha-1 = alpha impossible for lam!=0  ==> IRREGULAR unless lam=0")
