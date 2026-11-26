"""o06_hvv_bound.py -- O06 curvature side: exact H_vv, Phi'' floor, composition.

Model A (H1A-EXACT Sec 2.1 convention, single plaquette SU(2), flat chart X=log U):
    xi in R^3, |xi|<2pi ; B = exp(i theta sigma_3/2)
    S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  u = xi_3/rho
    J = sin^2(rho/2)/rho^2 ,  Phi = -log J        (Phi = -log|c_SU2/c_flat|^2)
    conditioned slice xi_3 = 0 (v = e_3, a_0 = 0), measure e^{-S} J 2pi rho drho
Deterministic Gauss-Legendre quadrature, no sampling.  Run: py -3.12 w7/o06_hvv_bound.py
"""
import sys
import numpy as np

TP = 2.0 * np.pi
NQ = 400  # Gauss-Legendre nodes


def Phi(rho):
    return -2.0 * np.log(np.sin(rho / 2.0) / rho)


def Jf(rho):
    return np.sin(rho / 2.0) ** 2 / rho**2


def Phi2(rho):
    """Phi'' along e_3 at xi_3=0 : (dPhi/drho)/rho = 2/rho^2 - cot(rho/2)/rho."""
    return 2.0 / rho**2 - 1.0 / (np.tan(rho / 2.0) * rho)


def Hvv(rho, beta, th):
    """H_vv on the slice = beta cos(th/2) sin(rho/2)/(2 rho) = (beta/2)cos(th/2) sqrt(J)."""
    return beta * np.cos(th / 2.0) * np.sin(rho / 2.0) / (2.0 * rho)


def quad(f, beta, th, n=NQ):
    x, w = np.polynomial.legendre.leggauss(n)
    rho = 0.5 * TP * (x + 1.0)
    jac = 0.5 * TP * w
    S = beta * (1.0 - np.cos(th / 2.0) * np.cos(rho / 2.0))
    W = np.exp(-S) * Jf(rho) * TP * rho * jac
    return float(np.dot(W, f(rho))), float(np.sum(W))


def force(rho, beta, th):
    """G_v = dS/dxi_3 at xi_3=0 (exact): beta sin(rho/2) sin(th/2)/rho."""
    return beta * np.sin(rho / 2.0) * np.sin(th / 2.0) / rho


def moments(beta, th, n=NQ):
    Z, _ = quad(lambda r: np.ones_like(r), beta, th, n)
    EH, _ = quad(lambda r: Hvv(r, beta, th), beta, th, n)
    Ep, _ = quad(Phi2, beta, th, n)
    Ef, _ = quad(lambda r: force(r, beta, th), beta, th, n)
    Ef2, _ = quad(lambda r: force(r, beta, th) ** 2, beta, th, n)
    EH, Ep, Ef, Ef2 = EH / Z, Ep / Z, Ef / Z, Ef2 / Z
    return dict(EH=EH, Ep=Ep, Ef=Ef, var=Ef2 - Ef**2, W=Ep, Z=Z)


def fd_H(beta, th, rho, s=1e-5):
    def S(xi):
        r = np.linalg.norm(xi)
        return beta * (1.0 - np.cos(r / 2) * np.cos(th / 2)
                       + np.sin(r / 2) * np.sin(th / 2) * xi[2] / r)
    x = np.array([rho, 0.0, 0.0])
    e = np.array([0.0, 0.0, 1.0])
    return (S(x + s * e) - 2 * S(x) + S(x - s * e)) / s**2


def fd_Phi(rho, s=1e-5):
    f = lambda xi: Phi(np.linalg.norm(xi))
    x = np.array([rho, 0.0, 0.0])
    e = np.array([0.0, 0.0, 1.0])
    return (f(x + s * e) - 2 * f(x) + f(x - s * e)) / s**2


def su2(rng):
    """hermitian traceless a.sigma (the *coordinate* A with chart X = i A)."""
    a = rng.normal(size=4)
    a /= np.linalg.norm(a)
    s = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]])]
    return a[0] * np.eye(2, dtype=complex) + sum(a[i + 1] * s[i] for i in range(3))


def expx(X):
    """exp of an antihermitian traceless X (= i A)."""
    ev, U = np.linalg.eigh(-1j * X)
    return U @ np.diag(np.exp(1j * ev)) @ U.conj().T


def modelB(nsamp=120, beta=2.5, seed=7, s=1e-4):
    """4-link single plaquette, direction v constant in the chart (v = i sigma_3/2)."""
    rng = np.random.default_rng(seed)
    v = 1j * np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex) / 2.0
    hs, als = [], []
    for _ in range(nsamp):
        X = [1j * (0.6 * su2(rng)) for _ in range(4)]

        def hol(scale):
            L = [expx(X[l] + scale * v) for l in range(4)]
            return L[0] @ L[1] @ L[2].conj().T @ L[3].conj().T

        Sf = lambda sc: beta * (1.0 - 0.5 * np.real(np.trace(hol(sc))))
        hs.append((Sf(s) - 2.0 * Sf(0.0) + Sf(-s)) / s**2)
        U0, Up, Um = hol(0.0), hol(s), hol(-s)
        Wp = (Up - Um) / (2 * s) @ np.linalg.inv(U0)
        Hp = -1j * Wp
        Hp = 0.5 * (Hp + Hp.conj().T)
        K = Hp @ Hp
        tK = np.real(np.trace(K))
        if tK > 1e-14:
            als.append(np.real(np.trace(U0 @ K)) / tK)
    return np.array(hs), np.array(als)


p = lambda *a: print(*a, flush=True)

p("=" * 78)
p("O06  H_vv exact form / Phi'' floor / composition    py3.12 numpy, GL n=%d" % NQ)
p("=" * 78)

p("\n[1] exact identity  H_vv|slice = beta*cos(th/2)*sin(rho/2)/(2 rho) = (beta/2)cos(th/2)*sqrt(J)")
for (b, th, r) in [(2.5, 0.0, 1.7), (2.5, 1.1, 1.7), (2.5, np.pi, 1.7), (6.0, 0.3, 2.9)]:
    fd, ex = fd_H(b, th, r), Hvv(r, b, th)
    p(f"  beta={b:<4} th={th:8.5f} rho={r:<4} FD={fd: .12e} formula={ex: .12e} rel={abs(fd-ex)/max(abs(ex),1e-30):.1e}")

p("\n[2] Phi'' = 2/rho^2 - cot(rho/2)/rho  (= g'(rho)/rho, radial Hessian transverse)")
for r in [1e-3, 0.5, np.pi, 3.0, 6.0]:
    p(f"  rho={r:<8.5g} FD={fd_Phi(r): .10e} formula={float(Phi2(np.array([r]))[0]): .10e}")
rr = np.linspace(1e-6, TP * (1 - 1e-9), 400001)
vv = Phi2(rr)
p(f"  inf Phi'' = {np.nanmin(vv):.10f} at rho={rr[np.nanargmin(vv)]:.4g}   1/6={1/6:.10f}  N/12(N=2)=%.10f" % (2 / 12))
p(f"  Phi'' at rho=2,3,4,5,6: {[round(float(Phi2(np.array([x]))[0]),6) for x in [2,3,4,5,6]]}")

p("\n[3] conditioned-slice moments (v=e_3, a_0=0): E[H_vv|0], E[Phi''|0], Var(force|0)")
p("  beta      th       E[H_vv|0]     E[Phi''|0]    Sum           Var(G_v|0)    Sum-1/6")
tab = {}
for b in [2.5, 4.30661, 6.0]:
    for th in [0.0, np.pi / 2, np.pi]:
        m = moments(b, th)
        s = m["EH"] + m["Ep"]
        tab[(b, th)] = (m, s)
        p(f"  {b:<8.5g} {th:<7.5f} {m['EH']: .6e} {m['Ep']: .6e} {s: .6e} {m['var']: .6e} {s-1/6: .3e}")

m0, s0 = tab[(2.5, 0.0)]
mp, sp = tab[(2.5, np.pi)]
p(f"\n  theta=pi (Savvidy pt): E[H_vv|0]={mp['EH']:.3e} (c=cos(pi/2)=0 exact)"
  f"  E[Phi''|0]={mp['Ep']:.6f}  Var={mp['var']:.6f}  floor/var={(1/6)/mp['var']:.3f}")
p(f"  theta=0  (aligned)   : E[H_vv|0]={m0['EH']:.6f}  E[Phi''|0]={m0['Ep']:.6f}"
  f"  Var={m0['var']:.6f}  Var-1/6={m0['var']-1/6: .3e}")
p(f"  theta=pi: E[H+Phi'']-1/6 = {sp-1/6: .6e}  (floor holds with equality iff H_vv=0 & Phi''=1/6)")

p("\n[3b] quadrature stability (n=200 vs 400), beta=2.5 th=pi")
for n in [200, 400]:
    mm = moments(2.5, np.pi, n)
    p(f"   n={n:<4} E[Phi'']={mm['Ep']:.12f}  Var={mm['var']:.12f}")

p("\n[4] model B (4 links, constant chart direction, beta=2.5): is H_vv indefinite?")
hs, als = modelB()
p(f"  n={len(hs)}  min H_vv={hs.min(): .6e}  max={hs.max(): .6e}  mean={hs.mean(): .6e}"
  f"  frac<0={np.mean(hs < 0):.3f}")
p(f"  alignment Re tr(U_p K_p)/tr K_p : min={als.min(): .6f}  max={als.max(): .6f}"
  f"  mean={als.mean(): .6f}  frac<0={np.mean(als < 0):.3f}")
p("  => H_vv = beta*sum_p (1/2)Re tr(U_p K_p) - beta*sum_p (1/2)Re tr(U_p W'_p),  K_p=H_p^2>=0")
p("     K_p>=0 = positivity ; U_p = geometric (sign-carrying) factor ;")
p("     negative alignment => H_vv<0 pointwise => no positive pointwise floor.")
sys.stdout.flush()
