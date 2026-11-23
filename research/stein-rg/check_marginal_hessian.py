"""Does the *global* Hessian lower bound survive marginalization (blocking)?

Model (globally strongly convex for c>|gamma|):
    F(x,y) = (c/2)(x^2+y^2) + gamma*x*y + (delta/4)(x^4+y^4)
    F'' = [[c+3 delta x^2, gamma],
           [gamma,          c+3 delta y^2]]
    lambda := inf min-eig F'' = c - |gamma|   (attained at x=y=0)

Marginal G(x) = -log int exp(-F(x,y)) dy.  Test whether G''(x) >= lambda.

Formula (exact):  G''(x) = E[d2xx F | x] - Var(d_x F | x).
If G'' can drop below lambda, strong log-concavity is NOT preserved, and the
obstruction is exactly the variance (non-Gaussianity) term.
"""
import numpy as np

def Fval(x, y, c, gam, delta):
    return 0.5 * c * (x * x + y * y) + gam * x * y + 0.25 * delta * (x ** 4 + y ** 4)

def min_eig(x, y, c, gam, delta):
    a = c + 3.0 * delta * x * x
    d = c + 3.0 * delta * y * y
    b = gam
    tr, det = a + d, a * d - b * b
    disc = np.sqrt(np.maximum(tr * tr - 4.0 * det, 0.0))
    return 0.5 * (tr - disc)

def Gx(x, c, gam, delta, ys):
    e = -Fval(x, ys, c, gam, delta)
    m = e.max()
    return -(np.log(np.exp(e - m).sum()) + m)

def Gpp(x, c, gam, delta, ys, h=2e-3):
    return (Gx(x + h, c, gam, delta, ys) - 2 * Gx(x, c, gam, delta, ys)
            + Gx(x - h, c, gam, delta, ys)) / h ** 2

def formula(x, c, gam, delta, ys):
    e = -Fval(x, ys, c, gam, delta)
    w = np.exp(e - e.max()); w /= w.sum()
    dx = c * x + gam * ys + delta * x ** 3
    d2 = c + 3.0 * delta * x * x
    Ed2 = np.sum(w * d2)
    Var = np.sum(w * dx ** 2) - np.sum(w * dx) ** 2
    return Ed2 - Var

ys = np.linspace(-12.0, 12.0, 400001)

for (c, gam, delta) in [(1, .3, 1.), (1, .3, 5.), (1, .3, 20.),
                        (1, .3, 100.), (2, .8, 25.), (1, .8, 50.)]:
    lam = c - abs(gam)
    xs = np.linspace(-3, 3, 61)
    gs = np.array([Gpp(x, c, gam, delta, ys) for x in xs])
    worst, wx = gs.min(), xs[gs.argmin()]
    gap = lam - worst
    tag = "COUNTEREXAMPLE (G'' < lambda)" if gap > 2e-3 else "preserved"
    print(f"c={c} gam={gam} delta={delta}: lambda={lam:.4f}  minG''={worst:.4f}"
          f" at x={wx:+.2f}  (lambda-minG''={gap:+.4f})  -> {tag}")

x0, c, gam, delta = 0.7, 1.0, 0.3, 5.0
print("\nformula check: G''(fd) =", round(Gpp(x0, c, gam, delta, ys), 6),
      " E[d2xx]-Var =", round(formula(x0, c, gam, delta, ys), 6))
