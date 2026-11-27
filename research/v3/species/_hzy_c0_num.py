# HZY-C0 numerics: deterministic quadrature only (zero MC).
import math

def simpson(f, a, b, n=20000):
    if n % 2: n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4.0 * f(a + i * h)
    for i in range(2, n, 2):
        s += 2.0 * f(a + i * h)
    return s * h / 3.0

def I0(z):
    # modified Bessel I0 via power series
    s = 1.0
    term = 1.0
    z2 = (z / 2.0) ** 2
    for k in range(1, 400):
        term *= z2 / (k * k)
        s += term
        if term < 1e-18 * s:
            break
    return s

def moments_single_plaq(beta):
    # density on x in [-1,1]:  e^{beta x} sqrt(1-x^2)
    def w(x):
        return math.exp(beta * x) * math.sqrt(max(0.0, 1.0 - x * x))
    Z = simpson(w, -1.0, 1.0)
    m1 = simpson(lambda x: x * w(x), -1.0, 1.0) / Z
    m2 = simpson(lambda x: x * x * w(x), -1.0, 1.0) / Z
    return m1, m2, m2 - m1 * m1, Z

def moments_link_tilt(c, phi):
    # u in S^3, density propto exp(c*(cos phi*u0 + sin phi*u1))
    # marginal of t=u0: p(t) ~ e^{lam t} I0(mu sqrt(1-t^2)), lam=c cos phi, mu=c sin phi
    lam = c * math.cos(phi)
    mu = c * math.sin(phi)
    def w(t):
        return math.exp(lam * t) * I0(mu * math.sqrt(max(0.0, 1.0 - t * t)))
    Z = simpson(w, -1.0, 1.0)
    m1 = simpson(lambda t: t * w(t), -1.0, 1.0) / Z
    m2 = simpson(lambda t: t * t * w(t), -1.0, 1.0) / Z
    return m2 - m1 * m1

print("== single plaquette / 2D YM : Var(q) = Var(x), x=Re tr U /2 ==")
for beta in [0.0, 0.5, 1.0, 2.0, 2.5, 3.0]:
    m1, m2, var, Z = moments_single_plaq(beta)
    print(f"  beta={beta:4.2f}  <x>={m1:+.8f}  <x^2>={m2:.8f}  Var(q)={var:.8f}")

print()
print("== vMF covariance at c = 6 beta : a=(1-m2)/3 (transverse), b=m2-m1^2 (along h) ==")
for beta in [1.0, 2.0, 2.5, 3.0]:
    c = 6.0 * beta
    m1, m2, var, Z = moments_single_plaq(c)   # same density shape (x in [-1,1])
    a = (1.0 - m2) / 3.0
    b = m2 - m1 * m1
    print(f"  beta={beta:4.2f} c={c:5.2f}: m1={m1:.8f} m2={m2:.8f}  a={a:.8f}  b={b:.8f}  min={min(a,b):.8f}")

print()
print("== single-link conditional floor  v(beta) = min_phi Var(u0), |c|=2(d-1)beta, d=4 => c=6beta ==")
for beta in [0.5, 1.0, 2.0, 2.5, 3.0]:
    c = 6.0 * beta
    best = None
    rows = []
    for k in range(0, 91, 5):
        phi = math.pi * k / 180.0
        v = moments_link_tilt(c, phi)
        rows.append((k, v))
        if best is None or v < best[1]:
            best = (k, v)
    print(f"  beta={beta:4.2f} c={c:5.2f}:  phi scan (deg): " +
          " ".join(f"{k}:{v:.5f}" for k, v in rows))
    print(f"      -> min at phi={best[0]} deg, v={best[1]:.8f}")
