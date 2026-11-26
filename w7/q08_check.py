"""scratch: first-derivative identity check for the Q08 reduction (kept, tiny)."""
import numpy as np
import importlib.util

spec = importlib.util.spec_from_file_location("q08", r"D:\Dev\math\w7\q08_largea.py")
# do not import (it would rerun the whole report); copy the needed pieces instead
TWOPI = 2.0 * np.pi
_GL = {}


def gl(n):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = (0.5 * (x + 1.0), 0.5 * w)
    return _GL[n]


def integ(f, lo, hi, n=3000):
    t, wt = gl(n)
    return (hi - lo) * float(np.dot(wt, f(lo + (hi - lo) * t)))


def W(r):
    return 4.0 * np.sin(0.5 * r) ** 2 / r


def PH(r):
    return np.sin(0.5 * r) / r


def PHp(r):
    return (0.5 * r * np.cos(0.5 * r) - np.sin(0.5 * r)) / r ** 2


def Wp(r):
    return (2.0 * r * np.sin(r) - 4.0 * np.sin(0.5 * r) ** 2) / r ** 2


def I_of(a, b, n=3000):
    return integ(lambda r: W(r) * np.exp(-b * a * PH(r)), abs(a), TWOPI, n)


def mom(a, b, n=3000):
    e = lambda r: W(r) * np.exp(-b * a * PH(r))
    I = integ(e, abs(a), TWOPI, n)
    return I, integ(lambda r: e(r) * PH(r), abs(a), TWOPI, n) / I, \
        integ(lambda r: e(r) * PH(r) ** 2, abs(a), TWOPI, n) / I


def G_of(a, b, n=3000):
    return b - np.log(I_of(a, b, n))


h = 1e-3
for a in (0.5, 1.5, np.pi, 4.0, 5.5):
    b = 2.5
    Gp = (-G_of(a + 2 * h, b) + 8 * G_of(a + h, b) - 8 * G_of(a - h, b) + G_of(a - 2 * h, b)) / (12 * h)
    I, m1, m2 = mom(a, b)
    kap = W(a) * np.exp(-b * a * PH(a)) / I
    print(f"a={a:6.3f} I={I:.8f} m1={m1:.8f} var={m2-m1*m1:.8f} kap={kap:.8f} "
          f"Gp_st={Gp:+.8f} kap+b*m1={kap+b*m1:+.8f}  w={W(a):.6f} w'={Wp(a):.6f}")
