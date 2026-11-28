import numpy as np
from scipy.special import iv

# Direct comparison of the two M_f constructions used in rounds 13 and 14.
# Round 13: loop over doubled indices a,c (0..J); d doubled, step 2; c_d = e^-b(I_d - I_{d+2})
# Round 14: loop over doubled indices; d in j-units step 1; c_k = e^-b(I_{2k} - I_{2k+2})
# Both should equal  (M_f)_{j'j} = sum_{j''=|j-j'|}^{j+j'} c_{j''}

def M13(J, b):
    M = np.zeros((J + 1, J + 1))
    for a in range(J + 1):
        for c in range(J + 1):
            s = 0.0
            for d in range(abs(a - c), a + c + 1):
                if (d - abs(a - c)) % 2 == 0:
                    s += np.exp(-b) * (iv(d, b) - iv(d + 2, b))
            M[a, c] = s
    return M

def M14(J, b):
    M = np.zeros((J + 1, J + 1))
    c = np.array([np.exp(-b) * (iv(2 * k, b) - iv(2 * k + 2, b)) for k in range(J + 1)])
    for a in range(J + 1):
        for bb in range(J + 1):
            ja, jb = a / 2.0, bb / 2.0
            lo, hi = abs(ja - jb), ja + jb
            s = 0.0; d = lo
            while d <= hi + 1e-12:
                k = int(round(d))
                if 0 <= k < len(c):
                    s += c[k]
                d += 1.0
            M[a, bb] = s
    return M

J = 10; b = 1.5
A = M13(J, b); B = M14(J, b)
print("=== round-13 vs round-14 M_f ===")
print("  max |A-B| = %.3e" % np.abs(A - B).max())
print("  A[0,:4] =", A[0, :4], " B[0,:4] =", B[0, :4])
print("  A[1,:4] =", A[1, :4], " B[1,:4] =", B[1, :4])
print()

def min_minor(M, J):
    mn = np.inf; arg = None
    for i in range(J):
        for j in range(i + 1, J + 1):
            for k in range(J):
                for l in range(k + 1, J + 1):
                    dt = M[i, k] * M[j, l] - M[i, l] * M[j, k]
                    if dt < mn: mn = dt; arg = (i, j, k, l)
    return mn, arg

for name, M in (("M13", A), ("M14", B)):
    mn, arg = min_minor(M, J)
    print("  %s : min 2x2 minor = %+.6e at %s" % (name, mn, arg))

print()
print("=== which is the correct (M_f)_{j'j}?  direct numerical check via Haar integral ===")
# (M_f)_{j'j} = int chi_{j'} f chi_j dU ; compute by quadrature for a few pairs
def charval(j, om):
    s = np.sin(om); num = np.sin((2 * j + 1) * om)
    return np.where(np.abs(s) < 1e-12, (2 * j + 1), num / np.where(np.abs(s) < 1e-12, 1.0, s))
om = np.linspace(0, np.pi, 400001); meas = (2 / np.pi) * np.sin(om) ** 2
f = np.exp(-b * (1 - np.cos(om)))
def entry(jp, j):
    return np.trapezoid(charval(jp, om) * f * charval(j, om) * meas, om)
for (jp, j) in [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (1.0, 0.0), (1.0, 0.5)]:
    a = int(round(2 * jp)); c = int(round(2 * j))
    print("  (j'=%3.1f, j=%3.1f) : quad=%.10f  M13=%.10f  M14=%.10f"
          % (jp, j, entry(jp, j), A[a, c], B[a, c]))
