import numpy as np
from scipy.special import iv

# Question: is the "peaked positive weight decreases the spectral ratio" phenomenon
# GENERAL for positive kernels, or special to our K?
#   Given a positive-definite kernel operator  (Kf)(x) = int K(x,y) f(y) mu(dy)  with K >= 0,
#   and a positive weight w peaked at one point, compare
#       spec ratio of  K_op        vs   M_w K_op M_w.
# Test on random positive kernels on a discretized interval with the (2/pi)dtheta measure.

rng = np.random.default_rng(7)
N = 60
th = (np.arange(N) + 0.5) * np.pi / N
dth = np.pi / N
mu = (2 / np.pi) * dth          # weight for sums

def spec_ratio(A):
    w = np.sort(np.linalg.eigvalsh((A + A.T) / 2))[::-1]
    if w[0] <= 0: return np.nan
    return w[1] / w[0]

print("=== random positive kernels: does a peaked w decrease t2/t1 ? ===")
fails = 0; tot = 0
for trial in range(40):
    # random positive kernel (entrywise positive, PSD)
    B = rng.uniform(0.05, 1.0, (N, N))
    B = (B + B.T) / 2
    Kmat = B @ B.T / N                       # PSD, entrywise positive
    # peaked positive weight
    lam = rng.uniform(0.5, 3.0)
    w = np.exp(lam * np.cos(th))
    W = np.diag(w)
    A0 = Kmat * mu
    A1 = W @ A0 @ W
    r0, r1 = spec_ratio(A0), spec_ratio(A1)
    tot += 1
    if not (r1 <= r0 + 1e-10):
        fails += 1
print("  tested %d random positive kernels ; 'ratio decreased' FAILURES = %d" % (tot, fails))

print()
print("=== the actual K: closed form check and ratio monotonicity in lambda ===")
beta = 2.5
def K_closed(a, b):
    return (np.exp(beta * np.cos(a - b)) - np.exp(beta * np.cos(a + b))) / (2 * beta)
TH, PH = np.meshgrid(th, th, indexing='ij')
Kg = K_closed(TH, PH)
print("  K closed form: min=%.3e max=%.3e  (K>=0 since cos(a-b)>=cos(a+b))" % (Kg.min(), Kg.max()))
print("  symmetry K(a,b)=K(b,a): %s" % (np.abs(Kg - Kg.T).max() < 1e-12))
print()
for lam in (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0):
    w = np.exp(lam * np.cos(th))
    A0 = Kg * mu
    A1 = np.diag(w) @ A0 @ np.diag(w)
    print("  lam=%.2f : r(K_op)=%.10f   r(M_w K M_w)=%.10f   decreased=%s"
          % (lam, spec_ratio(A0), spec_ratio(A1), spec_ratio(A1) <= spec_ratio(A0) + 1e-12))
