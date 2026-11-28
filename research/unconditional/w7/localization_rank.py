import numpy as np
from scipy.special import iv

# Mechanism check (round 22): as lambda grows, the weight w = e^{lam cos} localizes near theta=0,
# where  K(th,ph) = (1/(2b))(e^{b cos(th-ph)} - e^{b cos(th+ph)})  ~  e^b sin(th) sin(ph)
# i.e. K is asymptotically RANK-1 on the localized region.
# Quantify: participation ratio  PR = (sum t_i)^2 / sum t_i^2   of spec(T),
# and the leading-eigenvector concentration.

beta = 2.5
N = 1000
th = (np.arange(N) + 0.5) * np.pi / N
dth = np.pi / N
mu = (2 / np.pi) * dth

def K_closed(a, b, beta):
    return (np.exp(beta * np.cos(a - b)) - np.exp(beta * np.cos(a + b))) / (2 * beta)

TH, PH = np.meshgrid(th, th, indexing='ij')
Kg = K_closed(TH, PH, beta)
A0 = (Kg * mu)
A0 = (A0 + A0.T) / 2

print("=== rank-1 asymptotics on the localized region ===")
for eps in (0.3, 0.2, 0.1, 0.05):
    m = th < eps
    if m.sum() < 3: continue
    Kloc = Kg[np.ix_(m, m)]
    s = np.linalg.svd(Kloc, compute_uv=False)
    s2 = s ** 2
    print("  eps=%.2f : rank-1 fraction  s1^2/sum s^2 = %.10f  (top eigenvalues %s)"
          % (eps, s2[0] / s2.sum(), np.array2string(s[:3], precision=6)))
print("  (K(th,ph) ~ e^beta sin(th) sin(ph) for small th,ph => rank 1)")

print()
print("=== participation ratio of spec(T) vs lambda ===")
def spec(beta, lam):
    w = np.exp(lam * np.cos(th))
    T = np.diag(w) @ A0 @ np.diag(w)
    ev = np.sort(np.linalg.eigvalsh((T + T.T) / 2))[::-1]
    return ev

for lam in (0.0, 0.25, 0.5, 1.0, 1.5, 2.5, 4.0, 6.0):
    ev = spec(beta, lam)
    pos = np.maximum(ev, 0)
    PR = pos.sum() ** 2 / (pos ** 2).sum()
    # concentration of the top eigenvector
    w = np.exp(lam * np.cos(th))
    T = np.diag(w) @ A0 @ np.diag(w)
    evs, U = np.linalg.eigh((T + T.T) / 2)
    x1 = U[:, -1]
    frac_low = (2 / np.pi) * dth * np.sum(x1[th < 1.0] ** 2)
    print("  lam=%.2f : r=t2/t1=%.8f  PR=%.4f  mass(x1 on th<1)=%.4f"
          % (lam, ev[1] / ev[0], PR, frac_low))

print()
print("=== check the rank-1 prediction: is t1 ~ e^{2lam} * const ? ===")
for lam in (0.5, 1.0, 2.0, 3.0, 4.0):
    ev = spec(beta, lam)
    print("  lam=%.1f : t1=%.6e   t1/e^{2lam}=%.6e" % (lam, ev[0], ev[0] / np.exp(2 * lam)))
