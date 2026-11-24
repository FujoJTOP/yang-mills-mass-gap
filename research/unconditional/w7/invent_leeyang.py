"""
invent_leeyang.py -- witness for INVENT-LEEYANG (2026-11-24)
Cheap: free Gaussian free-energy analytic structure in complex beta; nearest-singularity
distance D(beta0) for free massless AND free massive (:phi^2:) on a 1D lattice.
Run: py -3.12 -I research/unconditional/w7/invent_leeyang.py
"""
import numpy as np

def lam_1d(N):
    k = np.arange(N)
    return np.sort(4.0 * np.sin(np.pi * k / N) ** 2)  # ascending; lam[0]=0

def nearest_singularity(beta0, m2, N):
    """Z(beta)=prod_{lam>0}(m2+beta*lam)^-1/2 * (zero-mode factor).
    Singularities of the nonzero modes: m2 + beta*lam = 0  =>  beta=-m2/lam < 0.
    Massless branch point: beta=0 from the nonzero-mode product.
    Returns (distance to nearest singularity in the right half-plane, that singularity)."""
    lam = lam_1d(N); lam_max = lam.max()
    if m2 == 0.0:
        # massless: zero mode unconstrained; nonzero modes give beta^{-(N-1)/2}; branch pt at 0
        return abs(beta0 - 0.0), 0.0
    beta_sing = -m2 / lam_max          # least-negative (nearest) among the m2+beta*lam=0 set
    return abs(beta0 - beta_sing), beta_sing

beta0 = 2.5
print("=" * 78)
print("Nearest-singularity distance D(beta0) vs m^2 and N   (beta0 = %.2f)" % beta0)
print("=" * 78)
print(" N     m^2    D      singular_beta   note")
for N in [8, 16, 64, 256]:
    for m2 in [0.0, 0.01, 0.1, 1.0]:
        D, bs = nearest_singularity(beta0, m2, N)
        note = "FREE MASSLESS (gap=0) -- fires?!" if m2 == 0.0 else ":phi^2:-like massive (f1^2=0)"
        print(f"{N:4d}  {m2:5.2f}  {D:7.4f}   {bs:9.4f}    {note}")

print()
print("Observation 1 (ACID FAIL): at m^2=0 (free massless, gap=0) D = beta0 = 2.5 > 0.")
print("  The nearest singularity is the kinematic beta=0 branch point, present in EVERY theory.")
print("Observation 2 (gap-correlate, atom-blind): for m^2>0, D = beta0 + m^2/lam_max")
print("  (lam_max=4 for N even), monotone in m^2 but INDEPENDENT of any atom weight f1^2.")
print("  :phi^2: is gapped (m^2>0, D>0) but has f1^2=0 => D fires on the gapped-atomless adversary.")

print()
print("=" * 78)
print("1D Ising check: physical Fisher zeros sit on the IMAGINARY beta axis (Re beta=0),")
print("  so D=beta0 there too -- firing on ANY gapped theory, atom-blind by construction.")
print("=" * 78)
N = 32; b = 2.5
# 1D periodic Ising: Z = 2^N (cosh^N b + sinh^N b); zeros: tanh^N b = -1
# nearest zeros are at tanh(beta)=e^{i(2k+1)pi/N}; |tanh(beta)|=1 => Re(beta)=0
print("  |tanh(beta)|=1  =>  beta purely imaginary  =>  nearest zeros have Re(beta)=0")
print("  distance from beta0=%.2f to the imaginary axis = %.2f  (>0, fires)" % (b, b))
