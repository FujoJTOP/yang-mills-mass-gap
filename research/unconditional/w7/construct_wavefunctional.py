"""
construct_wavefunctional.py -- witness for
research/unconditional/CONSTRUCT-WAVEFUNCTIONAL.md  (2026-11-24)

QUESTION: does a CONTROLLED-VARIANCE 0++-type wavefunctional close the Temple
sandwich?  Calibrate on a genuinely interacting, solvable model (1D anharmonic
oscillator H = -1/2 d^2/dx^2 + 1/2 w^2 x^2 + lam x^4, w=1), then read 3+1D.

Temple (1928): for normalized trial psi, rho=<H>, sigma^2=<H^2>-rho^2,
E0 <= rho < E1, and W >= E1 an upper bound to E1:
        E1 >= rho - sigma^2/(W - rho).
It certifies a GAP (E1 >= E0^L) iff the Temple ratio
        sigma^2 / ((W - rho)(rho - E0^L))  < 1.

Exact spectrum : dense ED in the harmonic-oscillator basis (no factorial).
Trial functions: spectral (FFT) position space -- smooth, no overflow.

Run:  py -3.12 -I research/unconditional/w7/construct_wavefunctional.py
"""
import numpy as np


# ---------------------------------------------------------------- exact spectrum
def exact_spectrum(lam, w=1.0, N=200):
    n = np.arange(N)
    a = np.zeros((N, N))
    for i in range(N - 1):
        a[i, i + 1] = np.sqrt(i + 1)          # a_{n,n+1}=sqrt(n+1)
    x = (a + a.T) / np.sqrt(2 * w)
    x4 = (x @ x) @ (x @ x)
    H = np.diag(w * (n + 0.5)) + lam * x4
    ev, U = np.linalg.eigh(H)
    return H, ev, U


# ---------------------------------------------------------------- spectral grid
class Grid:
    def __init__(self, lam, w=1.0, Lg=12.0, Ng=2048):
        self.w, self.lam = w, lam
        self.x = np.linspace(-Lg, Lg, Ng, endpoint=False)
        self.dx = 2 * Lg / Ng
        self.k = 2 * np.pi * np.fft.fftfreq(Ng, self.dx)
        self.V = 0.5 * w * w * self.x ** 2 + lam * self.x ** 4

    def H(self, psi):
        d2 = np.fft.ifft(-(self.k ** 2) * np.fft.fft(psi)).real
        return -0.5 * d2 + self.V * psi

    def rv(self, psi):
        psi = psi / np.sqrt(np.sum(psi * psi) * self.dx)
        Hp = self.H(psi)
        rho = np.sum(psi * Hp) * self.dx
        h2 = np.sum(Hp * Hp) * self.dx
        return rho, max(h2 - rho * rho, 0.0)


def rv_basis(H, c):
    c = c / np.linalg.norm(c)
    Hc = H @ c
    rho = c @ Hc
    return rho, max(Hc @ Hc - rho ** 2, 0.0)


# ---------------------------------------------------------------- analysis
def analyse(lam, label):
    H, ev, U = exact_spectrum(lam)
    E0, E1, E2 = ev[0], ev[1], ev[2]
    gap = E1 - E0
    g = Grid(lam)
    print("=" * 76)
    print(f"{label}:  lam = {lam}")
    print(f"  exact  E0={E0:.6f}  E1={E1:.6f}  E2={E2:.6f}  gap={gap:.6f}")

    sgrid = np.linspace(0.45, 1.9, 30)

    # generic Gaussian families
    ev_e = [g.rv(np.exp(-g.x ** 2 / (2 * s ** 2))) for s in sgrid]     # even
    ev_o = [g.rv(g.x * np.exp(-g.x ** 2 / (2 * s ** 2))) for s in sgrid]  # odd
    rho_e_min = min(e[0] for e in ev_e)
    var_e_min = min(e[1] for e in ev_e)
    W = min(o[0] for o in ev_o)          # upper bound on E1 (odd trial >= E1)
    E0L = E0                              # calibration: exact lower bound
    print(f"  generic even-Gaussian: min<H>={rho_e_min:.6f} (var floor={var_e_min:.3e})"
          f"   |  odd-Gaussian min<H>={W:.6f}  => W(=E1^upper)")
    print(f"  E0^L = {E0L:.6f}")

    # Temple scan over mixed trials  a*even(s1) + b*odd(s2)
    best = None
    s_scan = np.linspace(0.5, 1.8, 14)
    ts = np.linspace(0.03, 0.6, 15)
    for s1 in s_scan:
        fe = np.exp(-g.x ** 2 / (2 * s1 ** 2))
        for s2 in s_scan:
            fo = g.x * np.exp(-g.x ** 2 / (2 * s2 ** 2))
            for t in ts:
                psi = np.cos(t) * fe + np.sin(t) * fo
                rho, var = g.rv(psi)
                if not (E0L < rho < W):
                    continue
                tb = rho - var / (W - rho)
                ratio = var / ((W - rho) * (rho - E0L))
                if best is None or tb > best[0]:
                    best = (tb, ratio, rho, var, t, s1, s2)
    if best:
        tb, ratio, rho, var, t, s1, s2 = best
        print(f"  BEST generic mixed: rho={rho:.5f} var={var:.3e} (t={t:.2f})")
        print(f"     Temple E1>={tb:.6f}  (gap above E0 = {tb-E0L:.3e})  ratio={ratio:.3f}"
              f"  [{'CERTIFIES' if tb>E0L else 'NO'} gap]")
    else:
        print("  no generic mixed trial fell in E0<rho<E1 -> generic families miss the Temple regime")

    # controlled variance: perturb around the EXACT eigenvector (HO basis, no grid)
    print("  controlled variance  psi = normalize(U[:,1]+eps*U[:,k]):")
    for k in (0, 2, 3):
        for eps in (0.3, 0.1, 0.03):
            c = U[:, 1].copy()
            c[k] += eps
            rho, var = rv_basis(H, c)
            if not (E0L < rho < W):
                continue
            tb = rho - var / (W - rho)
            ratio = var / ((W - rho) * (rho - E0L))
            print(f"     k={k} eps={eps:5.2f}: var={var:.2e}  Temple={tb:.6f}"
                  f"  gap={tb-E0L:.2e}  ratio={ratio:.4f}")
    return dict(E0=E0, E1=E1, gap=gap, W=W, var_floor=var_e_min, best=best)


r0 = analyse(0.0, "FREE oscillator (acid baseline)")
r1 = analyse(0.5, "INTERACTING anharmonic oscillator")

print("=" * 76)
print("SUMMARY")
print(f"  FREE (lam=0):     generic-Gaussian residual var floor = {r0['var_floor']:.2e}")
print("                    (==0: exact ground IS Gaussian -> controllable variance is trivial/vacuous)")
print(f"  INTERACTING (0.5):generic-Gaussian residual var floor = {r1['var_floor']:.2e}")
print("                    (>0: no Gaussian is an eigenstate -> a FLOOR on any Gaussian trial's variance)")
if r1["best"]:
    tb, ratio = r1["best"][0], r1["best"][1]
    print(f"  interacting best generic Temple deficit to E1 = {r1['E1']-tb:.2e}, ratio = {ratio:.3f}")
print("  => Temple closes the sandwich ONLY as var->0, i.e. as the trial -> the exact eigenstate.")
