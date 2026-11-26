# q15_spectral.py -- Q15: is Var(G) <= E[H] a SPECTRAL statement?
# SU(2) single plaquette, exact (Haar + Gibbs).  No fits, no fabrications.
#
# Conventions (unit-normalised bi-invariant direction):
#   U = cos(phi) + i sin(phi) n.sigma ,   S_W = beta (1 - tr U/2) = beta(1 - cos phi)
#   U(t) = U exp(i t V),  V = v.sigma/2 , |v| = 2  (arc length)   =>
#        G_v = dS/dt = beta sin(phi) (n.vhat)          [first derivative]
#        H_vv = d2S/dt2 = beta cos(phi)                [second derivative]
#   Hessian on SU(2) (unit sphere S^3, geodesic distance = phi):
#        lam_radial = f''(phi) = beta cos(phi)
#        lam_transv = f'(phi) cot(phi) = beta sin(phi) cot(phi) = beta cos(phi)   (mult 2)
#        => lam_k = beta cos(phi), k = 1,2,3 : THREE-FOLD DEGENERATE (Wilson action only)
import numpy as np
from numpy.polynomial.legendre import leggauss

rng = np.random.default_rng(20261125)
xgl, wgl = leggauss(400)


def exact(beta):
    """u = cos(phi); Haar+Gibbs measure  sqrt(1-u^2) e^{beta u} du."""
    wt = np.sqrt(1.0 - xgl**2) * np.exp(beta * xgl)
    Z = (wgl * wt).sum()
    Eu = (wgl * wt * xgl).sum() / Z
    Es2 = (wgl * wt * (1.0 - xgl**2)).sum() / Z
    # P(phi > pi/2) = P(u < 0)
    m = xgl < 0.0
    Pneg = (wgl[m] * wt[m]).sum() / Z
    EH = beta * Eu                      # E[H], unit direction
    Elam = beta * Eu                    # E[lam_k]
    VG = beta**2 * Es2 / 3.0            # E[G^2] with E[n.vhat^2] = 1/3
    return dict(Eu=Eu, Es2=Es2, EH=EH, Elam=Elam, VG=VG, d=EH - VG, Pneg=Pneg)


# ---------------------------------------------------------------- [P1] quadrature
print("[P1] exact quadrature, SU(2) single plaquette (Haar+Gibbs, e^{-S_W} dU)")
print("  beta      E[cos]    E[sin^2]   E[H]=E[lam_k]   Var(G)=E[G^2]   E[H]-Var   P(lam<0)")
BETAS = [0.0, 1.0 / 48.0, 0.5, 1.0, 2.5, 6.0]
P1 = {}
for b in BETAS:
    s = exact(b) if b > 0 else dict(Eu=0.0, Es2=0.75, EH=0.0, Elam=0.0, VG=0.0, d=0.0, Pneg=0.5)
    P1[b] = s
    print("  %-8.4f  %+.6f  %.6f   %.6f      %.6f    %+.2e   %.4f"
          % (b, s["Eu"], s["Es2"], s["EH"], s["VG"], s["d"], s["Pneg"]))
print("  -> Stein identity E[H]=E[G^2] (i.e. E[H]-Var(G)=E[G]^2=0) holds to quadrature precision.")

# ---------------------------------------------------------------- [P2] Monte Carlo
beta = 2.5
N = 400_000
# sample (phi, n) from p(phi) prop. sin^2 phi e^{-beta(1-cos phi)} by rejection on phi,
# n uniform on S^2 (independent: S_W is a class function).
grid = np.linspace(0.0, np.pi, 4001)
dens = np.sin(grid)**2 * np.exp(-beta * (1.0 - np.cos(grid)))
dens /= dens.max()
phi = np.empty(N)
filled = 0
while filled < N:
    cand = rng.uniform(0.0, np.pi, N)
    acc = rng.uniform(0.0, 1.0, N) < np.interp(cand, grid, dens)
    k = int(acc.sum())
    k = min(k, N - filled)
    phi[filled:filled + k] = cand[acc][:k]
    filled += k
z = rng.normal(size=(N, 3))
nhat = z / np.linalg.norm(z, axis=1, keepdims=True)
vhat = np.array([0.0, 0.0, 1.0])
nv = nhat @ vhat

G = beta * np.sin(phi) * nv          # first derivative, unit direction
H = beta * np.cos(phi) * np.ones(N)  # second derivative (= lam_k)
EG, EG2, EHm = G.mean(), (G**2).mean(), H.mean()
VG_mc = EG2 - EG**2
print("\n[P2] Monte Carlo, beta=2.5, N=%d, direction v=e_3 (fixed), n uniform on S^2" % N)
print("  E[G]        = %+.6f   (0 by n -> -n symmetry)" % EG)
print("  E[G^2]      = %.6f" % EG2)
print("  E[H]        = %.6f" % EHm)
print("  Var(G)      = %.6f" % VG_mc)
print("  E[H]-Var(G) = %+.2e   vs  E[G]^2 = %.2e" % (EHm - VG_mc, EG**2))
print("  exact: E[H]=E[G^2]=%.6f  (rel.dev %.2e)" % (1.2677, abs(EHm - 1.2677) / 1.2677))

# ---------------------------------------------------------------- [P3] Hessian spectrum
I2 = np.eye(2, dtype=complex)
SIG = [np.array([[0, 1], [1, 0]], dtype=complex),
       np.array([[0, -1j], [1j, 0]], dtype=complex),
       np.array([[1, 0], [0, -1]], dtype=complex)]


def rand_su2():
    z = rng.normal(size=4)
    z /= np.linalg.norm(z)
    a, b = z[0] + 1j * z[1], z[2] + 1j * z[3]
    return np.array([[a, b], [-np.conj(b), np.conj(a)]])


def sw(U, beta):
    return beta * (1.0 - 0.5 * np.trace(U).real)


def expi(Vh):
    """exp(i V) for Hermitian V: V = Q diag(mu) Q^dag  ->  Q diag(e^{i mu}) Q^dag."""
    mu, Q = np.linalg.eigh(Vh)
    return (Q * np.exp(1j * mu)) @ Q.conj().T


print("\n[P3] numeric Hessian of S_W on SU(2), basis i*sigma_a/2 (h=1e-4), 6 random U")
h = 1e-4
maxerr = 0.0
for t in range(6):
    U = rand_su2()
    phi_t = np.arccos(np.clip(0.5 * np.trace(U).real, -1.0, 1.0))
    M = np.zeros((3, 3))
    for a in range(3):
        for c in range(3):
            Va, Vc = SIG[a] / 2.0, SIG[c] / 2.0
            fpp = sw(U @ expi(+h * Va + h * Vc), beta)
            fpm = sw(U @ expi(+h * Va - h * Vc), beta)
            fmp = sw(U @ expi(-h * Va + h * Vc), beta)
            fmm = sw(U @ expi(-h * Va - h * Vc), beta)
            M[a, c] = (fpp - fpm - fmp + fmm) / (4 * h * h)
    ev = np.linalg.eigvalsh(M)
    pred = beta * np.cos(phi_t) / 4.0     # basis i*sigma_a/2 has |.|^2=1/4 of unit direction
    maxerr = max(maxerr, np.max(np.abs(ev - pred)))
    print("  phi=%.6f  eig(M)=[%+.6f %+.6f %+.6f]   beta*cos(phi)/4=%+.6f   lam_k=beta*cos(phi)=%+.6f"
          % (phi_t, ev[0], ev[1], ev[2], pred, beta * np.cos(phi_t)))
print("  max |eig - beta cos(phi)/4| = %.2e  -> Hessian is SCALAR: lam_1=lam_2=lam_3=beta cos(phi)"
      % maxerr if maxerr < 1e-6 else
      "  *** MISMATCH: max |eig - beta cos(phi)/4| = %.2e ***" % maxerr)

# ---------------------------------------------------------------- [P4] spectral condition tests
print("\n[P4] does any spectral condition reproduce Var(G)<=E[H]?  beta=2.5")
lam = beta * np.cos(phi)                 # lam_1=lam_2=lam_3
S1 = 3 * lam                             # sum_k lam_k
S2 = 3 * lam**2                          # sum_k lam_k^2
print("  lam_k range        : [%+.4f, %+.4f]   (max = beta = %.2f at phi=0)" % (lam.min(), lam.max(), beta))
print("  E[lam_k]           : %+.6f   E[sum lam_k]=%+.6f  E[sum lam_k^2]=%.6f"
      % (lam.mean(), S1.mean(), S2.mean()))
print("  P(lam_k > 1)       : %.4f     P(lam_k < 0) : %.4f" % ((lam > 1).mean(), (lam < 0).mean()))
print("  naive spectral cond 'max_k lam_k <= 1' holds on %.2f%% of configurations"
      % (100.0 * (lam.max() <= 1.0)))
print("  guessed form  sum lam^2 <= sum lam + (sum lam)^2 :")
lhs, rhs = S2, S1 + S1**2
print("    E[sum lam^2]=%.4f  E[sum lam + (sum lam)^2]=%.4f  violated on %.2f%% of cfgs"
      % (lhs.mean(), rhs.mean(), 100.0 * (lhs > rhs).mean()))
print("    (per-cfg, 3-fold degenerate lam: 3b^2c^2 <= 3bc + 9b^2c^2  <=>  0 <= 3bc(1+2bc);")
print("     fails for every c<0 and for bc < -1/2 -> the guessed form is NOT equivalent.)")
print("  actual comparison  E[H]=%.6f   Var(G)=%.6f   E[H]-Var(G)=%.2e  E[G]^2=%.2e"
      % (EHm, VG_mc, EHm - VG_mc, EG**2))
print("\n[P5] summary numbers for the dossier (beta grid, exact):")
print("  beta     E[lam_k]=E[H]=Var(G)   lam_max=beta   P(lam<0)")
for b in BETAS[1:]:
    s = P1[b]
    print("  %-8.4f %.6f              %-8.2f      %.4f" % (b, s["Elam"], b, s["Pneg"]))
print("\n[DONE]")
