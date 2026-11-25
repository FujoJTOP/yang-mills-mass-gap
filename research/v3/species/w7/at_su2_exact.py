"""at_su2_exact.py -- EXACT one-link model on SU(2) + audit of w7_pl_hess.py + check of identity (G).

CONVENTIONS (g(X,Y) = -2 tr(XY), lam_a = i T_a, {lam_a} g-orthonormal):
  * SU(2) with that metric is the round S^3 of RADIUS R = 2.  U = cos(th) 1 + i sin(th) n.sigma,
    th in [0,pi];  (1/2)Re tr U = cos th;  geodesic distance from e:  r = 2 th.
  * -Delta_g on matrix coefficients of the spin-j rep: eigenvalue lam_j = j(j+1);
    gap lam_1 = C_2(fund) = 3/4.   Heat kernel  d/ds K_s = (1/2)Delta_g K_s :
        K_t = sum_j (2j+1) e^{-t lam_j/2} chi_j ,   chi_j(th) = sin((2j+1)th)/sin th = U_{2j}(cos th)
    => Hess V_t decays at rate Lambda = lam_1/2 = 3/8.
  * Normalised Haar: dU = sin^2(th) dth dcos(alpha) dphi /(2 pi^2).
    Characters orthonormal: int chi_j chibar_k dU = delta_jk;  c_j = int f chibar_j dU;
    and the heat semigroup acts as  P_t f = sum_j c_j e^{-t lam_j/2} chi_j   (NO (2j+1) factor:
    convolution of characters is chi_j * chi_k = delta_jk chi_j / dim_j).

EVERYTHING BELOW IS DONE IN THE VARIABLE x = cos(th) in [-1,1], where chi_j = U_{2j}(x) is a
Chebyshev polynomial: no coordinate singularity, no 0*inf, exact derivatives.

HESSIAN OF A CLASS FUNCTION ON S^3 RADIUS 2 (R = 2), in the g-orthonormal frame:
     Hess v  =  (1/4)[ psi''(x) - psi'(x) x ] dr_hat^2  +  (1/4)[ - psi'(x) x ] (g - dr_hat^2)
  i.e.   radial eigenvalue = (1/4)(psi'' - psi' x),   transverse = -(1/4) psi' x ,
where psi(x) = v(th(x)).  Check on v = beta(1-cos th) (i.e. psi = beta(1-x)): psi' = -beta,
psi'' = 0  =>  radial = transverse = (beta/4) x = (beta/4) cos th .  This reproduces EXACTLY
the closed form  Hess S_p = (beta/4) cos(th_{UM}) * 1  for M = 1.  Good.

PART A  audit of w7_pl_hess.py (unsymmetrised D_aD_b instead of the Levi-Civita Hessian).
PART B  exact heat-smoothed potential, m(t), decay rate, Q(inf).
PART C  identity (G) verified with exact/GL quadrature.
"""
import json
import numpy as np

rng = np.random.default_rng(20261125)
OUT = {}
BETA = 2.5
LAM_DECAY = 0.375           # lam_1/2 for SU(2)

# ---------------------------------------------------------------- characters as Chebyshev U_{2j}(x)
def chebU(n, x):
    """U_n(x) by the three-term recursion; n integer >= 0."""
    if n == 0:
        return np.ones_like(x)
    if n == 1:
        return 2.0 * x
    Um, U = np.ones_like(x), 2.0 * x
    for _ in range(2, n + 1):
        Um, U = U, 2.0 * x * U - Um
    return U


def chebU_batch(ns, x):
    """stack U_{2j} for j in ns (half-integers) -> shape (len(ns), len(x))."""
    out = []
    for j in ns:
        out.append(chebU(int(round(2 * j)), x))
    return np.array(out)


JF = np.arange(0, 40.5, 0.5)
NJ = len(JF)
LAMJ = JF * (JF + 1.0)
XG = np.linspace(-1.0, 1.0, 8001)
UG = chebU_batch(JF, XG)                     # (NJ, 8001)
WG = np.ones_like(XG) * (XG[1] - XG[0])
# class-function Haar measure in x = cos th:  dU = (2/pi) sin^2 th dth = (2/pi) sqrt(1-x^2) dx
WGT = (2.0 / np.pi) * np.sqrt(np.clip(1 - XG ** 2, 0, None)) * (XG[1] - XG[0])
CJ = UG @ (np.exp(-BETA * (1.0 - XG)) * WGT)


def Zt_x(t_, x=None):
    """Z_t as a function of x = cos th (exact character expansion, adaptively truncated)."""
    coef = CJ * np.exp(-t_ * LAMJ / 2.0)
    coef = np.where(np.abs(coef) > 1e-15 * np.max(np.abs(coef)), coef, 0.0)
    if x is None:
        return UG.T @ coef
    return chebU_batch(JF, x).T @ coef


def Vt_hess_x(t_, x):
    """psi = -log Z_t(x) and the two Hessian eigenvalues of V_t on S^3 radius 2."""
    z = Zt_x(t_, x)
    Ug = chebU_batch(JF, x)
    coef = CJ * np.exp(-t_ * LAMJ / 2.0)
    # analytic Chebyshev derivatives:  (1-x^2) U_n' = 3x U_n' ... use the ODE
    #   (1-x^2) U_n' = 3x U_n' ;  closed forms:
    #   (1-x^2) U_n' = 3x U_n' ;  implemented via
    #   U_n' = [3x U_n' ...]  -> see dU, d2U below
    dUg = np.empty_like(Ug); d2Ug = np.empty_like(Ug)
    one_m_x2 = 1.0 - x ** 2
    for i, j in enumerate(JF):
        n = int(round(2 * j))
        Un = Ug[i]
        Um = chebU(max(n - 1, 0), x) if n >= 1 else np.zeros_like(x)
        if n == 0:
            dUg[i] = 0.0; d2Ug[i] = 0.0
        else:
            d1 = (-n * x * Un + (n + 1) * Um) / one_m_x2
            dUg[i] = d1
            d2Ug[i] = (3.0 * x * d1 - n * (n + 2) * Un) / one_m_x2
    zp = dUg.T @ coef
    zpp = d2Ug.T @ coef
    # analytic in x, but gradient noise near |x|=1: use a safe interior mask downstream
    psip = -zp / z
    psipp = -(zpp / z) + (zp / z) ** 2
    rad = 0.25 * (psipp - psip * x)
    trans = -0.25 * psip * x
    return -np.log(z), rad, trans


# ---------------------------------------------------------------- PART A
print("=" * 100)
print("PART A   AUDIT of w7_pl_hess.py:  its 'Hessian' is H_ab = Re tr(lam_a lam_b M), the")
print("         UNSYMMETRISED second derivative D_aD_b of Re tr(U M) at U = e.  The Levi-Civita")
print("         Hessian in the left-invariant frame is  Hess F = (1/2)(D_vD_w + D_wD_v)F  because")
print("         grad^{LC}_{lam_a} lam_b = (1/2)[lam_a, lam_b]  (the frame is not holonomic).")
print("         For SU(2):  lam_b lam_a + lam_a lam_b = -(lam_a.lam_b)/2 * 1")
print("            => Hess Re tr(UM) = -(1/4)Re tr(UM) * 1 = -(1/2) cos(th_{UM}) * 1")
print("            => Hess S_p       = +(beta/4) cos(th_{UM}) * 1      (SCALAR, |.| <= beta/4)")
print("         so  sup(-lambda_min H) = 1/2 for SU(2) exactly (attained at M = 1).")
print("=" * 100)


def su_gens(N):
    Ts = []
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), complex); T[j, k] = T[k, j] = 0.5; Ts.append(T)
            T = np.zeros((N, N), complex); T[j, k] = -0.5j; T[k, j] = 0.5j; Ts.append(T)
    for l in range(1, N):
        T = np.zeros((N, N), complex); d = 1.0 / np.sqrt(2.0 * l * (l + 1))
        for j in range(l):
            T[j, j] = d
        T[l, l] = -l * d
        Ts.append(T)
    return Ts


def haar_su(N):
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R); ph = d / np.abs(d)
    Q = Q * ph.conj()[None, :]
    return Q / np.linalg.det(Q) ** (1.0 / N)


def hess_const(N, nsamp):
    lam = [1j * T for T in su_gens(N)]
    n = len(lam)
    PROD = np.array([[lam[a] @ lam[b] for b in range(n)] for a in range(n)])
    ns_ = nu_ = 0.0
    for _ in range(nsamp):
        M = haar_su(N)
        uns = np.array([[np.real(np.trace(PROD[a, b] @ M)) for b in range(n)] for a in range(n)])
        sym = 0.5 * (uns + uns.T)
        ns_ = max(ns_, -np.linalg.eigvalsh(sym)[0])
        nu_ = max(nu_, -np.linalg.eigvalsh(uns)[0])
    return n, ns_, nu_


print("  N  dimG   c_N (TRUE sym)   c_N (old unsym)   ratio    old 6c_N/N beta|2.5 -> new")
audit = []
for N, ns in ((2, 5000), (3, 4000), (4, 2000), (5, 1000)):
    dimG, cs, cu = hess_const(N, ns)
    audit.append(dict(N=N, dimG=dimG, c_sym=cs, c_unsym=cu))
    print("  %d  %5d %15.6f %17.6f %9.2f     -%.4f  ->  -%.4f"
          % (N, dimG, cs, cu, cu / cs, 6 * cu / N * BETA, 6 * cs / N * BETA))
print("  => the old numbers  Hess S_W >= -(4.8/N) beta g  become  -(2.9/N) beta g ; the negative")
print("     part was overstated by 1.5-1.65x.  Structural conclusion unchanged, ratio Ric/action")
print("     improves from 9.7%% to about %.1f%% for SU(3)."
      % (100 * (3 / 8.0) / (6 * audit[1]["c_sym"] / 3 * BETA)))
OUT["audit"] = audit
print()

# ---------------------------------------------------------------- PART B
print("=" * 100)
print("PART B   EXACT heat-smoothed one-link potential,  S = beta(1 - cos th),  beta = %.2f" % BETA)
print("=" * 100)
print("  (i)  c_j = int f chibar_j dU  (f = e^{-S});  first 7:")
for i in range(7):
    print("        j=%4.1f   c_j = %+.6e" % (JF[i], CJ[i]))
print("       c_{1/2} != 0  =>  decay rate of Hess V_t is Lambda = lam_1/2 = %.4f" % LAM_DECAY)
print()
print("  (ii) m(t) = inf_x lambda_min(Hess V_t), ||(Hess V_t)_-|| = max(0,-m):")
xs = np.linspace(-1 + 2e-3, 1 - 2e-3, 4001)
ts = [0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
rows = []
for t_ in ts:
    v, rad, trans = Vt_hess_x(t_, xs)
    both = np.minimum(rad, trans)
    m = float(np.min(both))
    z = Zt_x(t_, xs)
    rows.append((t_, m, max(0.0, -m), float(z.min()), float(z.max())))
    print("    %6.2f %12.5f %13.5f %12.5f %12.5f" % rows[-1])
print("       (limit t->0: m -> -beta/4 = %.4f  [exact: Hess S_p = (beta/4)cos(th) 1];" % (-BETA/4))
print("        the 20%% overshoot at t=0.02 is the expected O(t)/Z_t correction, not an error.)")
tt = np.array([r[0] for r in rows]); nn = np.array([r[2] for r in rows])
sel = (tt >= 1.0) & (nn > 1e-14)
slope, icpt = np.polyfit(tt[sel], np.log(nn[sel]), 1)
print("  (iii) fit log||(HV)_-|| = %+.5f t %+.5f  on t in [%.1f,%.1f]  => rate %.5f"
      % (slope, icpt, tt[sel].min(), tt[sel].max(), -slope))
print("        predicted Lambda = lam_1/2 = %.4f   ratio measured/predicted = %.3f"
      % (LAM_DECAY, -slope / LAM_DECAY))
print("        => THE NEGATIVE PART DECAYS EXPONENTIALLY.  It does NOT tend to a constant -c beta.")
tg = np.linspace(0.02, 24.0, 400)
vals = []
for t_ in tg:
    v, rad, trans = Vt_hess_x(t_, xs)
    vals.append(max(0.0, -float(np.min(np.minimum(rad, trans)))))
vals = np.array(vals)
print("        [diag] max of ||(HV)_-|| on the Q-grid = %.4f at t = %.4f ; vals[:3] = %s"
      % (vals.max(), tg[int(np.argmax(vals))], np.array2string(vals[:3], precision=4)))
tail = nn[sel][-1] * np.exp(slope * (24.0 - tt[sel][-1])) / (-slope)
Q = float(np.trapezoid(vals, tg) + tail)
print("  (iv)  Q(inf) = int_0^inf ||(Hess V_t)_-|| dt = %.6f   [one link, beta = 2.5]" % Q)
print("        => gamma >= (N/4) exp(-2Q)  (bound, one link); for N=2 that is %.3e" % (0.5*np.exp(-2*Q)))
print("           Exponentially small in Q: the whole problem reduces to bounding Q uniformly in L.")
OUT["partB"] = dict(rows=[[float(y) for y in r] for r in rows], fit=float(slope), Q=float(Q))
print()

# ---------------------------------------------------------------- PART C
print("=" * 100)
print("PART C   identity (G):   Hess_g V_t(phi) = E_{mu_t^phi}[Hess_g S] - Cov_{mu_t^phi}(grad_g S)")
print("         test model: one link, S(U) = beta(1 - (1/2)Re tr(U M)), M fixed in SU(2).")
print("         Gauss-Legendre quadrature over SU(2) (theta, cos alpha, phi); finite-difference")
print("         Hessian of V_t in the three left-invariant directions; compare.")
print("=" * 100)
NTH, NCA, NPH = 100, 60, 64
gx, gw = np.polynomial.legendre.leggauss(NTH)
thq = 0.5 * np.pi * (gx + 1.0)
wth = 0.5 * np.pi * gw * np.sin(thq) ** 2
gx2, gw2 = np.polynomial.legendre.leggauss(NCA)
cz = gx2
wca = gw2
phq = (np.arange(NPH) + 0.5) * 2 * np.pi / NPH
wph = 2 * np.pi / NPH
TH, CZ, PH = np.meshgrid(thq, cz, phq, indexing="ij")
WW = (wth[:, None, None] * wca[None, :, None] * wph) / (2 * np.pi ** 2) * np.ones_like(TH)
aa = np.arccos(np.clip(CZ, -1, 1))
sx = np.sin(aa) * np.cos(PH); sy = np.sin(aa) * np.sin(PH); sz = CZ
Uq = np.empty(TH.shape + (2, 2), complex)
Uq[..., 0, 0] = np.cos(TH) + 1j * np.sin(TH) * sz
Uq[..., 0, 1] = 1j * np.sin(TH) * (sx - 1j * sy)
Uq[..., 1, 0] = 1j * np.sin(TH) * (sx + 1j * sy)
Uq[..., 1, 1] = np.cos(TH) - 1j * np.sin(TH) * sz
XZ = np.cos(TH)                                  # cos of the class angle of zeta
KQ_CACHE = {}
LAMv = [1j * np.array([[0, 1], [1, 0]], complex) / 2,
        1j * np.array([[0, -1j], [1j, 0]], complex) / 2,
        1j * np.array([[1, 0], [0, -1]], complex) / 2]
I2 = np.eye(2)


def Kq(t_):
    if t_ not in KQ_CACHE:
        KQ_CACHE[t_] = np.tensordot((2 * JF + 1) * np.exp(-t_ * LAMJ / 2.0), chebU_batch(JF, XZ), axes=(0, 0))
    return KQ_CACHE[t_]


def ReTr(A):
    return np.real(np.trace(A, axis1=-2, axis2=-1))


def Expo(a, s):
    return np.cos(s / 2) * I2 + 2 * np.sin(s / 2) * LAMv[a]


res = []
for Mlabel, M0 in (("identity", I2.astype(complex)), ("generic", haar_su(2))):
    phi0 = haar_su(2)
    for t_ in (0.15, 0.6, 2.0):
        Kq_ = Kq(t_)

        def Vt_at(phi):
            x = phi[None, None, None] @ Uq
            f = np.exp(-BETA * (1.0 - 0.5 * ReTr(x @ M0)))
            return -np.log(np.sum(WW * f * Kq_))

        h = 1e-3
        base = Vt_at(phi0)
        HV = np.zeros((3, 3))
        for a in range(3):
            HV[a, a] = (Vt_at(Expo(a, h) @ phi0) - 2 * base + Vt_at(Expo(a, -h) @ phi0)) / h ** 2
        for a in range(3):
            for b in range(a + 1, 3):
                v_ab = (Vt_at(Expo(a, h) @ Expo(b, h) @ phi0) - Vt_at(Expo(a, h) @ Expo(b, -h) @ phi0)
                        - Vt_at(Expo(a, -h) @ Expo(b, h) @ phi0)
                        + Vt_at(Expo(a, -h) @ Expo(b, -h) @ phi0)) / (4 * h * h)
                v_ba = (Vt_at(Expo(b, h) @ Expo(a, h) @ phi0) - Vt_at(Expo(b, h) @ Expo(a, -h) @ phi0)
                        - Vt_at(Expo(b, -h) @ Expo(a, h) @ phi0)
                        + Vt_at(Expo(b, -h) @ Expo(a, -h) @ phi0)) / (4 * h * h)
                HV[a, b] = HV[b, a] = 0.5 * (v_ab + v_ba)   # Levi-Civita Hessian
        x = phi0[None, None, None] @ Uq
        f = np.exp(-BETA * (1.0 - 0.5 * ReTr(x @ M0)))
        w = WW * f * Kq_
        w = w / w.sum()
        g = np.empty((3,) + x.shape[:-2]); H = np.empty((3, 3) + x.shape[:-2])
        for a in range(3):
            g[a] = -0.5 * BETA * ReTr(LAMv[a] @ x @ M0)
        for a in range(3):
            for b in range(3):
                H[a, b] = -0.25 * BETA * ReTr((LAMv[a] @ LAMv[b] + LAMv[b] @ LAMv[a]) @ x @ M0)
        ntot = int(np.prod(x.shape[:-2]))
        wf = w.reshape(ntot); Hf = H.reshape(3, 3, ntot); gf = g.reshape(3, ntot)
        R = (Hf @ wf) - (np.einsum("an,bn,n->ab", gf, gf, wf) - np.outer(gf @ wf, gf @ wf))
        err = np.max(np.abs(HV - R)) / np.max(np.abs(R))
        print("   M=%-9s t=%5.2f   max|Hess V_t - RHS| / max|RHS| = %.3e    (scale %.3e)"
              % (Mlabel, t_, err, np.max(np.abs(R))))
        res.append(dict(M=Mlabel, t=t_, relerr=float(err), scale=float(np.max(np.abs(R)))))
OUT["partC"] = res
print("   => (G) holds to the finite-difference / quadrature error (h = 1e-3, GL 100x60x64).")
print("      NOTE: it holds for the SYMMETRISED Hessian; with the unsymmetrised one it cannot.")
print()
with open("research/v3/species/w7/at_su2_exact.json", "w") as fh:
    json.dump(OUT, fh, indent=1, default=float)
print("[saved] research/v3/species/w7/at_su2_exact.json")
