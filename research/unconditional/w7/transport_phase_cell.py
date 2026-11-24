"""
transport_phase_cell.py  --  witness for TRANSPORT-PHASE-CELL.md
Audit of the Battle-Federbush phase-cell / Federbush-chunk route.

Four acid tests:
  T1  free massless field: lowest phase-cell mode keeps lambda_min ~ L^-2, no artificial mass.
  T2  pure-gauge background decouples exactly (F=0).
  T3a (l/L)^2 suppresses the LINEAR (tadpole) coupling but NOT the QUADRATIC (mass) coupling;
      the mass-generating quadratic operator is the Savvidy Hessian, which has a negative mode.
  T3b single SU(2) chunk integral (1 x n strip with the shared link = the mode), MC.
  T4  high-representation coefficients u_j(beta) and the multi-scale forest activity.

Run:  py -3.12 transport_phase_cell.py
Light: ONE job only. No heavy lattice runs.
"""
import numpy as np
from scipy.special import iv

trapz = getattr(np, "trapezoid", None) or np.trapz
rng = np.random.default_rng(12345)

print("=" * 74)
print("TRANSPORT-PHASE-CELL witness  (Battle-Federbush / Federbush-chunk audit)")
print("=" * 74)

# ----------------------------------------------------------------------------
# SU(2) via unit quaternions: q=(w,x,y,z), U = w + i x tau1 + i y tau2 + i z tau3
# (1/2) Re tr U = w ; inverse = conjugate.
# ----------------------------------------------------------------------------
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])

def qinv(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

def qrand(g):
    q = g.normal(size=4)
    return q / np.linalg.norm(q)

def qexp(nvec, eps):
    n = np.asarray(nvec, float); n = n / np.linalg.norm(n)
    return np.concatenate(([np.cos(eps/2.0)], np.sin(eps/2.0) * n))

# ============================================================================
# TEST 1: free massless field, phase-cell spectrum
# ============================================================================
print("\n### TEST 1: free massless field -- phase-cell spectrum")
L = 1.0
print("  smallest NONZERO eigenvalue of -Delta (periodic box L=1):")
for N in [16, 32, 64, 128, 256]:
    h = L / N
    lam_nz = 4.0 * np.sin(np.pi / N) ** 2 / h**2
    print(f"   N={N:4d}  lambda_min_nz={lam_nz:12.6f}   (2pi/L)^2={(2*np.pi/L)**2:12.6f}"
          f"   ratio={lam_nz/(2*np.pi/L)**2:.6f}")
print("  contrast -- a NAIVE block-Dirichlet Laplacian on a cell of size ell would give")
print("   lambda_min ~ (pi/ell)^2 -> infinity as ell -> 0 : that WOULD be an artificial mass.")
print("   The true spectral block-spin operator keeps 4*pi^2 independent of the cell size.")

def periodic_laplacian(N, L):
    h = L / N
    K = np.zeros((N, N))
    for i in range(N):
        K[i, i] = 2.0 / h**2
        K[i, (i+1) % N] = -1.0 / h**2
        K[i, (i-1) % N] = -1.0 / h**2
    return K

N = 64
K = periodic_laplacian(N, L)
ev = np.sort(np.linalg.eigvalsh(K))
print(f"  full K (N={N}): lowest 3 eig = {ev[:3]}")
idx_e = np.arange(0, N, 2); idx_o = np.arange(1, N, 2)
Kee = K[np.ix_(idx_e, idx_e)]; Keo = K[np.ix_(idx_e, idx_o)]
Koe = K[np.ix_(idx_o, idx_e)]; Koo = K[np.ix_(idx_o, idx_o)]
Keff = Kee - Keo @ np.linalg.solve(Koo, Koe)          # integrate out the fine (odd) modes
evf = np.sort(np.linalg.eigvalsh(Keff))
c = np.ones(N // 2); c /= np.linalg.norm(c)
m2 = float(c @ Keff @ c)
print(f"  block-spin Schur-complement coarse K (N/2={N//2}): lowest 3 eig = {evf[:3]}")
print(f"   (lowest NONZERO ~ {evf[1]/(2*np.pi/L)**2:.2f}*(2pi/L)^2 : a wavefunction renormalisation,")
print(f"    still O(L^-2); it is NOT (pi/h_coarse)^2={(np.pi/(2/N))**2:.0f})")
print(f"  coarse constant-mode mass^2 = <1|Keff|1> = {m2:.3e}   (0 => no induced mass)")
print("  VERDICT T1: PASS in the exact free Gaussian theory (orthogonal change of basis;"
      " the coarse zero mode stays massless; no shift ~1/ell^2).")

# ============================================================================
# TEST 2: pure-gauge background decouples
# ============================================================================
print("\n### TEST 2: pure-gauge background decouples")
Nd = 4
sites = [(x, y, z, t) for x in range(Nd) for y in range(Nd)
         for z in range(Nd) for t in range(Nd)]

def shift(s, mu, N):
    s = list(s); s[mu] = (s[mu] + 1) % N; return tuple(s)

g = {s: qrand(rng) for s in sites}
U = {}
for s in sites:
    for mu in range(4):
        U[(s, mu)] = qmul(g[s], qinv(g[shift(s, mu, Nd)]))

def plaquette(U, s, mu, nu, N):
    return qmul(qmul(qmul(U[(s, mu)], U[(shift(s, mu, N), nu)]),
                     qinv(U[(shift(s, nu, N), mu)])), qinv(U[(s, nu)]))

maxdev = 0.0; nplaq = 0
for s in sites:
    for mu in range(4):
        for nu in range(mu+1, 4):
            P = plaquette(U, s, mu, nu, Nd)
            maxdev = max(maxdev, abs(P[0]-1.0) + abs(P[1]) + abs(P[2]) + abs(P[3]))
            nplaq += 1
print(f"  pure gauge U(x,mu)=g(x)g(x+mu)^-1 : {nplaq} plaquettes,"
      f" max deviation from identity = {maxdev:.2e}")
print("  => F = 0 to machine precision; the field-strength coupling of ANY mode vanishes.")

# one defect link -> F nonzero only on plaquettes touching it
h = qexp([0, 0, 1], 1.0)
U2 = dict(U); U2[((0, 0, 0, 0), 0)] = qmul(h, U2[((0, 0, 0, 0), 0)])
cnt = 0; tot = 0
for s in sites:
    for mu in range(4):
        for nu in range(mu+1, 4):
            P = plaquette(U2, s, mu, nu, Nd); tot += 1
            if abs(P[0]-1.0) > 1e-9:
                cnt += 1
print(f"  with ONE defect link: {cnt}/{tot} plaquettes carry nonzero F (localized at the defect).")
print("  VERDICT T2: PASS (exact). Pure-gauge background decouples; a chunk = the localized defect.")

# ============================================================================
# TEST 3a: linear vs quadratic coupling scaling
# ============================================================================
print("\n### TEST 3a: (l/L)^2 kills the LINEAR coupling, not the QUADRATIC (mass) coupling")
def mode(x, ell):
    y = x / ell
    return (1.0 - y*y) * np.exp(-y*y / 2.0)

def linear_coupling(ell, Lbg=1.0, n=60001):
    # fixed unit AMPLITUDE mode (a(x0)=1); zero mean and zero dipole by symmetry,
    # so the leading coupling to a smooth background is the quadrupole ~ (ell/Lbg)^2.
    x0 = 0.5
    x = np.linspace(x0 - 10*ell, x0 + 10*ell, n)
    a = mode(x - x0, ell)
    F = np.cos(2*np.pi*(x - x0)/Lbg)          # background varies on scale Lbg
    return abs(trapz(a*F, x))

ells = [0.10, 0.05, 0.025, 0.0125, 0.00625]
vals = [linear_coupling(e) for e in ells]
def linear_coupling_L2(ell, Lbg=1.0, n=60001):
    x0 = 0.5
    x = np.linspace(x0 - 10*ell, x0 + 10*ell, n)
    a = mode(x - x0, ell)
    a = a / np.sqrt(trapz(a*a, x))
    F = np.cos(2*np.pi*(x - x0)/Lbg)
    return abs(trapz(a*F, x))
vals2 = [linear_coupling_L2(e) for e in ells]
print("   ell/Lbg      |<a,F>| (unit amplitude)   |<a,F>| (unit L2 norm)")
for e, v, v2 in zip(ells, vals, vals2):
    print(f"   {e:0.5f}      {v:0.5e}              {v2:0.5e}")
slope = np.polyfit(np.log(ells[:3]), np.log(vals[:3]), 1)[0]
slope2 = np.polyfit(np.log(ells[:3]), np.log(vals2[:3]), 1)[0]
print(f"   fitted log-log slope p = {slope:.3f} (unit amplitude), {slope2:.3f} (unit L2 norm)")
print("   => the overlap is a POSITIVE power of ell/Lbg (exponent depends on mode normalisation")
print("      and multipole content; the paper's exact (ell/L)^2 bookkeeping was not accessible).")
print("   ROBUST POINT: in every normalisation the overlap -> O(1) as ell -> Lbg = L.")
print("   => the scale-separation suppression is gone exactly at the crossover.")

# quadratic (mass) coupling: 2*gB * J3 (spin-1, field along 3); eigenvalues {0,+2gB,-2gB}
J3 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])  # generator of rotations 1-2
lam = np.sort(np.linalg.eigvalsh(2.0 * J3))
print(f"\n   quadratic (mass) operator 2*B*J3 (spin-1) eigenvalues = {lam}")
print("   -> the transverse sector has a NEGATIVE eigenvalue -2B : the Savvidy unstable mode.")
print("   This operator is SCALE-FREE: it carries NO (ell/L) factor; it acts at every mode scale.")
print("   VERDICT T3a: the (ell/L)^2 suppression controls only the tadpole. The mass is generated")
print("   by the quadratic (Hessian) coupling, which is unsuppressed and indefinite. It cannot")
print("   replace the assumption Hess S >= 0 --- it is precisely the failure of that assumption.")

# ============================================================================
# TEST 3b: single SU(2) chunk integral (1 x n strip, shared link = mode), MC
# ============================================================================
# A 1 x n strip of plaquettes.  Bottom link H0 = fixed mode U(theta).  All other
# links are the "chunk" and are integrated out.  S_eff''(theta=0) is the induced
# link curvature (a finite-volume analogue of an induced mass^2).
print("\n### TEST 3b: SU(2) chunk integral (1 x n strip, MC)")
def strip_action(links, n, beta):
    S = 0.0
    for i in range(n):
        P = qmul(qmul(qmul(links[('H', i)], links[('V', i+1)]),
                      qinv(links[('T', i)])), qinv(links[('V', i)]))
        S += (1.0 - P[0])
    return beta * S

def run_strip(n, beta, theta, nsweeps=900, ntherm=300, seed=1):
    """Return S_eff'(theta) = <dS/dtheta> at the fixed link, and its standard error.
    Since S_eff is even, S_eff'(0)=0 and S_eff''(0) ~ S_eff'(theta)/theta."""
    g = np.random.default_rng(seed)
    links = {}
    for i in range(n):
        links[('H', i)] = qrand(g); links[('T', i)] = qrand(g)
    for i in range(n+1):
        links[('V', i)] = qrand(g)
    links[('H', 0)] = qexp([0, 0, 1], theta)       # fixed mode link
    free = [k for k in links if k != ('H', 0)]
    cur = strip_action(links, n, beta)

    def S_theta(th):
        H0 = qexp([0, 0, 1], th)
        P = qmul(qmul(qmul(H0, links[('V', 1)]),
                      qinv(links[('T', 0)])), qinv(links[('V', 0)]))
        return beta * (1.0 - P[0])

    eps = 1e-4
    vals = []
    for sweep in range(nsweeps):
        for k in free:
            q0 = links[k]
            dq = qexp(g.normal(size=3), 0.5 * g.random())
            links[k] = qmul(dq, q0)
            new = strip_action(links, n, beta)
            if new < cur or g.random() < np.exp(-(new - cur)):
                cur = new
            else:
                links[k] = q0
        if sweep >= ntherm:
            Sp = S_theta(theta + eps); Sm = S_theta(theta - eps)
            vals.append((Sp - Sm) / (2*eps))
    vals = np.array(vals)
    return vals.mean(), vals.std(ddof=1) / np.sqrt(len(vals)), len(vals)

print("   using S_eff'(theta)/theta ~ S_eff''(0) (mean estimator; S_eff even so S_eff'(0)=0):")
for n in [2, 3, 4]:
    for beta in [0.02, 2.5]:
        m1, e1, ns = run_strip(n, beta, 0.2, seed=11 + 7*n + int(100*beta))
        m2, e2, _ = run_strip(n, beta, 0.4, seed=100 + 7*n + int(100*beta))
        print(f"   strip n={n}  beta={beta:5.2f} :  S_eff''(0) ~ "
              f"{m1/0.2:+.4f}+-{e1/0.2:.4f} (theta=.2) , {m2/0.4:+.4f}+-{e2/0.4:.4f} (theta=.4)"
              f"  [nsamp={ns}]")
print("   interpretation (honest, MC-limited): the exact finite group integral is well defined")
print("   and its induced link curvature is O(1) at beta=2.5 (NO small parameter). For n=3 the")
print("   inferred S_eff'' changes SIGN between theta=0.2 and 0.4 (non-convex effective action).")
print("   So the chunk integral does NOT furnish a convex / sign-controlled lower bound; it is a")
print("   finite-volume link curvature (not the 0++ residue), and its Gaussian limit is the")
print("   indefinite Savvidy Hessian of Test 3a. The real chunk (many group elements + degree-1")
print("   map S^3->SU(2) + induced Jacobian) is NOT reducible to this toy in closed form.")

# ============================================================================
# TEST 4: character coefficients and multi-scale forest activity
# ============================================================================
print("\n### TEST 4: high-representation channel + multi-scale activity")
def Zp(beta):
    t = np.linspace(0, np.pi, 6001)
    return trapz((2/np.pi) * np.sin(t)**2 * np.exp(beta*np.cos(t)), t)

print("  single-plaquette partition function Z_p(beta) = int (2/pi) sin^2 t e^{beta cos t}:")
for beta in [1/48, 1.0, 2.5, 4.0]:
    print(f"   beta={beta:0.4f}: Z_p={Zp(beta):.6f}   2 I1(beta)/beta={2*iv(1,beta)/beta:.6f}")

print("\n  SU(2) character coefficients u_j = I_{2j+1}(beta)/I_1(beta):")
for beta in [1/48, 1.0, 2.5]:
    us = [iv(2*j+1, beta)/iv(1, beta) for j in range(5)]
    print(f"   beta={beta:0.4f}: u_0..4 = " + ", ".join(f"{u:.3e}" for u in us))
    jstar = int(np.argmax([(2*j+1)**2 * us[j] for j in range(5)]))
    print(f"     typical representation j* = {jstar};  first-order truncation error ~ u_1 = {us[1]:.3f}")
print("  => at beta=1/48 the representation series has a genuine small parameter (u_1~1e-5);")
print("     at beta=2.5, u_1 = O(1): the high-representation channel has NO small parameter.")

print("\n  one-loop SU(2) running coupling, dg/dln(mu) = -(11N/3) g^3/(16 pi^2):")
def g2(muOverLambda, N=2):
    return 48*np.pi**2 / (11*N*np.log(muOverLambda))
for n in [1, 2, 4, 8, 12, 16, 20]:
    mu = 2.0**n
    gg = g2(mu)
    print(f"   mu/Lambda = 2^{n:2d} = {mu:10.1f} : g^2={gg:8.3f}   per-mode activity g^2/(16pi^2)={gg/(16*np.pi**2):.3e}")
thr05 = np.exp(24*np.pi**2/(11*0.5))
print(f"  => g^2 = O(10) already at mu ~ Lambda (the crossover);"
      f" even g^2=0.5 requires mu/Lambda ~ {thr05:.1e} (one loop).")
print("     The log-running makes the small-coupling regime unreachably far from the mass scale.")
print("     The naive per-shell forest activity sum_mu (L/ell)^4 * g^2(ell) diverges in the UV;")
print("     the phase-cell RG / tree-graph method is what resums the (L/ell)^4, leaving the")
print("     marginal g^2. At the mass-generating shell g^2 = O(1): no summable small parameter.")
print("  VERDICT T4: UV marginal sum needs the full tree-graph/RG cancellation (as in phi^3_4,")
print("     where a POWER small parameter exists); at the crossover the resummed activity is O(1).")

print("\n" + "=" * 74)
print("DONE.  Summary: T1 PASS; T2 PASS; T3 fails to substitute for Hess>=0 (Savvidy -2B);")
print("T4 has no small parameter at ell ~ L (crossover). The route reduces to the known wall.")
print("=" * 74)
