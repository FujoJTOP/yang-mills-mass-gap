# bn_residual.py -- BN-RESIDUAL witness. N1/N2/N3/N4 of Round 3 (beyond-nu).
# py -3.12 -I bn_residual.py   (numpy/scipy, single process, light)
#
# Blocks:
#   [A] N2: 2D Ising row-transfer gap vs strip width L  (convergence from above)
#   [B] N2: TFIM ring gap vs L at h/J<1 (tunneling -> exponential collapse)
#   [C] N2: reflection positivity => chi(V1 u V2) - chi(V1) - chi(V2) = 2*Cross >= 0
#           + the task's subadditivity chi(2L) <= 2 chi(L) is FALSE for d>=2 (extensive)
#           + free massless: intensive version a-bar(L) unbounded (acid-clean)
#   [D] N1: the five constants; (L,a)-inconsistency of any dimensionless modulus
#   [E] N1: finite-T decay lemma  nu([0,eps)) <= M0 (1+theta) exp(-(c-eps)T) ; T -> inf
#   [F] N3: Fekete accounting -- what strength of deficit is needed (and is impossible)

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh

np.set_printoptions(precision=6, suppress=False, linewidth=140)


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ---------------------------------------------------------------- [A] 2D Ising
def ising_row_gap(L, K):
    """gap of the row-to-row transfer matrix of the 2D Ising model on a strip of
    width L (PBC in the row), free in time, in units of the row spacing."""
    n = 1 << L
    st = np.arange(n)
    s = ((st[:, None] >> np.arange(L)) & 1) * 2.0 - 1.0          # (n,L) spins
    b = np.exp((K / 2.0) * np.sum(s * np.roll(s, -1, axis=1), axis=1))
    A = np.exp(K * (s @ s.T))
    T = b[:, None] * A * b[None, :]
    T = 0.5 * (T + T.T)
    if n <= 1024:
        ev = np.linalg.eigvalsh(T)
        l0, l1 = ev[-1], ev[-2]
    else:
        ev = eigh(T, subset_by_index=[n - 2, n - 1], driver="evr",
                  eigvals_only=True)
        l0, l1 = ev[-1], ev[-2]
    return -np.log(l1 / l0), l0


def ising_exact_mass(K):
    """exact infinite-volume mass of the 2D Ising model, m = 2|K - K*|."""
    Ks = -0.5 * np.log(np.tanh(K))
    return 2.0 * abs(K - Ks), Ks


hdr("[A] N2: 2D Ising row-transfer gap mu1(L)  vs strip width L  (PBC)")
for K in (0.30, 0.60):
    m_ex, Ks = ising_exact_mass(K)
    print(f"\n  K = {K}   (K_c = 0.440687, K* = {Ks:.6f})   exact m_inf = 2|K-K*| = {m_ex:.6f}")
    print("     L     mu1(L)        mu1 - m_inf     monotone?")
    prev = None
    for L in range(2, 11):
        g, l0 = ising_row_gap(L, K)
        mono = "" if prev is None else ("DOWN" if g < prev else "UP  <-- !!")
        print(f"    {L:2d}   {g:.8f}   {g - m_ex:+.3e}    {mono}")
        prev = g
    if K < 0.440687:
        print("    -> disordered phase: gap converges to m_inf FROM ABOVE; "
              "mu1 NON-INCREASING in L.")
    else:
        print("    -> ordered phase: the two lowest levels are the two broken")
        print("       sectors, split by the tunnelling factor; gap -> 0 "
              "EXPONENTIALLY in L.")
    print("    -> in BOTH cases mu1(L) is NON-INCREASING: the task's hypothesis")
    print("       (monotone NON-DECREASING) is refuted by an exact solvable model.")

# ---------------------------------------------------------------- [B] TFIM ring
def tfim_gap(L, J, h):
    n = 1 << L
    st = np.arange(n)
    z = ((st[:, None] >> np.arange(L)) & 1) * 2.0 - 1.0
    diag = -J * np.sum(z * np.roll(z, -1, axis=1), axis=1)
    rows = [np.arange(n)]
    cols = [np.arange(n)]
    vals = [diag]
    for i in range(L):
        rows.append(np.arange(n))
        cols.append(st ^ (1 << i))
        vals.append(np.full(n, -h))
    H = sparse.coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n, n)).tocsr()
    ev = np.sort(eigsh(H, k=3, which="SA", tol=0, return_eigenvectors=False))
    return ev[1] - ev[0]


hdr("[B] N2: TFIM ring gap vs L   H=-J sum sz sz - h sum sx  (tunneling regime)")
for rr in (0.5, 2.0):
    print(f"\n  h/J = {rr}")
    print("     L      gap           log(gap)/L     monotone?")
    prev = None
    for L in range(4, 13):
        g = tfim_gap(L, 1.0, rr)
        mono = "" if prev is None else ("DOWN" if g < prev * 1.0 + 1e-15 else "UP  <-- !!")
        print(f"    {L:2d}   {g:.6e}   {np.log(g)/L:+.6f}    {mono}")
        prev = g
print("  -> h<J: gap collapses exponentially in L (tunneling); in the ordered")
print("     phase there is NO volume-monotonicity in the +L direction.")

# ------------------------------------------------------- [C] RP superadditivity
def green_torus(Lt, d, M):
    k = 2 * np.pi * np.fft.fftfreq(Lt)
    lam = np.full((Lt,) * d, M * M, dtype=float)
    for ax in range(d):
        sh = [1] * d
        sh[ax] = Lt
        lam += 4.0 * np.sin(k.reshape(sh) / 2.0) ** 2
    return 1.0 / lam                                   # hat C(k)


def chi_box(chat, Lt, V):
    """chi(V) = Var_mu( sum_{x in V} phi_x ) = Lt^-d sum_k chat(k) |1_V(k)|^2."""
    d = chat.ndim
    ind = np.zeros((Lt,) * d)
    sl = tuple(slice(0, n) for n in V)
    ind[sl] = 1.0
    Ih = np.fft.fftn(ind)
    return float(np.real(chat * np.abs(Ih) ** 2).sum() / Lt ** d)


def Cbox(chat, Lt):
    """C(z) as a function of the periodic difference (inverse FFT)."""
    return np.real(np.fft.ifftn(chat)) * chat.size / chat.size  # = ifftn * 1


hdr("[C] N2: RP => susceptibility is SUPER-additive under doubling (task's "
    "subadditivity is RP-forbidden)")
for d, Lt, M in ((1, 256, 0.20), (2, 128, 0.20), (4, 32, 0.40)):
    chat = green_torus(Lt, d, M)
    m = Lt // 8
    half = tuple([m] + [Lt] * (d - 1))
    full = tuple([2 * m] + [Lt] * (d - 1))        # a "doubled" box along axis 0
    c1 = chi_box(chat, Lt, half)
    c2 = chi_box(chat, Lt, full)
    cross = 0.5 * (c2 - 2.0 * c1)
    print(f"\n  d={d}  L_torus={Lt}  M={M}   box m={m} -> 2m={2*m}")
    print(f"    chi(V )  (m  ={half}) = {c1:.8f}")
    print(f"    chi(2V)  (2m ={full}) = {c2:.8f}")
    print(f"    cross = (chi(2V)-2chi(V))/2 = {cross:+.8f}   >=0 ? {cross >= -1e-10}")
    print(f"    chi(2V)/chi(V) = {c2 / c1:.6f}")
    print(f"    task's subadditivity  chi(2V) <= 2 chi(V) ?  {c2 <= 2 * c1}")

hdr("[C'] free MASSLESS (zero mode removed): a-bar(R)=sum_{|z|<=R} C(z) is "
    "unbounded  (acid-clean)")
for d in (2, 4):
    Lt = 128 if d == 2 else 32
    chat = green_torus(Lt, d, 0.0)
    chat.flat[0] = 0.0                       # remove the IR zero mode
    C = np.real(np.fft.ifftn(chat))
    C = np.fft.fftshift(C)
    ctr = np.array(C.shape) // 2
    print(f"\n  d={d}: a-bar(R) = sum over |z|_inf <= R of C(z)   (free massless)")
    vals = []
    for R in (2, 4, 8):
        sl = tuple(slice(c - R, c + R + 1) for c in ctr)
        v = C[sl].sum()
        vals.append((R, v))
        print(f"    R={R:2d}   a-bar = {v:.6f}   a-bar/R^2 = {v / R ** 2:.6f}")
    if len(vals) == 3:
        r = (vals[2][1] - vals[1][1]) / (vals[1][1] - vals[0][1])
        print(f"    growth factor a-bar(2R)/a-bar(R) = {r:.4f}  "
              f"(2^d for a purely extensive object = {2 ** d})")
    print("    -> a-bar grows without bound (massless). A functional that is")
    print("       finite/decaying on this object PASSES the free-massless acid test.")
    # control: the same object for a MASSIVE theory must saturate
    chatM = green_torus(Lt, d, 0.40)
    chatM.flat[0] = 0.0
    CM = np.fft.fftshift(np.real(np.fft.ifftn(chatM)))
    print(f"    control (M=0.40, same box sums):")
    for R in (2, 4, 8):
        sl = tuple(slice(c - R, c + R + 1) for c in ctr)
        print(f"      R={R:2d}   a-bar = {CM[sl].sum():.6f}")
    print("      -> the massive object SATURATES; the massless one does not.")
    print("      CAVEAT: zero-mode-subtracted torus propagator, valid for R << L.")

# --------------------------------------------------------- [D] N1: five sources
hdr("[D] N1: the five rate sources -- values, beta/(L,a)-behaviour")
from scipy.integrate import quad
cin = quad(lambda t: (1 - np.cos(t)) / t if t > 0 else 0.0, 0, 2 * np.pi,
           limit=400, points=[0])[0]
print(f"""
  source                         value            beta-dep?   (L,a)-dep?   acid(free)   type
  (1) W_geo>=0 (Prekopa const)   floor 0 (sign)   none        none         PASS=0       sign
  (2) BL floor HessPhi|_0=N/12   N/12  (SU2 1/6)  none        none         same as free point-eval
  (3) trace anomaly 11N/3        b0=11N/3         none        V-extensive  PASS!=0      coefficient
  (4) strong-coupling radius     1/48             threshold   none         regime only  threshold
  (5) exact Haar anchor          2Cin(2pi)={2*cin:.8f}         none        none         same as free anchor
""")
print(f"  check: Cin(2pi) = {cin:.10f},  2*Cin(2pi) = {2*cin:.10f}")
print("  check: N/12 for N=2,3,4 = ",
      [f"{n/12:.6f}" for n in (2, 3, 4)])
print("  (1) == (2):  W_geo's floor IS the Prékopa/Haar Hessian floor N/12 "
      "(AUDIT-R1 A1).  => only FOUR independent sources.")
print("""
  (L,a)-inconsistency argument [T]:
    a modulus c that is a pure dimensionless number yields  m1*a >= c  on every
    lattice. But along the physical line the gap in lattice units VANISHES:
        m1*a = m_phys * a -> 0   as  a -> 0 .
    Hence a fixed c>0 is inconsistent: it would force m1*a >= c > 0 forever,
    i.e. m_phys >= c/a -> infinity. Any (L,a)-consistent floor must have a
    modulus that itself scales like a (or carries a physical scale).""")
for aa in (0.50, 0.25, 0.10, 0.05):
    print(f"    a={aa:5.2f}   m1*a (physical gap 0.7 in fm^-1) = {0.7*aa:.4f}"
          f"   -> a fixed modulus c=0.5 is {'OK ' if 0.7*aa>=0.5 else 'VIOLATED'}")

# --------------------------------------------- [E] N1: finite-T decay lemma
hdr("[E] N1: finite-T decay => floor, WITH the exact rate; T=inf is the gap")
print("""
  Lemma (T).  nu>=0, C(t)=int e^{-s t} dnu, M0=C(0). Suppose
        C(t) <= M0 (1+theta) e^{-c t}      for all t in [0,T].
  Then for every eps in [0,c):
        nu([0,eps))  <=  M0 (1+theta) exp(-(c-eps) T) .
  Proof: C(t) >= nu([0,eps)) e^{-eps t}, plug t=T.
  Corollary: T<inf gives nu([0,eps)) <= eta only for
        T >= log(M0(1+theta)/eta)/(c-eps) ,
  and nu([0,eps))=0 (a floor at eps) requires T=inf. So the "rate" one needs is
  theta(T), the DEFECT of the decay bound -- and theta==0 (all T) IS the gap.
""")
for eta in (1e-2, 1e-4, 1e-6, 1e-10):
    for eps in (0.0, 0.5):
        c = 1.0
        T = np.log(1.0 / eta) / (c - eps)
        print(f"    eta={eta:8.0e}  eps={eps:3.1f}  c={c}  ->  T >= {T:8.3f}")
print("    -> T ~ (1/(c-eps)) log(1/eta):  divergence as eps->0 (the floor is")
print("       the eps=0 endpoint) and as eta->0 (exactness).  No finite T.")
print("\n  witness: two-point family nu_h = (1-h) delta_{s=1} + h delta_{s=0}")
print("     h        m1^2=inf supp   r(t)=-(1/t)log C(t)/M0  at t=1,2,4,8")
for h in (0.0, 1e-6, 1e-2, 0.1):
    for t in (1.0, 2.0, 4.0, 8.0):
        C = (1 - h) * np.exp(-1.0 * t) + h
        r = -np.log(C) / t
        print(f"    {h:8.0e}  {0.0 if h > 0 else 1.0:5.2f}          "
              f"r({t:4.1f})={r:.8f}")
    print("      " + "-" * 60)

# ----------------------------------------------------------- [F] N3/Fekete
hdr("[F] N2/N3: Fekete accounting -- correct version + the exact missing quantity")
print("""
  LEMMA (Fekete accounting, [T]).  a(L)>=0, deficit delta(L):=2a(L)-a(2L).
  Then exactly (divide a(2L)=2a(L)-delta(L) by 2L):
        a(L_n)/L_n = a(L_0)/L_0 - sum_{k<n} delta(L_k)/(2 L_k),  L_k = 2^k L_0.
  (i)   delta >= 0 (subadditive): a(L)/L non-increasing along dyadics, so
        Fekete gives a(L)/L -> inf >= 0  ==>  a(L) <= C L : LINEAR, not bounded.
  (ii)  a <= A bounded  <=>  a(L)/L -> 0  <=>  the slope-drops satisfy
        sum_k [ a(L_k)/L_k - a(L_{k+1})/L_{k+1} ]  =  a(L_0)/L_0  (EXACT equality).
        I.e. the "deficit" bookkeeping is an IDENTITY: any deficit input strong
        enough to conclude boundedness is EQUIVALENT to the conclusion.  There
        is no such thing as "subadditivity + an upper bound => bounded" with the
        two inputs weaker than the output.
  (iii) delta = c (Theta(1)): sum_k c/(2 L_k) < c/L_0 < inf  ==> the slope only
        drops by c/L_0: a(L_n) ~ (a_0/L_0 - c/L_0) L_n, STILL LINEAR/UNBOUNDED.
  (iv)  delta = c L: sum_k c/2 = n c/2 grows without bound ==> a(L_n)/L_n < 0
        for n > 2a_0/(c L_0): the hypothesis SELF-REFUTES.
  ==> There is NO admissible strength of "negative linear term": (iii) is too
      weak and (iv) is self-refuting.  The only remaining option is that the
      deficit equals the extensive part exactly -- i.e. an UPPER (reverse)
      doubling inequality, which is what H1 is.
""")
a0, L0, c = 8.0, 1.0, 1.0
print("     n    a(L_n) [delta=0]  a(L_n) [delta=c]  a(L_n) [delta=cL]"
      "   a/L [delta=c]")
for n in range(0, 7):
    Ln = L0 * 2 ** n
    a_sub = a0 * Ln
    a_c = (a0 - c * (1 - 2.0 ** (-n))) * Ln            # slope a0/L0 - c/L0*(1-2^-n)
    a_lin = a0 * 2 ** n - n * c * L0 * 2 ** (n - 1)
    print(f"    {n:2d}   {a_sub:16.3f}   {a_c:16.3f}   {a_lin:16.3f}"
          f"   {a_c / Ln:10.4f}")
print("""
  THE CORRECT VOLUME STATEMENT (what actually has to be proved, [T]).
  Write a doubling as V + V' across a plane; chi(V u V') = chi(V)+chi(V')+2Cross,
        Cross := sum_{x in V, y in V'} C(x-y)  (the plane-crossing correlation).
  REFLECTION POSITIVITY gives  Cross = <Phi_V , theta Phi_V> = ||.||^2  >= 0
  (block [C] verifies the sign numerically), hence
        a-bar(2V) := chi(2V)/|2V| = a-bar(V) + Cross/|V|   >=  a-bar(V):
  the density is NON-DECREASING under doubling -- the sign RP forces.
  So bounding a-bar needs an UPPER bound on the RP-positive Cross:
        Cross(V,V') <= C |boundary| = O(L^{d-1})      ==>  Cross/|V| = O(1/L) -> 0.
  But Cross = O(L^{d-1}) IS "exponential clustering across the plane", i.e. the
  gap.  Precisely: the missing inequality is
        Cross(V,V') / |V|  ->  0 ,
  and by positivity this is equivalent to m1 > 0.  The "linear term" the task
  asks about is the extensive part 2*chi(V); the "deficit" is Cross; RP fixes
  its SIGN (positive), so no term of the right sign exists.
""")

# --------------------------------------------------- [G] N3: the sums of eps
hdr("[G] N3: two-scale transport -- what eps( beta ) the polymer route can give")
u25 = 1.01439
print(f"""
  Inputs (T1-POLYMER-RG, logged): at beta=2.5
     activity u = {u25:.5f}   (need u < 1/lambda_a = 0.01839 for KP; {u25/0.01839:.1f}x over)
     Ursell/KP ratio rho_Ur = 4.6e3 (need < 1)
     xi/a = 1.3  ==> block-spin steps available = log2(xi/a) in [-0.9, 0.4] = 0..1
  Consequence 1 (stronger than the brief): the Ursell bound on the step error is
     |eps| <= sum_{{n>=1}} u^n = u/(1-u) = {u25/(1-u25):.4f}
  which is NEGATIVE / meaningless: the geometric majorant has no value. So the
  polymer route does not give eps = O(1); it gives eps UNDEFINED.
  Consequence 2: even with a hypothetical eps = O(1) in lattice units,
     m1*a at beta=2.5 lies in [0.53,1.88] (T1 witness G), so eps/(m1*a) = O(1):
     no gain from one step, and there are 0-1 steps to run.
  Consequence 3: eps must be multiplied by the number of steps and compared with
  m1(beta_0) at beta_0<1/48; since the ladder beta_0 -> 2.5 has O(1) steps and
  each eps = O(1)*m1, the transport inequality m1(beta') >= m1(beta)-eps is
  VACUOUS rather than merely non-summable.
""")
print("  Ursell-majorant values u/(1-u) at the logged activities:")
for u in (0.0104, 0.480, 0.988, 1.014):
    val = u / (1 - u)
    print(f"    u={u:7.4f}  ->  eps_bound = {val:+12.4f}"
          f"   {'DIVERGES/MEANINGLESS' if u >= 1 else 'ok'}")

hdr("DONE")
