"""u06_zo.py -- U06: what IS `Z_O`, and is it computable at L=2, SU(2), beta=2.5?

Definition of record (INVENT-TRANSPORT-INEQ.md L10, ch06-tools.tex L34-36):
    f1^2(mu) = Z_O(mu,mu0)^2 f1^2(mu0),   Z_O(mu,mu0) = (beta_YM/2g)(mu0) / (beta_YM/2g)(mu) > 0
    gamma_O = -d ln(beta_YM/2g)/d ln mu = 2 b0 g^2 + O(g^4),  b0 = 11N/(48 pi^2)
i.e. Z_O is the *operator* (wave-function) renormalization constant of O = (1/4) tr F^2,
fixed by the trace anomaly theta^mu_mu = (beta_YM/2g) O (RG-invariant).

Conventions
-----------
* SU(2) Wilson:  S = beta * sum_p (1 - P_p),  P_p = (1/2) Re tr U_p,  beta = 2N/g^2 = 4/g^2.
* O = (1/4) tr F^2  =>  <O> a^4 = N (1 - <P>)/2 = (1 - <P>)   [N=2]  (tree-level matching).
* beta_eff = 1/48 (P3 strong-coupling base point);  beta_tgt = 2.5.

Three parts
-----------
A. Analytic `Z_O` at the (beta_eff -> 2.5) pair, three evaluation schemes (tree / 1-loop / strong
   coupling power) -- the *magnitude* is scheme-dependent, the *sign* is not.
B. L=2 SU(2) heat-bath + Metropolis cross-check at beta=2.5 (and a beta scan): what a single
   (L,beta) point actually measures, and why that is *not* Z_O.
C. Transport acid test: Z_O^2 * 0 = 0 (an exact homomorphism fixes 0); the floor Z_O^2 r0 moves
   with the scheme -> only the sign is physical.

Run:  py -3.12 w7/u06_zo.py
"""
import time
import numpy as np

RNG = np.random.default_rng(20261125)
N = 2
L = 2


# ======================================================================== A
def part_A():
    print("=" * 78)
    print("A.  ANALYTIC `Z_O` between beta_eff = 1/48 and beta_tgt = 2.5   (SU(2), N=2)")
    print("=" * 78)
    b0 = 11.0 * N / (48.0 * np.pi**2)
    b1 = 34.0 * N**2 / (3.0 * (16.0 * np.pi**2) ** 2)

    def g2(beta):
        return 2.0 * N / beta

    be, bt = 1.0 / 48.0, 2.5
    print(f"  b0 = 11N/(48 pi^2) = {b0:.6f}    b1 = {b1:.6e}")
    print(f"  g^2(1/48) = 2N/beta = {g2(be):.3f}     g^2(2.5) = {g2(bt):.4f}")
    print(f"  1-loop anomalous dimension at the target bare coupling:")
    print(f"    gamma_O = 2 b0 g^2(2.5) = {2*b0*g2(bt):.4f}    "
          f"[at 1/48: {2*b0*g2(be):.2f}]")

    # --- 2-loop asymptotic scaling: the scale ratio between the two beta's
    def lna(beta):                      # ln(a*Lambda), 2-loop, g^2 = 2N/beta
        s = 4.0 * N * b0 / beta         # = 2 b0 g^2
        return -(b1 / (2.0 * b0**2)) * np.log(s) - 1.0 / s

    ratio = float(np.exp(lna(be) - lna(bt)))       # a(1/48)/a(2.5) = mu_tgt/mu_eff
    print(f"  2-loop scaling: a(1/48)/a(2.5) = mu_tgt/mu_eff = {ratio:.1f}")

    # --- scheme 1: tree level (bare lattice coupling at both ends)
    Z_tree = g2(be) / g2(bt)
    # --- scheme 2: 1-loop running of the anomaly coefficient, g^2(mu0) = 192,
    #     1/g^2(mu) = 1/g^2(mu0) + 2 b0 ln(mu/mu0);  beta_YM/2g = -b0 g^2/2 (1-loop)
    #     => Z_O = g^2(mu0)/g^2(mu) exactly at this order.
    inv = 1.0 / g2(be) + 2.0 * b0 * np.log(ratio)
    Z_1loop = g2(be) * inv
    # --- scheme 3 (counter-example, NOT Z_O): R1 sec.3 -- the strong-coupling *residue ratio*
    #     is a power of the scale ratio (INVENT L13: dim-6 residue ~ (a_eff/a_tgt)^6 > 1)
    Z_sc = ratio**6

    print()
    print("  scheme        Z_O        Z_O^2      note")
    print("  " + "-" * 74)
    print(f"  tree         {Z_tree:9.3f}  {Z_tree**2:11.3e}   g^2(1/48)/g^2(2.5) = 192/1.6")
    print(f"  1-loop       {Z_1loop:9.3f}  {Z_1loop**2:11.3e}   g^2(mu0) run to 1/a(2.5), ln-ratio")
    print(f"  power(not Z) {Z_sc:9.3e}  {Z_sc**2:11.3e}   dim-6 residue ratio^6 (R1 s3; NOT Z_O)")
    print(f"  spread: Z_O^2(tree)/Z_O^2(1-loop) = {(Z_tree**2)/(Z_1loop**2):.2f}  "
          f"-- the VALUE is not a datum (scheme-dependent)")
    print(f"  Z_O > 0 in both genuine schemes: {Z_tree > 0 and Z_1loop > 0}   "
          f"[T] ratio of like-signed anomaly coefficients, beta_YM<0, g>0")
    print()
    print(f"  exact 1-loop identity: Z_O - 1 = 2 b0 g^2(mu0) ln(mu/mu0) = "
          f"{2*b0*g2(be)*np.log(ratio):.3f}   => Z_O = 1 + O(g^2)")
    print(f"  lattice tree reading: Z_O = g^2(mu0)/g^2(mu) = beta(mu)/beta(mu0) = "
          f"{bt/be:.0f}  (linear in beta)")
    return dict(b0=b0, b1=b1, Z_tree=Z_tree, Z_1loop=Z_1loop, Z_sc=Z_sc, ratio=ratio)


# ======================================================================== B
def qmul(p, q):
    """SU(2) quaternion product, U = p0 I + i p.sigma; component axis last."""
    p0, p1, p2, p3 = p[..., 0], p[..., 1], p[..., 2], p[..., 3]
    q0, q1, q2, q3 = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([
        p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3,
        p0 * q1 + q0 * p1 - (p2 * q3 - p3 * q2),
        p0 * q2 + q0 * p2 - (p3 * q1 - p1 * q3),
        p0 * q3 + q0 * p3 - (p1 * q2 - p2 * q1),
    ], axis=-1)


def qconj(p):
    return p * np.array([1.0, -1.0, -1.0, -1.0])


def shift(U, mu, sign=+1):
    return np.roll(U, -sign, axis=mu)


class Lat:
    """SU(2) L^4 lattice, U[mu] shape (L,L,L,L,4) quaternion links."""

    def __init__(self, beta, L=2):
        self.L, self.beta = L, beta
        self.U = np.zeros((4, L, L, L, L, 4))
        for mu in range(4):
            self.U[mu] = uniform_su2((L, L, L, L), RNG)
        idx = np.indices((L, L, L, L))
        self.par = (idx.sum(axis=0) % 2)          # checkerboard: sites of equal parity
        # have staples built only from links of the opposite parity -> valid block update

    # ---- action pieces -------------------------------------------------
    def plaquette(self, mu, nu):
        """U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag  (4 links)."""
        a = self.U[mu]
        b = shift(self.U[nu], mu)
        c = qconj(shift(self.U[mu], nu))
        d = qconj(self.U[nu])
        return qmul(qmul(qmul(a, b), c), d)

    def P(self):
        """<P> = mean over 6 planes of (1/2) Re tr U_p = quaternion 0-component."""
        vals = []
        for mu in range(4):
            for nu in range(mu + 1, 4):
                vals.append(self.plaquette(mu, nu)[..., 0])
        return float(np.mean(np.concatenate([v.ravel() for v in vals])))

    def staple(self, mu):
        """w such that the local Boltzmann weight is exp(beta * U.w)  [= exp((beta/2) Re tr(UV))].

        Verified against dS = -beta d(U.w) with dS from the exact plaquette sum (selftest).
        """
        w = np.zeros_like(self.U[mu])
        for nu in range(4):
            if nu == mu:
                continue
            # forward plaquette (base x):  U_mu(x) A,  A = U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
            A = qmul(qmul(shift(self.U[nu], mu), qconj(shift(self.U[mu], nu))),
                     qconj(self.U[nu]))
            w = w + np.stack([A[..., 0], -A[..., 1], -A[..., 2], -A[..., 3]], axis=-1)
            # backward plaquette (base x-nu):  Re tr(U^+ C),  C = U_nu(x-nu)^+ U_mu(x-nu) U_nu(x-nu+mu)
            C = qmul(qmul(qconj(shift(self.U[nu], nu, -1)), shift(self.U[mu], nu, -1)),
                     shift(shift(self.U[nu], nu, -1), mu))
            w = w + C                                     # contribution (C0, +C_vec)
        return w

    # ---- heat bath (checkerboard: one parity class at a time) ----------
    def heatbath_sweep(self):
        for mu in range(4):
            for p in (0, 1):
                w = self.staple(mu)
                kappa = self.beta * np.linalg.norm(w, axis=-1)
                what = w / np.maximum(np.linalg.norm(w, axis=-1), 1e-300)[..., None]
                t = sample_t(kappa, RNG)                   # density ~ sqrt(1-t^2) e^{kappa t}
                r = RNG.normal(size=what.shape)            # uniform random tangent direction
                r = r - np.sum(r * what, axis=-1, keepdims=True) * what
                r = r / np.maximum(np.linalg.norm(r, axis=-1), 1e-300)[..., None]
                u = t[..., None] * what + np.sqrt(np.maximum(1 - t**2, 0))[..., None] * r
                m = self.par == p
                self.U[mu][m] = u[m]

    # ---- independent check of the sampler's target distribution --------
    def dS_exact(self, x, mu, Unew):
        """S(Unew) - S(Uold) from the full plaquette sum (no staple involved)."""
        def S():
            t = 0.0
            for a in range(4):
                for b in range(a + 1, 4):
                    t += np.sum(1.0 - self.plaquette(a, b)[..., 0])
            return self.beta * t
        old = self.U[mu][tuple(x)].copy()
        self.U[mu][tuple(x)] = Unew
        sp = S()
        self.U[mu][tuple(x)] = old
        return sp - S()


def uniform_su2(shape, rng):
    q = rng.normal(size=shape + (4,))
    return q / np.linalg.norm(q, axis=-1, keepdims=True)


def rand_near_identity(shape, eps, rng):
    q = rng.normal(size=shape + (4,)) * eps
    q[..., 0] += 1.0
    return q / np.linalg.norm(q, axis=-1, keepdims=True)


def sample_t(kappa, rng):
    """t in [-1,1], density ~ sqrt(1-t^2) exp(kappa t): exponential proposal + accept."""
    out = np.empty(kappa.shape)
    todo = np.ones(kappa.shape, dtype=bool)
    kap = np.maximum(kappa, 1e-12)
    for _ in range(10000):
        if not todo.any():
            break
        k = kap[todo]
        # proposal: t = 1 + log(1 - u(1-e^{-2k}))/k  (truncated exponential on [-1,1])
        u = rng.random(k.shape)
        t = 1.0 + np.log(1.0 - u * (1.0 - np.exp(-2.0 * k))) / k
        ok = rng.random(k.shape) < np.sqrt(np.maximum(1.0 - t**2, 0.0))
        idx = np.flatnonzero(todo.ravel())
        out.ravel()[idx[ok]] = t[ok]
        todo.ravel()[idx[ok]] = False
    return out


def measure(beta, sweeps=3000, therm=400, mode="hb"):
    lat = Lat(beta)
    P = []
    for s in range(therm + sweeps):
        lat.heatbath_sweep()
        if s >= therm and (s - therm) % 2 == 0:
            P.append(lat.P())
    P = np.array(P)
    n = len(P)
    m = P.mean()
    jk = np.array([np.delete(P, i).mean() for i in range(0, n, max(1, n // 200))])
    err = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2))
    return m, err, P


def part_B(a):
    print()
    print("=" * 78)
    print("B.  L=2, SU(2): what one (L,beta) point measures  -- and why it is NOT Z_O")
    print("=" * 78)

    # validation 1: the staple reproduces the exact local dS from the plaquette sum
    lat = Lat(2.5)
    worst = 0.0
    for _ in range(20):
        x = tuple(int(v) for v in RNG.integers(0, 2, 4))
        mu = int(RNG.integers(0, 4))
        Uold = lat.U[mu][x].copy()
        Unew = qmul(rand_near_identity((), 0.4, RNG), Uold)
        Unew /= np.linalg.norm(Unew)
        w = lat.staple(mu)[x]
        worst = max(worst, abs(lat.dS_exact(x, mu, Unew) + 2.5 * (np.dot(Unew, w) - np.dot(Uold, w))))
    print(f"  [check 1] staple vs exact dS (beta=2.5, L=2): max |err| = {worst:.2e}"
          f"   -> the sampler targets the Wilson action")

    # validation 2: exact strong-coupling anchor  <P> = beta/4 + O(beta^2)
    m, e, _ = measure(0.1, sweeps=800, therm=150)
    print(f"  [check 2] beta=0.1 : <P> = {m:.5f} +- {e:.5f}   exact 1st order beta/4 = "
          f"{0.1/4:.5f}   ({(m-0.1/4)/e:+.2f} sigma)")

    mh, eh, _ = measure(2.5, sweeps=3000, therm=400)
    print(f"  [check 3] beta=2.5 : <P> = {mh:.4f} +- {eh:.4f}   "
          f"(ordered, as it must be at g^2=1.6 on 2^4)")

    # beta scan -> d ln<O>/d beta ; <O>a^4 = (1-<P>) for N=2
    betas = np.array([2.30, 2.40, 2.50, 2.60, 2.70])
    P = np.array([measure(b, sweeps=3000, therm=400)[0] for b in betas])
    O = 1.0 - P
    print()
    print("  beta      <P>        <O>a^4 = 1-<P>     ln(O)")
    for b, p, o in zip(betas, P, O):
        print(f"  {b:.2f}   {p:.5f}     {o:.5f}        {np.log(o):+.5f}")
    i = 2
    dlnO_db = (np.log(O[i + 1]) - np.log(O[i - 1])) / (betas[i + 1] - betas[i - 1])
    print(f"  local slope d ln<O>/d beta at beta=2.5 : {dlnO_db:+.3f}")

    # convert to d ln a via the 2-loop beta function: dlna/dbeta = (b1/2b0^2)/beta - 1/(4N b0)
    b0, b1 = a["b0"], a["b1"]
    dlna_db = (b1 / (2 * b0**2)) / 2.5 - 1.0 / (4 * N * b0)
    db_dlna = 1.0 / dlna_db
    d_eff = dlnO_db * db_dlna
    g2 = 2 * N / 2.5
    print(f"  2-loop: d ln a/d beta = {dlna_db:+.3f}  =>  d beta/d ln a = {db_dlna:+.3f}")
    print(f"  MEASURED scaling dimension  d ln<O>a^4/d ln a = {d_eff:+.2f}")
    print(f"  expected if O had its canonical dim and 1-loop anomalous dim:")
    print(f"    d = 4 + 2 b0 g^2(2.5) = {4 + 2*b0*g2:.3f}   "
          f"(canonical 4.000)")
    print(f"  => the L=2, beta=2.5 point reproduces NEITHER 4 nor 4+2b0g^2 "
          f"(off by a factor {abs(4/d_eff):.1f}) -> outside any scaling window")
    print()
    print("  Single-point content of L=2, beta=2.5:")
    print(f"    <O>a^4 = {1-mh:.4f}   (one number; no second scale, no Z_O)")
    print("    Z_O is by definition a RATIO of the anomaly coefficient at TWO scales")
    print("    (INVENT-TRANSPORT-INEQ L10) -> it is not a function of (L, beta) at all.")
    print("    What a beta-FAMILY measures is gamma_O = -dln(beta_YM/2g)/dln mu, i.e. the")
    print("    integrand of Z_O; on L=2 at beta=2.5 the integrand is not in scaling.")

    # transport acid test on the measured point
    print()
    print("  ACID TEST (transport):  f1^2 -> Z_O^2 f1^2  is a homomorphism, 0 -> 0.")
    for f0 in (0.0, 1e-3):
        print(f"    f1^2(mu0) = {f0:8.4f} -> tree {a['Z_tree']**2*f0:10.4f} | "
              f"1-loop {a['Z_1loop']**2*f0:10.4f} | power {a['Z_sc']**2*f0:12.3e}")
    print("    => the free field (f1^2 = 0) rides to 0 under every Z_O: transport cannot")
    print("       manufacture weight; only the SIGN survives the scheme change.")


if __name__ == "__main__":
    t0 = time.time()
    A = part_A()
    part_B(A)
    print()
    print(f"[done in {time.time()-t0:.1f} s]  verdict:  Z_O defined [T] (INVENT L10 / ch06 L34-36);")
    print("Z_O > 0 [T] (ratio of like-signed anomaly coefficients; Z_O^2 a square);  value =")
    print("scheme-dependent (tree 120 vs 1-loop 84.0, Z_O^2 1.44e4 vs 7.06e3);  NOT measurable at")
    print("(L=2, SU(2), beta=2.5) -- Z_O needs two scales; the missing link of t02 is Omega/(CP),")
    print("not Z_O>0; and Z_O is finite and nonzero as a->0 at fixed physical scales => not (c).")
