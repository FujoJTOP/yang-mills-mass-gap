"""Single-plaquette Hessian bound for the Wilson action on SU(N).

Purpose: test the structural premise that on a COMPACT group the Hessian of the
Wilson action is bounded BELOW by a field-INDEPENDENT constant times beta
(contrast: scalar phi^4 has V'' = 3*lam*phi^2 + mu, unbounded below as |phi|->inf).

Setup.  One plaquette, one distinguished link U1 = exp(X), X = sum_a x_a * lam_a,
lam_a = i*T_a with tr(T_a T_b) = delta_ab/2.  Then g(X,Y) := -2 tr(X Y) makes
{lam_a} orthonormal, i.e. the x-coordinates are an orthonormal frame for g.
With M := U2 U3 U4 U1 (fixed product of the other links, unitary),
    S_p = beta * (1 - (1/N) Re tr(U1 U2 U3 U4)) = beta*(1 - (1/N) Re tr(e^X M)),
and the Hessian at X = 0 is EXACT (no finite differences):
    H_ab = -(beta/N) * Re tr(lam_a lam_b M)      [second variation of -Re tr(e^X M)]
because tr(e^X M) = tr M + tr(XM) + (1/2)tr(X^2 M) + O(X^3) and every higher
term contributes 0 to the second derivative at X = 0.

We report, over Haar-random M in SU(N):
    Cn(N) := sup over M of ( - lambda_min(H) * N / beta )   (per-plaquette constant)
so that Hess S_p >= -(Cn/N) * beta * g  for every configuration, and a link in
d=4 belongs to 2(d-1)=6 plaquettes, hence
    Hess S_W >= -(6*Cn/N) * beta * g   uniformly in the configuration.
"""
import numpy as np

rng = np.random.default_rng(20261125)


def su_generators(N):
    """Generalized Gell-Mann T_a, hermitian, tr(T_a T_b) = delta_ab / 2."""
    Ts = []
    # off-diagonal symmetric
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex)
            T[j, k] = T[k, j] = 0.5
            Ts.append(T)
    # off-diagonal antisymmetric
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex)
            T[j, k] = -0.5j
            T[k, j] = 0.5j
            Ts.append(T)
    # diagonal (Cartan)
    for l in range(1, N):
        T = np.zeros((N, N), dtype=complex)
        d = 1.0 / np.sqrt(2.0 * l * (l + 1))
        for j in range(l):
            T[j, j] = d
        T[l, l] = -l * d
        Ts.append(T)
    return Ts


def haar_su(N, rng):
    """Haar-random SU(N) via QR of a complex Ginibre matrix."""
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R)
    ph = d / np.abs(d)
    Q = Q * ph.conj()[None, :]
    det = np.linalg.det(Q)
    return Q / det ** (1.0 / N)


def hess_constants(N, nsamp=4000):
    Ts = su_generators(N)
    lam = [1j * T for T in Ts]
    n = len(lam)
    # lam_a lam_b for all pairs, flattened
    prod = np.empty((n, n, N, N), dtype=complex)
    for a in range(n):
        for b in range(n):
            prod[a, b] = lam[a] @ lam[b]
    neg = 0.0
    pos = 0.0
    for _ in range(nsamp):
        M = haar_su(N, rng)
        H = np.empty((n, n))
        for a in range(n):
            for b in range(a, n):
                v = np.real(np.trace(prod[a, b] @ M))
                H[a, b] = v
                H[b, a] = v
        ev = np.linalg.eigvalsh(H)
        neg = max(neg, -ev[0])
        pos = max(pos, ev[-1])
    return n, neg, pos


print("Hess S_W  >=  -(6 * sup(-lam_min) / N) * beta * g   per link, d=4 (6 plaquettes per link)")
print("N  dim(G)  sup(-lam_min)   sup(lam_max)   coef4D=6*sup(-lmin)/N   coef4D*2.5    (1/2)Ric=N/8   ratio(act/geom)")
rows = []
for N in (2, 3, 4, 5):
    dim, neg, pos = hess_constants(N, nsamp=1500 if N >= 4 else 4000)
    coef = 6.0 * neg / N
    halfric = N / 8.0
    rows.append((N, dim, neg, pos, coef, coef * 2.5, halfric))
    print(f"{N}  {dim:5d}  {neg:12.6f}  {pos:12.6f}  {coef:20.6f}  {coef*2.5:12.6f}  {halfric:12.6f}  {coef*2.5/halfric:8.2f}")

print()
print("beta=2.5: lambda_dot_0 >= -(coef4D*beta) + (1/2)Ric  (the t=0 value the BD criterion must start from)")
for N, dim, neg, pos, coef, lb, halfric in rows:
    print(f"  SU({N}): lambda_dot_0 >= {(-lb + halfric):9.4f}   [action {(-lb):8.4f}, geometry {halfric:2.4f}]")
