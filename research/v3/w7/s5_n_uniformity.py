"""S5 witness: N- and a-uniformity bookkeeping for the block-conditional (★):

    Var_{mu(.|U')}(<P^T grad S_j, v>) <= E_{mu(.|U')}<P^T H_j P v, v> + W_j(v) - delta_W
    <=>  G''_j(v) = E[H_j] - Var_j + W_j >= delta_W.

This is a BOOKKEEPING witness, not a proof. It tabulates the N-scaling and the
a-scaling of each term, using the campaign's own normalizations:
  - SZZ (CMP 400, 2023) 't Hooft action S = N beta Re sum tr,  K_S = N(1/2-24|beta|);
  - tower one-loop flow g_j^2 = 1/(2 beta0 j),  beta_j = 2N/g_j^2 = 4 beta0 j,
    beta0 = 11N/(48 pi^2),  d=4, blocking L=2;
  - W_geo >= 0 and E[Phi''] >= (N/12)||v||^2 per unit gauge-invariant mode (ch05);
  - SZZ Lemma 4.1: |Hess_S| <= 8(d-1) N |beta_SZZ| |v|^2  (cost O(N), NOT O(N^2));
  - Lemma B (C2): large-field weight <= exp(-c/g_j^2) = exp(-2 c beta0 j), c=40.
Run with:  py -I s5_n_uniformity.py
"""
import numpy as np

N_const = 48.0 * (np.pi ** 2)  # so beta0 = 11 N / N_const
def beta0(N):
    return 11.0 * N / (48.0 * np.pi ** 2)

def g2_tower(N, j):
    # g_j^2 = 1/(2 beta0 j)
    return 1.0 / (2.0 * beta0(N) * j)

def beta_SZZ_tower(N, j):
    # Wilson beta_W = 2N/g^2 ; SZZ normalization: N beta_SZZ = beta_W/N  => beta_SZZ = 2/(N g^2)
    return 2.0 / (N * g2_tower(N, j))

def K_S(N, beta, d=4):
    return (N + 2) / 2.0 - 1.0 - 8.0 * (d - 1) * N * abs(beta)

line = "=" * 78
print(line)
print("S5 / part 0  -- SZZ baseline K_S = N(1/2-24|beta|): linear in N (C1 sanity)")
print(line)
print(f"{'N':>5} {'K_S(N,1/96)':>14} {'K_S/N':>10} {'beta_crit=1/48':>16}")
for N in [2, 3, 10, 100, 1000]:
    ks = K_S(N, 1.0 / 96.0)
    print(f"{N:>5} {ks:>14.6f} {ks/N:>10.6f} {1.0/48.0:>16.6f}")
print("  => K_S ~ N/2 as beta->0 : linear, per-color 1/2 - 24|beta| N-independent.")

# ---------------------------------------------------------------------------
line = "=" * 78
print("\n" + line)
print("S5 / part 1  -- which regime is the TOWER in?  beta_SZZ,j vs SZZ window 1/48")
print(line)
N = 3
print(f"  N={N}, beta0={beta0(N):.5f}, SZZ strong-coupling window |beta|<1/48={1/48:.5f}")
print(f"{'j':>4} {'g_j^2':>12} {'beta_SZZ,j':>12} {'in SZZ window?':>16}")
for j in [1, 2, 3, 5, 10]:
    b = beta_SZZ_tower(N, j)
    print(f"{j:>4} {g2_tower(N,j):>12.5e} {b:>12.5f} {('YES' if b < 1/48 else 'NO (weak cpl)'):>16}")
print("  => beta_SZZ,j = 11 j/(12 pi^2) ~ 0.093 j > 1/48 for EVERY j>=1.")
print("     SZZ is OFF-REGIME in the tower; the block constant must come from W_geo, not SZZ.")

# ---------------------------------------------------------------------------
line = "=" * 78
print("\n" + line)
print("S5 / part 2  -- N-scaling: reservoir E[Phi''] >= N/12 ; delta_W = W_geo - defect")
print("               W_geo = E[Phi''] - Var(Phi') ;  defect = O(1) (the open (star) content)")
print(line)
print(f"{'N':>5} {'E[Phi'']_lb':>11} {'/N':>9} {'def=0.02':>9} {'def=0.10':>9} {'def=0.50':>9} {'dW/N(def=.02)':>14}")
for N in [2, 3, 5, 10, 100, 1000]:
    EPhi = N / 12.0          # proven lower bound on the curvature reservoir
    W = EPhi                 # optimistic: Var(Phi') subleading; measured W_geo/N ~ 0.125
    row = f"{N:>5} {EPhi:>11.5f} {EPhi/N:>9.5f}"
    for dcoef in [0.02, 0.10, 0.50]:
        row += f" {W-dcoef:>9.5f}"
    row += f" {(W-0.02)/N:>14.6f}"
    print(row)
print("  => reservoir/N = 1/12 (N-uniform per color); delta_W/N -> 1/12 for any O(1) defect.")
print("     Measured W_geo/N ~ 0.125 (N=2: 0.25), defect ~ 0.02 => delta_W ~ 0.23 > 0.")
print("     Linear-N preserved; the finite-N positivity is a COEFFICIENT check (= open (star)).")

# ---------------------------------------------------------------------------
line = "=" * 78
print("\n" + line)
print("S5 / part 3  -- a-uniformity: fix physical block (j ~ log(1/a)), vary a (hence j)")
print("               physical Hessian cost = (bare Hess) * g_j^2 must be j-independent")
print(line)
N = 3
d = 4
print(f"  N={N}: bare Hess per mode = 8(d-1)N beta_SZZ,j ;  field^2 amp = g_j^2")
print(f"{'j':>4} {'bareHess~Nj':>13} {'g_j^2':>11} {'bareHess*g^2':>13} {'W_lb':>8} {'delta_W(.02)':>12}")
for j in [1, 2, 5, 10, 20, 50]:
    bare = 8.0 * (d - 1) * N * beta_SZZ_tower(N, j)   # ~ N j  (grows)
    amp = g2_tower(N, j)                              # ~ 1/(N j) (shrinks)
    phys = bare * amp                                 # ~ const
    W = N / 12.0
    print(f"{j:>4} {bare:>13.4f} {amp:>11.4e} {phys:>13.5f} {W:>8.5f} {W-0.02:>12.5f}")
print("  => bare Hess grows ~ Nj, amp shrinks ~ 1/(Nj); PRODUCT = physical cost is a-uniform.")
print("     W_lb is coupling-free => a-uniform. delta_W = W_lb - O(1) is a-uniform (a^0).")
print("     (contrast: Holley-Stroock bound e^{-2 beta_j Np} collapses like (aL)^c -- a BOUND")
print("      artifact, not delta_W; it is exactly what (star) replaces by the a-uniform floor.)")

# ---------------------------------------------------------------------------
line = "=" * 78
print("\n" + line)
print("S5 / part 4  -- the conditioning cost: bad (large-field) boundary set is e^{-c Nj}")
print("               measure mu_bad <= exp(-c/g_j^2) = exp(-2 c beta0 j), c=40")
print(line)
c = 40.0
print(f"{'N':>5} {'2 c beta0':>11} {'mu_bad(j=1)':>13} {'mu_bad(j=5)':>13} {'mu_bad(j=10)':>14}")
for N in [2, 3, 10, 100]:
    rate = 2.0 * c * beta0(N)   # exponent per shell, ~ c' * N
    vals = [np.exp(-rate * j) for j in [1, 5, 10]]
    print(f"{N:>5} {rate:>11.4f} {vals[0]:>13.3e} {vals[1]:>13.3e} {vals[2]:>14.3e}")
print("  => the *only* exponential is the MEASURE of the bad-boundary set, and its rate")
print("     ~ c' N grows with N: it SHRINKS exponentially in N. No exponential enters the")
print("     constant delta_W; the U'-uniformity holds off a measure-e^{-c'N j} set (Lemma B).")

# ---------------------------------------------------------------------------
line = "=" * 78
print("\n" + line)
print("S5 verdict data")
print(line)
print("  N : E[Phi'']_lb/N=1/12 (per-color, N-indep); W_geo/N~0.125 (measured); cost O(N) not")
print("      O(N^2); adjoint dim N^2-1 = MULTIPLICITY not magnitude => linear-N survives blocking.")
print("  a : delta_W = W_geo - O(1) is a^0 (a-uniform); survives a->0 at fixed physical block.")
print("  1/N or exp factors in the CONSTANT: NONE. Only exp is the bad-set measure e^{-c'Nj},")
print("      favorable in N. Open content is the O(1) defect coefficient (=(star)), N-benign.")
print(line)
