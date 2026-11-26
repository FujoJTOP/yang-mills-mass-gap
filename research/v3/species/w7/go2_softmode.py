# -*- coding: utf-8 -*-
"""GO2-SOFTMODE: compact-resolvent soft mode (Path A) vs relative-bound theorem (Path B).

Pure stdlib, < 2 s, NO Monte Carlo.  Everything is exact diagonal linear algebra
on a truncation of l^2(N); the truncation is chosen so that all reported
Hausdorff distances are EXACT for the infinite sets restricted to [0, D-2]
(the tail beyond N is identical in both spectra after the tail shift).

Construction A (Path A):
    H    = diag(0,1,2,...),  ker H = R e_0,  Delta_cont = 1
    H_a  = diag(0,1,...,N-1, eps, N, N+1, ...)      (transplant e_N -> eps, shift tail down by 1)
    A_a  = H_a - H = -(N-eps)|e_N><e_N| - P_{>N}
Checks: spectrum, Delta_a, Hausdorff d_H, strong (not norm) resolvent convergence,
relative form bound, cost inequality saturation, uniform trace bound.
Also: the "naive" transplant without tail shift (d_H = 1, shows the shift is what saves Hausdorff),
and a bounded-strength rank-one perturbation (no dip: the eps=1 threshold is not reachable
with bounded relative strength).
"""
import bisect
import math
import random

OUT = []


def say(s=""):
    OUT.append(s)


def diag_H(D):
    """H truncated to span(e_0..e_{D-1}):  eigenvalues 0..D-1."""
    return [float(n) for n in range(D)]


def diag_Ha(N, eps, D):
    """H_a truncated to span(e_0..e_{D-1}):  n for n<N; eps at index N; n-1 for n>N."""
    d = [float(n) for n in range(N)]
    d.append(float(eps))
    d.extend(float(n - 1) for n in range(N + 1, D))
    return d


def diag_naive(N, eps, D):
    """Transplant WITHOUT tail shift: n for n<N; eps at N; n for n>N (hole at N)."""
    d = [float(n) for n in range(N)]
    d.append(float(eps))
    d.extend(float(n) for n in range(N + 1, D))
    return d


def delta(d):
    pos = [x for x in d if x > 0.0]
    return min(pos) if pos else float("inf")


def ker_dim(d, tol=1e-12):
    return sum(1 for x in d if abs(x) < tol)


def mult_max(d, tol=1e-9):
    vals = sorted(d)
    best = 1
    run = 1
    for i in range(1, len(vals)):
        if abs(vals[i] - vals[i - 1]) <= tol:
            run += 1
            best = max(best, run)
        else:
            run = 1
    return best


def hausdorff(A, B):
    A = sorted(set(A))
    B = sorted(set(B))

    def dist(x, S):
        i = bisect.bisect_left(S, x)
        c = []
        if i < len(S):
            c.append(abs(S[i] - x))
        if i > 0:
            c.append(abs(S[i - 1] - x))
        return min(c)

    return max(max(dist(x, B) for x in A), max(dist(y, A) for y in B))


def rel_form_bound(N, eps, D):
    """sup_{psi} |<psi,A psi>| / ||H^{1/2}psi||^2  for A = -(N-eps)|e_N><e_N| - P_{>N}."""
    # diagonal entries of A and of H
    A = [0.0] * D
    A[N] = -(N - eps)
    for n in range(N + 1, D):
        A[n] = -1.0
    Hd = diag_H(D)
    best = 0.0
    for n in range(1, D):  # n=0 excluded: ker H, cannot be in form domain with H^{1/2}
        if Hd[n] > 0:
            best = max(best, abs(A[n]) / Hd[n])
    return best


def rel_op_bound(N, eps, D):
    """|| A (H+1)^{-1} ||  (diagonal)."""
    A = [0.0] * D
    A[N] = -(N - eps)
    for n in range(N + 1, D):
        A[n] = -1.0
    best = 0.0
    for n in range(D):
        best = max(best, abs(A[n]) / (n + 1.0))
    return best


def resolvent_apply(h, z, v):
    """((H - z)^{-1} v) for diagonal h."""
    return [v[i] / (h[i] - z) for i in range(len(h))]


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def main():
    D = 40000
    say("=" * 78)
    say("GO2-SOFTMODE verification  (pure stdlib, exact diagonal algebra, no MC)")
    say("=" * 78)

    # ---------- [1] Construction A: spectrum / Delta / Hausdorff ----------
    say()
    say("[1] Construction A: H = diag(0,1,2,...);  H_a = transplant e_N -> eps, tail shifted down")
    say("    N_a = round(1/a), eps_a = a.   Delta_cont = 1.")
    say()
    say("    a          N      eps        Delta_a    d_H(H_a,H)  d_H(naive)  ker_dim  max_mult")
    for a in (1e-1, 1e-2, 1e-3, 1e-4):
        N = int(round(1.0 / a))
        eps = a
        d = diag_Ha(N, eps, D)
        dn = diag_naive(N, eps, D)
        Hh = diag_H(D - 1)  # truncate at D-2 so the tail matches exactly
        dtr = [x for x in d if x <= D - 2.0 + 1e-12]
        dntr = [x for x in dn if x <= D - 2.0 + 1e-12]
        say("    %-10.1e %-6d %-10.1e %-10.1e %-11.1e %-11.1e %-8d %-8d"
            % (a, N, eps, delta(d), hausdorff(Hh, dtr), hausdorff(Hh, dntr),
               ker_dim(d), mult_max(d)))
    say()
    say("    => Delta_a = eps_a -> 0  while Delta_cont = 1:  liminf Delta_a = 0 < 1 = Delta_cont.")
    say("    => d_H = eps_a -> 0 for the shifted family (Hausdorff convergence HOLDS).")
    say("    => the naive transplant has d_H = 1 (hole at N): the tail shift is what saves Hausdorff.")
    say("    => every eigenvalue is simple (max_mult = 1) and ker H_a = R e_0 (ker_dim = 1).")

    # ---------- [2] strong resolvent convergence (not norm) ----------
    say()
    say("[2] Resolvent difference at z = -1:  D_a = (H_a+1)^{-1} - (H+1)^{-1}")
    random.seed(20261126)
    vecs = {}

    def unit(v):
        nrm = norm(v)          # hoisted: DO NOT recompute inside the comprehension (O(D^2))
        return [x / nrm for x in v]

    v = [1.0 if n < 20 else 0.0 for n in range(D)]
    vecs["uniform[0,19]"] = unit(v)
    v = [1.0 / (n + 1.0) for n in range(D)]
    vecs["1/(n+1)"] = unit(v)
    v = [random.gauss(0, 1) for _ in range(D)]
    vecs["random(seed)"] = unit(v)
    say()
    say("    a        N      ||D_a psi||: uniform     1/(n+1)      random     ||D_a||_op")
    for a in (1e-1, 1e-2, 1e-3, 1e-4):
        N = int(round(1.0 / a))
        eps = a
        ha = diag_Ha(N, eps, D)
        h = diag_H(D)
        rowa = []
        for name in ("uniform[0,19]", "1/(n+1)", "random(seed)"):
            psi = vecs[name]
            ra = resolvent_apply(ha, -1.0, psi)
            r = resolvent_apply(h, -1.0, psi)
            rowa.append(norm([ra[i] - r[i] for i in range(D)]))
        opn = max(abs(1.0 / (eps + 1.0) - 1.0 / (N + 1.0)),
                  max(abs(1.0 / n - 1.0 / (n + 1.0)) for n in range(N + 1, D)))
        say("    %-8.1e %-6d %-18.2e %-12.2e %-11.2e %-12.3f"
            % (a, N, rowa[0], rowa[1], rowa[2], opn))
    say()
    say("    => strong resolvent convergence holds (all test vectors -> 0),")
    say("       but ||D_a||_op -> 1: convergence is strong, NOT in norm.")

    # ---------- [3] relative bound: the threshold at 1 ----------
    say()
    say("[3] Relative bounds of A_a = H_a - H  (the sharp threshold is 1)")
    say()
    say("    a        N      eps_form = sup|<A>|/||H^1/2||^2    ||A(H+1)^{-1}||    1-eps_form")
    for a in (1e-1, 1e-2, 1e-3, 1e-4):
        N = int(round(1.0 / a))
        eps = a
        ef = rel_form_bound(N, eps, D)
        eo = rel_op_bound(N, eps, D)
        say("    %-8.1e %-6d %-30.6f %-20.6f %-10.2e" % (a, N, ef, eo, 1.0 - ef))
    say()
    say("    => eps_form = 1 - eps/N = 1 - a^2 -> 1:  the counterexample SATURATES the threshold.")

    # ---------- [4] cost inequality, saturated ----------
    say()
    say("[4] Cost inequality  (1-eps_a)<phi,H phi> <= E_a + b_a  on the soft eigenvector phi=e_N")
    say()
    say("    a        <phi,Hphi>=N   (1-eps_form)*N   E_a=eps   allowed E*/(1-eps_form)")
    for a in (1e-1, 1e-2, 1e-3):
        N = int(round(1.0 / a))
        eps = a
        ef = rel_form_bound(N, eps, D)
        say("    %-8.1e %-15d %-16.6e %-10.1e %-12.1f"
            % (a, N, (1.0 - ef) * N, eps, eps / (1.0 - ef)))
    say()
    say("    => (1-eps_form)*N = eps exactly: the cost inequality is TIGHT.")
    say("    => allowed energy grows like 1/(1-eps): escape needs (1-eps_a) -> 0, i.e. eps_a -> 1.")

    # ---------- [5] uniform trace bound is blind ----------
    say()
    say("[5] tr exp(-t H_a) is uniformly bounded while Delta_a -> 0  (mechanism 3 is blind)")
    say()
    say("    t       tr exp(-tH)      tr exp(-tH_a)  (a=1e-3)   Delta_a")
    for t in (0.5, 1.0, 2.0):
        a = 1e-3
        N = int(round(1.0 / a))
        ha = diag_Ha(N, a, D)
        h = diag_H(D)
        th = sum(math.exp(-t * x) for x in h)
        tha = sum(math.exp(-t * x) for x in ha)
        say("    %-7.1f %-16.6f %-20.6f %-10.1e" % (t, th, tha, delta(ha)))
    say()
    say("    => tr exp(-tH_a) = 1/(1-e^{-t}) + e^{-eps} - e^{-N} -> 1/(1-e^{-t}) + 1  (bounded),")
    say("       while Delta_a = eps_a -> 0.  A uniform trace bound does NOT force an edge atom.")

    # ---------- [6] bounded strength cannot dip (the eps<1 theorem in action) ----------
    say()
    say("[6] Bounded-strength rank-one perturbation cannot reach below Delta_cont")
    say("    H_a = H - c|e_N><e_N| with c = 1/2 FIXED, N -> infinity")
    say()
    say("    N        spectrum gap check:  Delta_a       relative bound c/N")
    for N in (10, 100, 1000):
        c = 0.5
        d = diag_H(D)
        d[N] = N - c  # rank-one attraction at a high level only
        say("    %-8d %-14.6f %-14.3e" % (N, delta(d), c / N))
    say()
    say("    => Delta_a = 1 = Delta_cont for all N (the dip would need c >= Delta_cont = 1);")
    say("       relative bound c/N -> 0.  This is Theorem COST on the safe side of the threshold.")

    # ---------- [7] GN1-A endpoint witness ----------
    say()
    say("[7] GN1-A endpoint witness: limsup sigma(H_a) cap (0, Delta_cont) = empty but liminf Delta_a < Delta_cont")
    say("    limsup is computed as  x in limsup  <=>  dist(x, sigma(H_a)) -> 0 as a -> 0")
    say("    (one-off spectral points eps_a of a finite grid are excluded; a union over a grid is the WRONG object).")
    grid = [1e-1, 1e-2, 1e-3, 1e-4]
    eps_grid = [round(aa, 12) for aa in grid]
    sets = []
    for aa in grid:
        NN = int(round(1.0 / aa))
        sets.append(set(round(x, 12) for x in diag_Ha(NN, aa, 5000)))
    union = set()
    for s in sets:
        union |= s
    finest = sets[-1]
    xs = [round(i * 0.001, 12) for i in range(0, 2001)]
    limsup = []
    for x in xs:
        d_fine = min(abs(x - y) for y in finest)          # distance at a = 1e-4
        one_off = min(abs(x - e) for e in eps_grid)       # is x a one-off eps_a of the grid?
        if d_fine < 1e-9 and one_off > 1e-9:
            limsup.append(x)
    inter = [x for x in limsup if 0.0 < x < 1.0]
    say("    grid a in %s, N = 1/a" % (grid,))
    say("    limsup sigma(H_a) on [0,2]  ~ %s" % ([round(x, 3) for x in limsup[:12]],))
    say("    limsup cap (0,1) = %s   (empty)" % ([round(x, 3) for x in inter],))
    say("    (union-over-grid cap (0,1) = %s -- the WRONG object, printed for contrast)"
        % (sorted(round(x, 4) for x in union if 0.0 < x < 1.0),))
    say("    liminf Delta_a = 0 < 1 = Delta_cont")
    say("    => GN1-A as stated ('(NEW) <=> limsup cap (0,Delta_cont) = empty') is FALSE at the 0-endpoint;")
    say("       GN1-A' with [0,Delta_cont) and sigma cap (0,inf) is the correct form.")

    say()
    say("=" * 78)
    say("ALL CHECKS DONE.  No MC.  exit 0.")
    say("=" * 78)

    text = "\n".join(OUT)
    print(text)
    with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
