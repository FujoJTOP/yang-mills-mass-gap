#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
c05_cbar.py -- independent audit of the claim  Cbar <= Delta = 20  and  Cbar = O(1) in N.
NO reading of a01's derivation.  Own derivation from:
  P3-DRESSED-MASS :63-68   KP (2.1)/(2.2), animal count (eDelta)^n, Delta=20
  P3-FULLREP-KP   :30-31   c_r <= beta^{|r|}/(|r|! dim r) e^beta          (Lemma 1.1)
  P3-FULLREP-KP   :70-72   |I(P)| <= prod_p dim r_p                       (Lemma 2.1)
  P3-FULLREP-KP   :76-78   |w(P)| <= (C_g u)^{|P|}, C_g <= 2N e^beta
  Z06-HOLEA       :15-16   c_f = beta/(2N),  u := c_f/N = beta/(2N^2)
  Z06-HOLEA       :83      |w(S)| <= prod_{p in S} c_f   (PROOF.md Lemma 2')
  Z12-FIXC        :206     beta1^fix = 2N^2/(Chat eDelta Cbar), Chat=4

Checks
  A) Delta = 20 exactly (plaquettes/link = 6, adjacency degree = 4*(6-1) = 20)
  B) Kirchhoff + AM-GM:  t(G) = (1/n) prod_{i>=2} lam_i <= (1/n)(2E/(n-1))^{n-1} <= Delta^n
     -- is 20 or 2*20 the right *tree* constant?  (exhaustive animals + random growth)
  C) exact SU(2) cluster weights -> per-plaquette ratio |w|^{1/n}/u  (lower bound on Cbar)
  D) N-scaling: c_f/u = N ; Phi(0)/u = 2N^2 ; beta1^fix(N) with Cbar=20 vs Cbar=2N^2
  E) arithmetic of a01's numbers (80 vs 7.0633, 11.3x, 1.8394e-3) and Gupta (14.21)
"""
import itertools
import random
from fractions import Fraction

import numpy as np

OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s)


# ---------------------------------------------------------------- A) Delta = 20
def plaquettes_of_torus(L):
    """plaquette = (site, mu<nu)  -> its 4 links = ((site, mu), (site, nu), ...)"""
    plaq = []
    for site in itertools.product(range(L), repeat=4):
        for mu in range(4):
            for nu in range(mu + 1, 4):
                plaq.append((site, mu, nu))
    return plaq


def plaq_links(p, L=8):
    """4 links of plaquette (site,mu,nu) as frozensets of endpoint sites (collision-free)"""
    site, mu, nu = p
    S = tuple(site)
    Smu = list(S)
    Smu[mu] = (Smu[mu] + 1) % L
    Smu = tuple(Smu)
    Snu = list(S)
    Snu[nu] = (Snu[nu] + 1) % L
    Snu = tuple(Snu)
    Smnu = list(Smu)
    Smnu[nu] = (Smnu[nu] + 1) % L
    Smnu = tuple(Smnu)
    return [frozenset({S, Smu}), frozenset({Smu, Smnu}),
            frozenset({Smnu, Snu}), frozenset({Snu, S})]


def check_Delta(L=6):
    plaq = plaquettes_of_torus(L)
    link2plaq = {}
    for p in plaq:
        for lk in set(plaq_links(p, L)):
            link2plaq.setdefault(lk, []).append(p)
    per_link = {len(v) for v in link2plaq.values()}
    adj = {p: set() for p in plaq}
    for lk, ps in link2plaq.items():
        for p in ps:
            for q in ps:
                if p is not q:
                    adj[p].add(q)
    degs = {len(v) for v in adj.values()}
    say(f"[A] 4D torus L={L}: plaquettes={len(plaq)}  plaquettes/link={sorted(per_link)}"
        f"  adjacency degree={sorted(degs)}  => Delta = {max(degs)}")
    return max(degs)


# ------------------------------------------------- B) Kirchhoff + McKay (AM-GM)
def spanning_trees(adj_list):
    """exact integer # spanning trees (Kirchhoff, Bareiss-free: Fraction det)"""
    n = len(adj_list)
    if n == 1:
        return 1
    Lap = [[0] * n for _ in range(n)]
    for i in range(n):
        Lap[i][i] = len(adj_list[i])
        for j in adj_list[i]:
            Lap[i][j] -= 1
    M = [[Fraction(Lap[i][j]) for j in range(n - 1)] for i in range(n - 1)]
    det = Fraction(1)
    for c in range(n - 1):
        piv = None
        for r in range(c, n - 1):
            if M[r][c] != 0:
                piv = r
                break
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        inv = M[c][c]
        for r in range(c + 1, n - 1):
            f = M[r][c] / inv
            if f:
                for k in range(c, n - 1):
                    M[r][k] -= f * M[c][k]
    return int(det)


def animals(seed, adj_of, nmax, seen, out, cur):
    if len(cur) == nmax:
        out.add(frozenset(cur))
        return
    for p in list(cur):
        for q in adj_of(p):
            if q not in cur:
                cur.append(q)
                animals(seed, adj_of, nmax, seen, out, cur)
                cur.pop()


def check_trees(nmax=6, nrand=400, L=8):
    """connected plaquette sets in the L^4 torus; exact spanning trees of adjacency graph"""
    plaq = plaquettes_of_torus(L)
    idx = {p: i for i, p in enumerate(plaq)}
    link2plaq = {}
    for p in plaq:
        for lk in set(plaq_links(p, L)):
            link2plaq.setdefault(lk, []).append(p)
    adj = {p: set() for p in plaq}
    for lk, ps in link2plaq.items():
        for p in ps:
            for q in ps:
                if p is not q:
                    adj[p].add(q)
    Delta = max(len(v) for v in adj.values())
    rows = []
    seed = plaq[0]
    for n in range(1, nmax + 1):
        S = set()
        animals(seed, lambda p: adj[p], n, None, S, [seed])
        # keep only sets containing seed (animals "through a fixed plaquette")
        worst_t = 0
        worst_mckay = 0
        for s in S:
            lst = list(s)
            al = [[idx[q] for q in sorted(adj[p], key=idx.get)] for p in lst]
            # restrict adjacency to vertices inside the set
            pos = {p: i for i, p in enumerate(lst)}
            al = [[pos[q] for q in adj[p] if q in pos] for p in lst]
            t = spanning_trees(al)
            E = sum(len(a) for a in al) // 2
            nn = len(lst)
            mckay = Fraction(1, nn) * Fraction(2 * E, nn - 1) ** (nn - 1) if nn > 1 else Fraction(1)
            worst_t = max(worst_t, t)
            worst_mckay = max(worst_mckay, mckay)
        rows.append((n, len(S), worst_t, float(worst_t) ** (1.0 / n), float(worst_mckay) ** (1.0 / n)))
        say(f"[B] n={n:2d} animals_through_fixed_plaq={len(S):6d} max_t(G)={worst_t:10d}"
            f"  t^(1/n)={float(worst_t) ** (1.0 / n):7.4f}"
            f"  mckay^(1/n)={float(worst_mckay) ** (1.0 / n):7.4f}  Delta={Delta}")
    # random larger clusters
    rng = random.Random(20261125)
    worst = 0.0
    for _ in range(nrand):
        n = rng.randint(7, 12)
        cur = [seed]
        cs = {seed}
        while len(cur) < n:
            p = rng.choice(cur)
            cand = [q for q in adj[p] if q not in cs]
            if not cand:
                break
            q = rng.choice(cand)
            cs.add(q)
            cur.append(q)
        if len(cur) < 3:
            continue
        lst = list(cs)
        pos = {p: i for i, p in enumerate(lst)}
        al = [[pos[q] for q in adj[p] if q in pos] for p in lst]
        t = spanning_trees(al)
        worst = max(worst, float(t) ** (1.0 / len(lst)))
    say(f"[B] random clusters n=7..12 ({nrand} draws): max t^(1/n) = {worst:.4f}  (Delta={Delta}, 2Delta={2*Delta})")
    return Delta, rows


# ---------------------------------------- C) exact SU(2) cluster weights in units of u
def mult_singlet(spins):
    """number of singlets in tensor product of SU(2) irreps (spins in half-integers*2)"""
    cur = {0: 1}  # spin 0 multiplicity map, spins as ints = 2j
    for s in spins:
        nxt = {}
        for t, m in cur.items():
            for u in range(abs(t - s), t + s + 1, 2):
                nxt[u] = nxt.get(u, 0) + m
        cur = nxt
    return cur.get(0, 0)


def closed_surface_multiplicities(T):
    """four-wall tube (1xT rectangular tube surface) in 4D: plaquette set, link incidences.
    returns (n_plaq, prod_link mult(1 in tensor of incident reps))  with reps alternating f/fbar
    so that every shared link carries (f, fbar) => mult = 1."""
    # T temporal walls: 4 walls x T plaquettes; cross-section = 1 plaquette (4 links)
    n_plaq = 4 * T
    # every link of a closed surface is shared by exactly 2 plaquettes, orientation opposite
    # -> mult(1 in f (x) fbar) = 1 at each link
    links = 8 * T - 0  # 4T plaquettes x 4 links / 2
    m = mult_singlet([1, 1])  # f (x) fbar  (spins 1/2,1/2)
    return n_plaq, links, m ** links


def cube_surface():
    """6 plaquettes of a cube (the minimal closed plaquette surface): 12 links, 2 plaq/link"""
    return 6, 12, mult_singlet([1, 1]) ** 12


def check_exact_weights(N=2):
    say("")
    say("[C] exact character-coefficient weights, w(P) = (prod_p c_{r_p}) * prod_links mult(1 in tensor)")
    say("    per-plaquette ratio |w|^{1/n}/u  with u = c_f/N  ==>  this IS the lower bound on Cbar")
    for name, (n, nl, mprod) in [("cube(6 plaq)", cube_surface()),
                                 ("tube T=1 (4 plaq)", closed_surface_multiplicities(1)),
                                 ("tube T=2 (8 plaq)", closed_surface_multiplicities(2)),
                                 ("tube T=3 (12 plaq)", closed_surface_multiplicities(3))]:
        # w = c_f^n * mprod ;  Cbar_plaq := (|w| / u^n)^{1/n} = (c_f/u) * mprod^{1/n}
        ratio = N * mprod ** (1.0 / n)
        say(f"    {name:20s} n={n:3d} links={nl:3d} mult_prod={mprod:6d}"
            f"  c_f/u=N={N}  |w|^(1/n)/u = {ratio:8.4f}")
    say("    => Cbar >= N exactly (mult_prod>=1), i.e. NOT O(1): Cbar >= 2 at N=2, before any tree factor")


# ------------------------------------------------------ D) N-scaling of the KP constants
def per_plaquette_bounds(N, beta=1e-3):
    """c_r dim r <= beta^{|r|}/|r|! e^beta  (Lemma 1.1); c_f = beta/(2N); u = c_f/N"""
    c_f = beta / (2 * N)
    u = c_f / N
    routeA = c_f            # |w| <= prod c_f           (Z06:83, PROOF.md Lemma 2')
    routeB = beta / 2.0     # |w| <= prod (c_f dim f)   (Lemma 2.1, dim f = N)  [leading |r|=1]
    Phi0 = 2 * (beta / 2.0)  # Phi(0) = sum_{r != 0} c_r dim r  ~  beta   (f and fbar)
    return dict(c_f=c_f, u=u, routeA=routeA, routeB=routeB, Phi0=Phi0,
                A_over_u=routeA / u, B_over_u=routeB / u, Phi_over_u=Phi0 / u)


def check_scaling():
    import math
    Delta = 20.0
    eDelta = math.e * Delta
    Chat = 4.0
    say("")
    say("[D] N-scaling of the per-plaquette activity, beta=1e-3 (u = c_f/N = beta/2N^2)")
    say(f"    {'N':>3} {'c_f/u':>7} {'Beta/(2u)':>10} {'Phi0/u':>8} {'Cbar<=B*D':>10} {'Cbar<=A*D':>10}"
        f" {'b1(Cbar=20)':>12} {'b1(Cbar=2N^2*D)':>16}")
    for N in (2, 3, 4, 5, 10):
        d = per_plaquette_bounds(N)
        b1_20 = 2 * N * N / (Chat * eDelta * 20.0)
        CbarB = d["B_over_u"] * Delta
        b1_correct = 2 * N * N / (Chat * eDelta * CbarB)
        say(f"    {N:3d} {d['c_f']/d['u']:7.2f} {d['B_over_u']:10.2f} {d['Phi_over_u']:8.2f}"
            f" {CbarB:10.1f} {d['A_over_u']*Delta:10.1f} {b1_20:12.4e} {b1_correct:16.4e}")
    say(f"    eDelta = e*20 = {eDelta:.4f}")
    say("    KEY: beta1^fix = 2N^2/(Chat eDelta Cbar).  With Cbar = O(1) it scales as N^2 -> WRONG.")
    say("         With Cbar = 2N^2 * Delta it is N-independent, and beta1^fix = 1/(Chat eDelta e^beta).")
    say("    NOTE: Cbar_sharp = Phi(0)/u = 2N^2 (no extra Delta: trees already counted by lambda_a=eDelta)")
    say("          Cbar <= 20 holds only for 2N^2 <= 20  <=>  N <= 3.16  =>  FAILS for N >= 4")
    for N in (2, 3, 4, 10):
        d = per_plaquette_bounds(N)
        csharp = d["Phi_over_u"]
        say(f"          N={N:2d}: Cbar_sharp = 2N^2 = {csharp:6.1f}  <=20? {csharp <= 20}"
            f"   Chat*Cbar_sharp = {Chat*csharp:7.1f}  vs 7.0633 -> shrink {max(1.0, Chat*csharp/7.0633):6.3f}x"
            f"   beta1 = {1.0/(Chat*eDelta):.6e}")


# ------------------------------------------------------------- E) arithmetic of a01
def check_arith():
    import math
    Delta = 20.0
    eDelta = math.e * Delta
    Chat = 4.0
    say("")
    say("[E] arithmetic")
    b1 = 8.0 / (Chat * eDelta * 20.0)
    say(f"    Ĉ*C̄ = {Chat*20:.0f} ; 7.0633 ;  ratio = {Chat*20/7.0633:.3f}"
        f" ;  beta1^fix(20, N=2) = 8/(4 eDelta 20) = {b1:.6e} ; (1/48)/b1 = {(1/48)/b1:.3f}")
    say(f"    check 1.8394e-3 ? {b1:.4e}   (1/48 = {1/48:.6e})")
    # honest window
    b1h = 1.0 / (Chat * eDelta)
    say(f"    honest  beta1 = 1/(Ĉ eDelta) = {b1h:.6e} = (1/48)/{(1/48)/b1h:.3f}  [N-independent]")
    # Gupta (14.21) cross-check
    beta = 1 / 12.0
    N = 2
    u = beta / (2 * N * N)
    say(f"    Gupta (14.21): -4 ln(beta/2N^2), N=2 beta_W=1/12 -> u={u:.6f}, 4ln(1/u)={4*math.log(1/u):.5f}"
        f"   (=DR file's 4ln(1/0.0104)= {4*math.log(1/0.0104):.5f})")
    # exact SU(2) character coefficients
    from math import comb, factorial
    def I(n, b, terms=80):
        # I_n(b) = sum_k (b/2)^{2k+n} / (k! (n+k)!)
        s = 0.0
        for k in range(terms):
            s += (b / 2.0) ** (2 * k + n) / (factorial(k) * factorial(n + k))
        return s
    for b in (0.1, 1 / 12, 0.5):
        c0 = 2 * I(1, b) / b
        u_ex = I(2, b) / I(1, b)
        say(f"    SU(2) beta={b:.5f}: u_exact=I2/I1={u_ex:.6f}  c_f/N={2*I(2,b)/b/2:.6f}"
            f"  c0={c0:.6f}  q_exact={u_ex} vs c_f/N_species={b/(2*2*2):.6f}")
    return b1


def main():
    say("=== c05_cbar.py : independent audit of Cbar <= 20 and Cbar = O(1) ===")
    check_Delta()
    check_trees(nmax=6, nrand=400)
    check_exact_weights(N=2)
    check_scaling()
    check_arith()
    with open("c05_cbar.log", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")
    say("\n[log written: c05_cbar.log]")


if __name__ == "__main__":
    main()
