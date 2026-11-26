#!/usr/bin/env python3
r"""
r4_kp.py -- R4 / "the 2.12%": is u(beta)=c1/c0 the KP criterion, or a coarse surrogate?

Deterministic, numpy + scipy.special.iv, single process, light compute.
Run:  py -3.12 -I r4_kp.py [Pmax]   >  r4_kp.log

Maps 1:1 onto research/v3/species/R4-PERCENT.md.

SECTIONS
  [A] exact SU(2) Wilson data: c_n = I_n - I_{n+2}, sum (n+1)c_n = e^beta,
      u_P3 = c1/c0 = 2I2/I1, u_T5 = I2/I1 = u_P3/2   (E1 factor, exact).
  [B] the four beta values; the four-wall mass m1*a = 4 ln(1/u_tile) in BOTH
      normalisations  (UC-COMPLETE uses u_P3; P3-FOURWALL / Gupta use u_tile = c_f/N).
  [C] EXACT polymer weight of a closed surface, by explicit Haar contraction:
      every internal link integral  int dU U_ab conj(U_cd) = (1/N) delta_ad delta_cb,
      each plaquette trace = 4 chain edges; #cycles C from DSU  =>  w = a1^P N^(C-L).
      Cross-checked against  w = N^chi * u_tile^P  (chi = V-L+P).
  [D] enumeration of closed orientable plaquette surfaces (every link exactly twice,
      opposite orientations) through a fixed plaquette; counts by |X|; entropy.
  [E] the exact KP sums  S(a) = sum_{X contains x} |phi(X)| e^{a|X|}  at beta=2.5,
      truncations |X|<=4 (EMPTY), <=6, <=8, <=10; min over a; comparison with a.
  [F] beta_c^KP: root of min_a S(a)/a = 1, cube multiplicity m in {1,2,8};
      plus the surrogate's own roots (beta_tube, beta_KP) for the contrast.

Labels in the .md: [T] theorem, [C] conditional, [V] numerical witness,
[NG] no-go for a named class, [O] exact missing step.  This witness is [V].
"""

import numpy as np
from scipy.special import iv
from scipy.optimize import brentq, minimize_scalar
import json, sys, time

N = 2  # SU(2)

# ---------------------------------------------------------------- [A] exact data
def Ib(n, b):
    return float(iv(int(n), b))

def ccoef(n, b):
    """exact SU(2) character coefficient, n = #boxes = 2j, dim = n+1"""
    return Ib(n, b) - Ib(n + 2, b)

def c0(b):
    return ccoef(0, b)

def u_p3(b):
    """per-plaquette character-coefficient ratio  c1/c0 = 2 I2/I1"""
    return ccoef(1, b) / c0(b)

def u_tile(b):
    """per-TILE activity for a surface: c_f/(N c0) = u_p3 / N   (Gupta 14.18/14.21)"""
    return u_p3(b) / N

def m1a(ut):
    return 4.0 * np.log(1.0 / ut)

def Phi0(b):
    """surrogate: Phi(0,beta) = sum_{n>=1}(n+1) c_n/c0 = (e^b - c0)/c0  (UC-COMPLETE (star))"""
    return (np.exp(b) - c0(b)) / c0(b)

# ---------------------------------------------------------------- geometry
def _t(x, d, k=1):
    y = list(x); y[d] += k; return tuple(y)

def make_plaq(cyc):
    """canonical plaquette id from a 4-corner cycle (tuples).  Returns (mincorner,mu,nu,s)."""
    m = min(cyc)
    i = cyc.index(m)
    nxt = cyc[(i + 1) % 4]
    prv = cyc[(i - 1) % 4]
    d_out = [k for k in range(4) if nxt[k] != m[k]][0]
    d_in = [k for k in range(4) if prv[k] != m[k]][0]
    mu, nu = min(d_out, d_in), max(d_out, d_in)
    s = +1 if d_out < d_in else -1
    return (m, mu, nu, s)

def plaq_links(pid):
    """cyclic list of (link_key, orient); link_key = (min end, dir).
       Convention: orient=+1 means the edge runs start -> start+e_dir,
       orient=-1 means it runs start -> start-e_dir (start = upper end)."""
    (x, mu, nu, s) = pid
    x = tuple(x)
    if s > 0:
        raw = [(x, mu, +1), (_t(x, mu), nu, +1),
               (_t(_t(x, mu), nu), mu, -1), (_t(x, nu), nu, -1)]
    else:
        raw = [(x, nu, +1), (_t(x, nu), mu, +1),
               (_t(_t(x, nu), mu), nu, -1), (_t(x, mu), mu, -1)]
    return [((site, d) if o > 0 else (_t(site, d, -1), d), o) for (site, d, o) in raw]

def plaq_corners(pid):
    (x, mu, nu, s) = pid
    return [x, _t(x, mu), _t(_t(x, mu), nu), _t(x, nu)]

def plaqs_through_link(key):
    """all 12 signed plaquettes containing the undirected link (site,mu)"""
    (site, mu) = key
    out = []
    for nu in range(4):
        if nu == mu:
            continue
        a = _t(site, nu, -1)
        faceA = [site, _t(site, mu), _t(_t(site, mu), nu), _t(site, nu)]
        faceB = [a, _t(a, mu), _t(_t(a, mu), nu), _t(a, nu)]
        pid = make_plaq(faceA)
        out += [pid, (pid[0], pid[1], pid[2], -pid[3])]
        pid = make_plaq(faceB)
        out += [pid, (pid[0], pid[1], pid[2], -pid[3])]
    return out

def cycle_count(pid_list):
    """exact Haar contraction -> (C, L, V, chi); weight = a1^P * N^(C-L)."""
    P = len(pid_list)
    if P == 0:
        return 0, 0, 0, 0
    nidx = lambda p, s, w: 8 * p + 2 * s + (0 if w == 'a' else 1)
    parent = list(range(8 * P))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    for p in range(P):
        for s in range(4):
            union(nidx(p, s, 'b'), nidx(p, (s + 1) % 4, 'a'))
    seen = {}
    for p, pid in enumerate(pid_list):
        for s, (key, o) in enumerate(plaq_links(pid)):
            seen.setdefault(key, []).append((p, s, o))
    L = 0
    for key, apps in seen.items():
        if len(apps) != 2 or apps[0][2] == apps[1][2]:
            return None
        (p1, s1, _), (p2, s2, _) = apps
        union(nidx(p1, s1, 'a'), nidx(p2, s2, 'b'))
        union(nidx(p1, s1, 'b'), nidx(p2, s2, 'a'))
        L += 1
    C = len({find(i) for i in range(8 * P)})
    V = len({corner for pid in pid_list for corner in plaq_corners(pid)})
    return C, L, V, V - L + P

def box_boundary(x0, dirs, lengths):
    """oriented boundary of the rectangular 3-box spanned by dirs (|dirs|=3), side lengths.
       Local right-handed order (dk, di, dj): low face cycle b->b+di->b+di+dj->b+dj,
       high face reversed.  Canonical ids come from make_plaq, so mu<nu is automatic."""
    x0 = tuple(int(v) for v in x0)
    fs = []
    for k in range(3):
        dk = dirs[k]
        di, dj = dirs[(k + 1) % 3], dirs[(k + 2) % 3]
        Li, Lj = lengths[(k + 1) % 3], lengths[(k + 2) % 3]
        for pos in (0, lengths[k]):
            base = x0
            for _ in range(pos):
                base = _t(base, dk)
            for i in range(Li):
                for j in range(Lj):
                    b = base
                    for _ in range(i):
                        b = _t(b, di)
                    for _ in range(j):
                        b = _t(b, dj)
                    c_lo = [b, _t(b, di), _t(_t(b, di), dj), _t(b, dj)]
                    cyc = c_lo if pos == 0 else [c_lo[0], c_lo[3], c_lo[2], c_lo[1]]
                    fs.append(make_plaq(cyc))
    return fs

# ---------------------------------------------------------------- [D] enumeration
def enumerate_closed(Pmax, time_cap=120.0):
    t0 = time.time()
    p0 = make_plaq([(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0), (0, 1, 0, 0)])
    out = {P: [] for P in range(1, Pmax + 1)}
    seen_sets = {(p0,)}
    timed_out = [False]

    def rec(pids, cnt, tup):
        if time.time() - t0 > time_cap:
            timed_out[0] = True; return
        n = len(tup)
        open_links = [k for k, v in cnt.items() if len(v) == 1]
        if not open_links:
            out[n].append(tup)
        if n >= Pmax:
            return
        if len(open_links) > 4 * (Pmax - n):
            return
        used = {(q[0], q[1], q[2]) for q in pids}   # unsigned plaquette ids already used
        cands = set()
        for k in open_links:
            for q in plaqs_through_link(k):
                if (q[0], q[1], q[2]) not in used:
                    cands.add(q)
        for q in cands:
            ls = plaq_links(q)
            ok = True
            for key, o in ls:
                cur = cnt.get(key, [])
                if len(cur) >= 2 or (len(cur) == 1 and cur[0] == o):
                    ok = False; break
            if not ok:
                continue
            newp = pids | {q}
            ntup = tuple(sorted(newp))
            if ntup in seen_sets:
                continue
            seen_sets.add(ntup)
            ncnt = {k: list(v) for k, v in cnt.items()}
            for key, o in ls:
                ncnt.setdefault(key, []).append(o)
            rec(newp, ncnt, ntup)

    rec({p0}, {k: [o] for k, o in plaq_links(p0)}, (p0,))
    return out, timed_out[0], time.time() - t0

# ---------------------------------------------------------------- main
def main():
    L = []; A = L.append
    A("=" * 78)
    A("r4_kp.py -- R4 : is the 2.12% real, or is the KP criterion used coarsely?")
    A("=" * 78)

    # ---------------- [A]
    A("\n[A] exact SU(2) Wilson data")
    betas = [1 / 48, 0.5, 1.0, 2.5]
    for b in betas:
        cs = np.array([ccoef(n, b) for n in range(60)])
        s = float(np.sum(np.arange(1, 61) * cs))
        A(f"  beta={b:<10.6f} c0={cs[0]:.12f} c1={cs[1]:.12f} "
          f"sum(n+1)c_n/e^b-1={s / np.exp(b) - 1:+.2e}  min_n c_n={cs.min():.3e}")

    # ---------------- [B]
    A("\n[B] the four numbers, both normalisations, and the four-wall mass")
    A(f"  {'beta':>10} {'u_P3=c1/c0':>13} {'u_tile=u_P3/2':>14} "
      f"{'m1a(u_P3)':>12} {'m1a(u_tile)':>13}")
    tab = []
    for b in betas:
        up, ut = u_p3(b), u_tile(b)
        A(f"  {b:>10.6f} {up:>13.10f} {ut:>14.10f} {m1a(up):>12.6f} {m1a(ut):>13.6f}")
        tab.append((b, up, ut, m1a(up), m1a(ut)))
    A("  UC-COMPLETE §2.2 reads m1a = 4 ln(1/u_P3) = -0.057151 at beta=2.5 as 'the tube reversed'.")
    A("  P3-FOURWALL §2.1/§4.2-4.3 fixes the per-tile weight at u_tile = c_f/N (Gupta 14.18/21):")
    A("  m1a = +2.7156.  The whole '2.12%' sits in that factor N=2.")

    b_tube = brentq(lambda bb: u_p3(bb) - 1.0, 1e-8, 60.0)
    b_kp = brentq(lambda bb: u_p3(bb) - 0.5, 1e-8, 60.0)
    b_f0 = brentq(lambda bb: Phi0(bb) - 1.0, 1e-8, 60.0)
    A(f"  surrogate roots: beta_tube = u_P3^-1(1)   = {b_tube:.9f}")
    A(f"                   beta_KP   = u_P3^-1(1/2) = {b_kp:.9f}")
    A(f"                   beta_f0   = Phi0^-1(1)   = {b_f0:.9f}")
    A(f"  u_P3(2.5)={u_p3(2.5):.10f}  u_tile(2.5)={u_tile(2.5):.10f}  "
      f"2.5/beta_tube={2.5 / b_tube:.5f}  ('the 2.12%')")
    A(f"  u_tile(beta)<1 for EVERY beta (->1 as beta->inf): u_tile(50)={u_tile(50):.6f}")
    A("  => 'the four-wall tube does not decay' has NO beta root in the correct normalisation.")

    # ---------------- [C]
    A("\n[C] exact polymer weight of a closed surface (explicit Haar contraction)")
    b = 2.5
    a1 = ccoef(1, b) / c0(b)
    # p0 sanity
    A(f"  p0 = {plaq_links(make_plaq([(0,0,0,0),(1,0,0,0),(1,1,0,0),(0,1,0,0)]))}")
    cube = box_boundary((0, 0, 0, 0), (0, 1, 2), (1, 1, 1))
    r = cycle_count(cube)
    w_cube = a1 ** len(cube) * N ** (r[0] - r[1])
    A(f"  unit cube (3-cell boundary, P=6): C={r[0]} L={r[1]} V={r[2]} chi={r[3]}")
    A(f"    w = a1^6 2^(C-L) = {w_cube:.12e}    N^chi u_tile^6 = {N**r[3]*u_tile(b)**6:.12e}")
    for T in (2, 3):
        tb = box_boundary((0, 0, 0, 0), (3, 0, 1), (T, 1, 1))
        r2 = cycle_count(tb)
        if r2 is None:
            A(f"  capped tube T={T}: not closed -- skipped"); continue
        w2 = a1 ** len(tb) * N ** (r2[0] - r2[1])
        A(f"  capped 4-wall tube T={T} (P={len(tb)}): C={r2[0]} L={r2[1]} V={r2[2]} "
          f"chi={r2[3]}  w={w2:.6e}  N^chi u_tile^P={N**r2[3]*u_tile(b)**len(tb):.6e}")
    A(f"  => at beta=2.5 the minimal closed polymer has w = {w_cube:.4e}, whereas the")
    A(f"     surrogate's |X|=1 term is >= 2 u_P3 e^a = {2*u_p3(2.5):.4f}.")
    A(f"     surrogate over-count = {2*u_p3(2.5)/w_cube:.0f} x ;  |X| shifted 1 -> 6.")

    # ---------------- [D]
    A("\n[D] exhaustive enumeration: closed oriented plaquette surfaces through p0")
    Pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    out, to, dt = enumerate_closed(Pmax, time_cap=90.0)
    A(f"    Pmax={Pmax}  wall={dt:.1f}s  timed_out={to}   (Pmax=7 is COMPLETE: 3.9 s)")
    counts, sumw = {}, {}
    for P in range(1, Pmax + 1):
        lst = out[P]
        if not lst:
            A(f"    |X|={P:2d}: 0"); counts[P] = 0; sumw[P] = 0.0; continue
        ws, chis = [], []
        for tup in lst:
            rr = cycle_count(list(tup))
            if rr is None:
                continue
            ws.append(a1 ** P * N ** (rr[0] - rr[1])); chis.append(rr[3])
        counts[P], sumw[P] = len(lst), float(np.sum(ws))
        A(f"    |X|={P:2d}: {len(lst):7d} surfaces  sum|w|={sumw[P]:.6e}  "
          f"max|w|={max(ws):.4e}  chi in {sorted(set(chis))[:5]}")
    A("    => NO closed polymer with |X|<=5 exists; the minimum is |X|=6 = the 3-cell")
    A("       boundary (= the T=1 four-wall tube WITH its two end caps), 4 of them through p0.")
    A("    constructed box family (not exhaustive, exact weights):")
    fam = {}
    for T in (1, 2, 3):
        tb = box_boundary((0, 0, 0, 0), (3, 0, 1), (T, 1, 1))
        rr = cycle_count(tb)
        w2 = a1 ** len(tb) * N ** (rr[0] - rr[1])
        fam[T] = (len(tb), w2)
        A(f"      box 1x1x{T}: P={len(tb):2d}  w={w2:.6e}  = N^2 u_tile^P "
          f"({N**2*u_tile(b)**len(tb):.6e})   rate per slice = u_tile^{'4'}")
    A(f"    => the capped tube's RATE (per unit height) is u_tile^4, i.e. m1a=4ln(1/u_tile)={m1a(u_tile(b)):.4f}")
    A("       -- u_tile, NOT u_P3.  This is the P3-FOURWALL/Gupta statement.")

    # ---------------- [E]
    A("\n[E] exact KP sums at beta=2.5 :  S(a) = sum_{X contains x} |phi(X)| e^{a|X|}")
    A("    KP (CMP 103 (1986) 491): convergence + exponential clustering iff")
    A("         exists a>0 :  S(a) <= a      [threshold a, NOT 1]")
    A("    convention: x = origin site.  Every surface through p0 contains x, so the")
    A("    per-plaquette sums are rigorous LOWER bounds for the per-site sums:")
    A("    a failure of the per-plaquette sum transfers to the criterion.")
    Ps = [P for P in range(1, Pmax + 1) if counts.get(P, 0) > 0]
    for Pcut in (4, 6, 7):
        Qs = [P for P in Ps if P <= Pcut]
        if not Qs:
            A(f"    |X|<={Pcut}: truncation EMPTY (no closed polymer exists) => S_<={Pcut}(a) = 0")
            continue
        S = lambda a, Qs=Qs: float(sum(sumw[P] * np.exp(a * P) for P in Qs))
        r = minimize_scalar(lambda a: S(a) / a, bounds=(1e-4, 3.0), method='bounded')
        A(f"    |X|<={Pcut} (per-plaquette, {len(Qs)} sizes): min_a S/a = {r.fun:.6e} "
          f"at a*={r.x:.6f}  [need <=1] => {'FAILS' if r.fun > 1 else 'holds'}")
    for m, tag in ((1, 'one cube'), (counts.get(6, 0), 'per-plaquette'), (32, 'per-site')):
        r = minimize_scalar(lambda a: m * w_cube * np.exp(6 * a) / a,
                            bounds=(1e-4, 3.0), method='bounded')
        A(f"    cube-only m={m:2d} [{tag:22s}]: min_a S/a = {r.fun:.6e} at a*={r.x:.6f}"
          f" => {'FAILS' if r.fun > 1 else 'holds'}")
    A("    => at beta=2.5 the MINIMAL polymer alone already violates KP by 1.1x (m=1),")
    A("       4.4x (per-plaquette) , 35.5x (per-site).  Truncation is a LOWER bound => no")
    A("       a>0 works: the KP criterion FAILS at beta=2.5 within the fundamental skeleton.")

    # ---------------- [F]
    A("\n[F] beta_c^KP : the root of  m * N^2 * u_tile(beta)^6 * 6e = 1")
    A("    (the truncated S is a lower bound for S, so this root is an UPPER bound for the")
    A("     true KP root: beta_c^KP <= value below.)")
    roots = {}
    for m, tag in ((1, 'one cube'), (counts.get(6, 0), 'per-plaquette'), (32, 'per-site')):
        ut_c = (1.0 / (m * N ** 2 * 6 * np.e)) ** (1 / 6)
        b_c = brentq(lambda bb: u_tile(bb) - ut_c, 1e-8, 80.0)
        roots[m] = b_c
        A(f"    m={m:2d} [{tag:13s}]: u_tile(beta_c)={ut_c:.6f} -> beta_c^KP <= {b_c:.6f}"
          f"   (2.5/beta_c={2.5/b_c:.4f})   u_P3(beta_c)={u_p3(b_c):.6f}")
    A(f"\n  CONTRAST")
    A(f"    beta_tube(surrogate, = u_P3^-1(1)) = {b_tube:.9f}   ('the 2.12%')")
    A(f"    beta_KP  (surrogate, = u_P3^-1(1/2)) = {b_kp:.9f}")
    A(f"    m=1 root {roots[1]:.6f} vs beta_tube {b_tube:.6f}: {abs(roots[1]-b_tube):.2e} apart")
    A("      -> the 2.12% is a NUMERICAL COINCIDENCE: 4 (N^2) x 6e (min of e^{6a}/a) x u_tile^6")
    A("         lands on u_P3=1.  The two derivations share no step.")
    A(f"    honest bracket:  {b_kp:.4f} < beta_c^KP <= {roots[32]:.4f}  ;  2.5 overshoots by >= "
      f"{2.5/roots[32]-1:.0%}")

    # ---------------- [G]
    A("\n[G] beta-scan of the exact (|X|=6, per-site) KP test  min_a S/a <= 1 ?")
    A(f"    {'beta':>8} {'u_tile':>10} {'w_cube':>12} {'min_a S/a (per-site)':>16} {'verdict':>8}")
    for bb in (0.5, 1.0, 1.044456862, 1.5, 1.778421, 2.0, 2.5, 3.0):
        ut = u_tile(bb); wc = N ** 2 * ut ** 6
        r = minimize_scalar(lambda a: 32 * wc * np.exp(6 * a) / a,
                            bounds=(1e-4, 3.0), method='bounded')
        A(f"    {bb:>8.5f} {ut:>10.6f} {wc:>12.5e} {r.fun:>16.4e} "
          f"{'FAILS' if r.fun > 1 else 'ok':>8}   (a*={r.x:.4f})")
    A("    => the exact criterion's small-size part turns over between beta=1.0445 and 1.5,")
    A("       far below the surrogate's beta_tube=2.4469.  At beta=1.0445 (= the surrogate's")
    A("       'ceiling') the exact sum still has a factor-2 margin: UC-CEIL is NOT a ceiling of the")
    A("       exact KP criterion, only of the surrogate family rho = lambda_a Phi.")

    txt = "\n".join(L)
    print(txt)
    open("r4_kp.log", "w", encoding="utf-8").write(txt + "\n")
    json.dump({"tab": tab, "beta_tube": b_tube, "beta_KP": b_kp, "beta_f0": b_f0,
               "counts": counts, "sumw": sumw, "w_cube": w_cube,
               "fam": {str(k): v for k, v in fam.items()},
               "roots": {str(k): v for k, v in roots.items()}},
              open("r4_kp.json", "w"), indent=1)

if __name__ == "__main__":
    main()
