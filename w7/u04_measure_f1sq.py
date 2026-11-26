#!/usr/bin/env python3
r"""
u04_measure_f1sq.py -- measure f1^2 = lim_{t->inf} e^{E1 t} C(t) for
    C(t) = <O(t) O(0)>_c ,   O = (1/4) tr F_{mu nu} F_{mu nu}   (0++ / A1++ channel)
on SU(2) Wilson lattice, L=2 (and an L=4 control).

====================================================================
1) Lattice form of O = 1/4 tr F^2   (clover, a=1, rep: U = c0*1 + i c.sigma)
====================================================================
Plaquette  U_{mu nu}(x) = U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
Clover     C_{mu nu}(x) = (1/4)[ P_{mu nu}(x) + P_{mu nu}(x-nu)^dag
                                 + P_{mu nu}(x-mu)^dag + P_{mu nu}(x-mu-nu) ]
           (average of the 4 same-plane plaquettes containing x)
Expansion  U_{mu nu} = 1 + i a^2 F_{mu nu} + O(a^4)  =>  C = 1 + i a^2 F + O(a^4)
With  tr(T^a T^b) = delta^{ab}/2, F = F^a T^a :
            F^a_{mu nu} = (2/a^2) Re c_a(C_{mu nu}) + O(a^2)
=>          O(x) = (1/4) sum_{mu<nu} sum_{a=1,2,3} (F^a_{mu nu})^2
                 = (1/a^4) sum_{mu<nu} sum_a (Re c_a(C_{mu nu}(x)))^2      [a=1 below]
  **pointwise non-negative**, bounded (|c_a| <= 1) => all moments finite.
Time-slice operator (zero momentum, A1++ projection):
            O(t) = (1/V_s) sum_{x_s} O(x_s, t),      V_s = L^3
Magnetic part only:
            O_B(t) = (1/V_s) sum_{x_s} sum_{i<j} sum_a (Re c_a(C_ij))^2
--------------------------------------------------------------------
2) Correlator / extraction (periodic L_t)
   C(t) = <O(t)O(0)> - <O><O> ; single state:
   C(t) = f1^2 (e^{-E t} + e^{-E(L_t-t)})/(1+e^{-L_t E})
   => C(0) = f1^2 ,  C(L_t/2) = f1^2/cosh(E L_t/2)
   L_t = 2:  f1^2 = C(0) ,  E1 = arccosh(C(0)/C(1))
   2 data / 2 params => 0 d.o.f.  (see U04-MEASF1.md, verdict (c))
--------------------------------------------------------------------
3) Statistics: block estimates, **winsorised mean (2.5%/97.5%)** as central value,
   error = 1.4826*MAD/sqrt(n_b).  Plain mean+-sigma printed only as a control.
Run:  py -3.12 u04_measure_f1sq.py [--quick]
Labels [T]/[C]/[E].   No QED.
"""
import math
import random
import sys
import time

L = 4                 # lattice size, set by --L
BETAS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.3]
THERM = 1200
MEAS = 3000
NBLK = 30
WSTEP = 0.60
SEED = 2026


# ------------------------------------------------------------------ quaternions
def qm(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + b0 * a1 + a2 * b3 - a3 * b2,
            a0 * b2 + b0 * a2 + a3 * b1 - a1 * b3,
            a0 * b3 + b0 * a3 + a1 * b2 - a2 * b1)


def qdag(a):
    return (a[0], -a[1], -a[2], -a[3])


# ------------------------------------------------------------------ lattice build
def build(Lv):
    V = Lv ** 4
    nbr = [[0] * V for _ in range(4)]
    sites = []
    for i in range(V):
        sites.append((i // Lv ** 3, (i // Lv ** 2) % Lv, (i // Lv) % Lv, i % Lv))
    for i, x in enumerate(sites):
        for mu in range(4):
            y = list(x)
            y[mu] = (y[mu] + 1) % Lv
            nbr[mu][i] = (y[0] * Lv + y[1]) * Lv * Lv + y[2] * Lv + y[3]
    Q = [[1.0, 0.0, 0.0, 0.0] for _ in range(4 * V)]
    lk = lambda i, mu: i * 4 + mu
    plaqs, inc = [], [[] for _ in range(4 * V)]
    for mu in range(4):
        for nu in range(mu + 1, 4):
            for i in range(V):
                j = nbr[mu][i]
                k = nbr[nu][i]
                ks = (lk(i, mu), lk(j, nu), lk(k, mu), lk(i, nu))
                plaqs.append((mu, nu, i, ks))
                for pos in range(4):
                    inc[ks[pos]].append((len(plaqs) - 1, pos))
    return V, sites, nbr, Q, plaqs, inc, lk


V, SITES, NBR, Q, PLAQS, INC, LK = build(L)
SPINS = (1, 1, -1, -1)          # +1: link enters as U ; -1: as U^dag (positions 2,3)
NP = len(PLAQS)


def staple(ks, pos):
    """cyclic product of the other three links of the plaquette, dagger per spin"""
    out = None
    for d in (1, 2, 3):
        j = (pos + d) % 4
        q = Q[ks[j]]
        if SPINS[j] < 0:
            q = qdag(q)
        out = q if out is None else qm(out, q)
    return out


def a0_of(u, S, sp):
    """1/2 Re tr(P) with the updated link u (spin sp) and staple S."""
    return (u[0] * S[0] - sp * (u[1] * S[1] + u[2] * S[2] + u[3] * S[3]))


def sweep(beta, rng):
    NL = len(Q)
    order = list(range(NL))
    rng.shuffle(order)
    gauss = rng.gauss
    rnd = rng.random
    exp = math.exp
    sq = math.sqrt
    acc = 0
    for k in order:
        u = Q[k]
        st = []
        for (p, pos) in INC[k]:
            st.append((staple(PLAQS[p][3], pos), SPINS[pos]))
        v0 = u[0] + gauss(0.0, WSTEP)
        v1 = u[1] + gauss(0.0, WSTEP)
        v2 = u[2] + gauss(0.0, WSTEP)
        v3 = u[3] + gauss(0.0, WSTEP)
        nrm = sq(v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3)
        if nrm < 1e-12:
            continue
        v = (v0 / nrm, v1 / nrm, v2 / nrm, v3 / nrm)
        d = 0.0
        for (S, sp) in st:
            d += a0_of(v, S, sp) - a0_of(u, S, sp)
        if d * beta > 0.0 or rnd() < exp(d * beta):
            u[0], u[1], u[2], u[3] = v
            acc += 1
    return acc / NL


# ------------------------------------------------------------------ measurement
def measure():
    P = [None] * NP
    a0s = [0.0] * NP
    for p in range(NP):
        ks = PLAQS[p][3]
        q = qm(qm(qm(Q[ks[0]], Q[ks[1]]), qdag(Q[ks[2]])), qdag(Q[ks[3]]))
        P[p] = q
        a0s[p] = q[0]
    idx = {}
    for p, (mu, nu, i, _ks) in enumerate(PLAQS):
        idx[(mu, nu, i)] = p
    ntsl = V // L
    full = [0.0] * L
    magn = [0.0] * L
    for i, x in enumerate(SITES):
        for mu in range(4):
            for nu in range(mu + 1, 4):
                c1 = c2 = c3 = 0.0
                for su, sv in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    y = list(x)
                    if su < 0:
                        y[mu] = (y[mu] - 1) % L
                    if sv < 0:
                        y[nu] = (y[nu] - 1) % L
                    j = (y[0] * L + y[1]) * L * L + y[2] * L + y[3]
                    q = P[idx[(mu, nu, j)]]
                    if su * sv > 0:
                        c1 += q[1]; c2 += q[2]; c3 += q[3]
                    else:                       # reversed traversal = dagger
                        c1 -= q[1]; c2 -= q[2]; c3 -= q[3]
                c1 *= 0.25; c2 *= 0.25; c3 *= 0.25
                val = c1 * c1 + c2 * c2 + c3 * c3
                t = x[3]
                full[t] += val
                if nu < 3:               # both indices spatial (time direction = 3)
                    magn[t] += val
    up = sum(a0s) / NP
    return [full[t] / ntsl for t in range(L)], [magn[t] / ntsl for t in range(L)], up


# ------------------------------------------------------------------ robust stats
def rob(v, frac=0.025):
    n = len(v)
    s = sorted(v)
    k = max(0, int(frac * n))
    lo, hi = s[k], s[n - 1 - k]
    w = [min(max(x, lo), hi) for x in v]
    m = sum(w) / n
    med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    mad = sorted(abs(x - med) for x in v)[n // 2]
    return m, 1.4826 * mad / math.sqrt(n), med, lo, hi


def naive(v):
    n = len(v)
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / max(1, n - 1))
    return m, sd / math.sqrt(n)


def run(beta, seed, therm, meas, nblk, verbose=True):
    rng = random.Random(seed)
    for i in range(len(Q)):
        a = [rng.gauss(0, 1) for _ in range(4)]
        n = math.sqrt(sum(x * x for x in a))
        Q[i][0], Q[i][1], Q[i][2], Q[i][3] = (x / n for x in a)
    acc = 0.0
    for _ in range(therm):
        acc += sweep(beta, rng)
    bs = max(1, meas // nblk)
    keys = ("O0", "O1", "O0O1", "O0sq", "B0", "B1", "B0B1", "B0sq", "up")
    blocks = {k: [] for k in keys}
    cur = {k: 0.0 for k in keys}
    n = 0
    cnt = 0
    for _ in range(meas):
        sweep(beta, rng)
        of, ob, up = measure()
        cur["O0"] += of[0]; cur["O1"] += of[1]
        cur["O0O1"] += of[0] * of[1]; cur["O0sq"] += of[0] * of[0]
        cur["B0"] += ob[0]; cur["B1"] += ob[1]
        cur["B0B1"] += ob[0] * ob[1]; cur["B0sq"] += ob[0] * ob[0]
        cur["up"] += up
        n += 1
        if n % bs == 0:
            for k in keys:
                blocks[k].append(cur[k] / bs)
                cur[k] = 0.0
    if verbose:
        print(f"    acceptance(therm) = {acc/max(1,therm):.3f}  blocks={len(blocks['O0'])}",
              flush=True)
    return blocks


def analyse(bl):
    def cen(k):
        return rob(bl[k])
    o0, e_o0, _, _, _ = cen("O0")
    o1, e_o1, _, _, _ = cen("O1")
    o0s, e_o0s, _, _, _ = cen("O0sq")
    o0o1, e_o0o1, _, _, _ = cen("O0O1")
    b0, _, _, _, _ = cen("B0")
    b1, _, _, _, _ = cen("B1")
    b0s, _, _, _, _ = cen("B0sq")
    b0b1, _, _, _, _ = cen("B0B1")
    up, eup, mdup, _, _ = cen("up")
    C0 = o0s - o0 * o0
    C1 = o0o1 - o0 * o1
    CB0 = b0s - b0 * b0
    CB1 = b0b1 - b0 * b1
    # robust errors by the delta method on the block series
    eC0 = math.hypot(e_o0s, 2 * abs(o0) * e_o0)
    eC1 = math.hypot(math.hypot(e_o0o1, abs(o0) * e_o1), abs(o1) * e_o0)
    # also block-level direct estimates (independent cross-check)
    C0b = [s - m * m for s, m in zip(bl["O0sq"], bl["O0"])]
    C1b = [p - x * y for p, x, y in zip(bl["O0O1"], bl["O0"], bl["O1"])]
    CB1b = [p - x * y for p, x, y in zip(bl["B0B1"], bl["B0"], bl["B1"])]
    r0 = rob(C0b)
    r1 = rob(C1b)
    rb1 = rob(CB1b)
    return dict(O0=o0, C0=C0, eC0=eC0, C0b=r0, C1=C1, eC1=eC1, C1b=r1,
                CB0=CB0, CB1=CB1, CB1b=rb1, up=up, naive_up=naive(bl["up"]), eup=eup)


def main():
    global L, V, SITES, NBR, Q, PLAQS, INC, NP, BETAS
    for a in sys.argv[1:]:
        if a.startswith("--L="):
            L = int(a.split("=")[1])
        if a.startswith("--betas="):
            BETAS = [float(x) for x in a.split("=")[1].split(",")]
    therm = THERM if L <= 2 else max(200, THERM // 4)
    meas = MEAS if L <= 2 else max(300, MEAS // 5)
    V, SITES, NBR, Q, PLAQS, INC, LK = build(L)
    NP = len(PLAQS)
    t0 = time.time()
    print(f"[E] SU(2) Wilson   L={L}^4  V={V}  plaq={NP}  beta={BETAS}", flush=True)
    print(f"[E] O(t)=sum_(mu<nu)sum_a(Re c_a(clover))^2/V_s ;  C(t)=<O(t)O(0)>-<O>^2", flush=True)
    print(f"[E] therm={therm} meas={meas} nblk={NBLK} step={WSTEP} seed={SEED}", flush=True)
    print(f"[E] extraction (periodic cosh, L_t={L}): f1^2=C(0), E1=(2/L)arccosh(C(0)/C(t=L/2))",
          flush=True)
    out = {}
    tt = L // 2
    for beta in BETAS:
        print(f"\n[C] ---- beta={beta} ----", flush=True)
        bl = run(beta, SEED + int(round(beta * 1000)), therm, meas, NBLK)
        r = analyse(bl)
        C0 = r["C0"]; C1 = r["C1"]
        r["C0"] = r["C0b"][0]; r["eC0"] = r["C0b"][1]
        C0 = r["C0"]
        C1 = r["C1b"][0]; eC1 = r["C1b"][1]
        r["C1"] = C1; r["eC1"] = eC1
        CB1 = r["CB1b"][0]
        print(f"    <O>        = {r['O0']:.6e}", flush=True)
        print(f"    C(0)       = {C0:.6e}  +- {r['eC0']:.2e}   (winsor/blk)", flush=True)
        print(f"    C({tt})       = {C1:.6e}  +- {eC1:.2e}   (winsor/blk)", flush=True)
        print(f"    C({tt})/C(0) = {(C1/C0 if C0 != 0 else float('nan')):.6f}", flush=True)
        print(f"    C_B(0)={r['CB0']:.6e}  C_B({tt})={CB1:.6e} +- {r['CB1b'][1]:.2e}"
              f"  [magnetic only; beta=0 acid test]", flush=True)
        if C1 > 0 and C0 > C1:
            ratio = C0 / C1
            E1 = (2.0 / L) * math.acosh(ratio)
            eE1 = (2.0 / L) / math.sqrt(max(ratio ** 2 - 1, 1e-300)) * (
                r["eC0"] / C0 + eC1 / C1)
        else:
            ratio = float("inf"); E1 = float("nan"); eE1 = float("nan")
        print(f"    f1^2=C(0)  = {C0:.6e}  +- {r['eC0']:.2e}", flush=True)
        print(f"    E1         = {E1:.4f}  +- {eE1:.4f}", flush=True)
        print(f"    <1/2 Re tr U_p> = {r['up']:.6f} +- {r['eup']:.2e} (robust)"
              f" | {r['naive_up'][0]:.6f} +- {r['naive_up'][1]:.2e} (sigma, control)", flush=True)
        print(f"    [t={time.time()-t0:.1f}s]", flush=True)
        out[beta] = dict(O=r["O0"], C0=C0, eC0=r["eC0"], C1=C1, eC1=eC1,
                         CB0=r["CB0"], CB1=CB1, f1=C0, E1=E1, eE1=eE1,
                         up=r["up"], eup=r["eup"], ratio=(C1 / C0 if C0 else float("nan")))
    print(f"\n=== SUMMARY  L={L} ===", flush=True)
    hdr = "beta     <O>            C(0)           C(1)          C1/C0      E1      f1^2=C(0)      CB(1)"
    print(hdr, flush=True)
    for beta in BETAS:
        r = out[beta]
        print(f"{beta:4.1f}  {r['O']:.6e}  {r['C0']:.6e}  {r['C1']:.6e}  {r['ratio']:9.6f}"
              f"  {r['E1']:7.4f}  {r['f1']:.6e}  {r['CB1']:.6e}", flush=True)
    print(f"\n[E] total wall {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
