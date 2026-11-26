#!/usr/bin/env python3
# P3-CBAR: is `c_f/u = N` an identity or a bound?  Does re-defining u buy a 2x window?
# py -3.12, mpmath only (no numpy needed).  Read-only: writes nothing but stdout.
import sys, itertools, math
from mpmath import mp, mpf, besseli, log, exp, findroot, sqrt

mp.dps = 30
N   = mpf(2)
beta = mpf(1) / 12
DEG = 20                     # 4D: plaquettes per link adjacency, deg = 4(6-1) = 20
E_DELTA = mp.e * DEG         # e*Delta = 54.3656365692...
I48 = mpf(1) / 48

def H(t): print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)

# ---------------------------------------------------------------- A
H("[A] c_f/u = N : definitional identity?  (species + exact SU(2))")
print(f"e*Delta            = {mp.nstr(E_DELTA, 12)}")
for n in (2, 3, 4, 10):
    cf = beta / (2 * n); u = cf / n
    print(f"  N={n:<3} c_f=b/(2N)={mp.nstr(cf,10):<14} u=c_f/N={mp.nstr(u,10):<14} c_f/u={mp.nstr(cf/u,10)}")
# exact SU(2): c_f = <tr U_p> = I2(b)/I1(b)   (weight e^{(b/N)Re tr U} = e^{b cos th})
cf_ex = besseli(2, beta) / besseli(1, beta)
u_species = beta / (2 * N * N)
u_z06 = cf_ex / N
print(f"  SU(2) exact  c_f = I2/I1 = {mp.nstr(cf_ex,12)}   beta/(2N) = {mp.nstr(beta/(2*N),12)}  (rel.diff {mp.nstr(abs(cf_ex-beta/(2*N))/(beta/(2*N)),4)})")
print(f"  SU(2) exact  u = c_f/N   = {mp.nstr(u_z06,12)}   b/(2N^2)  = {mp.nstr(u_species,12)}  (rel.diff {mp.nstr(abs(u_z06-u_species)/u_species,4)})")
print(f"  => c_f/u = {mp.nstr(cf_ex/u_z06,12)}  (identity: u was DEFINED as c_f/N)")

# ---------------------------------------------------------------- B
H("[B] the OTHER slot: Phi(0)/u = 2N^2  (rep-summed plaquette activity, KP (★★))")
for n in (2, 3, 4, 10):
    u = beta / (2 * n * n)
    print(f"  N={n:<3} Phi(0)=b={mp.nstr(beta,8):<12} Phi(0)/u = {mp.nstr(beta/u,10)}   (2N^2 = {2*n*n})")

# ---------------------------------------------------------------- C
H("[C] u-renormalisation test: is (Cbar, u) -> (Cbar/alpha, alpha*u) inert?")
print("  Ubar := eD * Cbar * u = eD * sup_p |w_p|   <-- contains NO u.\n")
print(f"  {'u-def':<22}{'Cbar=c_f/u':>14}{'Ubar':>16}{'delta*a=ln(1/Ubar)':>22}{'beta1':>14}")
cf_tr = beta / (2 * N)
defs = [("c_f/(2N^2)  [species]", cf_tr / (2 * N * N)),
        ("c_f/(N c0), c0=1    ", cf_tr / N),
        ("c_f  [tr-unnorm]    ", cf_tr),
        ("c_f*2N^2            ", cf_tr * 2 * N * N)]
for name, u in defs:
    Cbar = cf_tr / u
    Ubar = E_DELTA * Cbar * u
    b1 = 2 * N / E_DELTA                        # solve eD*W(b)=1 with W=Cbar*u(b)=c_f(b)
    print(f"  {name:<22}{mp.nstr(Cbar,10):>14}{mp.nstr(Ubar,14):>16}{mp.nstr(log(1/Ubar),16):>22}{mp.nstr(b1,10):>14}")
print("  => every row has the SAME Ubar / delta / beta1: pure relocation (verdict (c)).")

# ---------------------------------------------------------------- D
H("[D] window as a function of the ONE physical input W := sup_p |w_p|")
print("  criterion: eD * W(beta) < 1  <=>  beta < beta1 ;  Cbar = W/u ; u:=c_f/N (fixed)\n")
reads = [("W = u   (Z06 per-tile, (1/N)tr)", lambda b: b / (2 * N * N), 1),
         ("W = c_f (PROOF 2' bound, tr)   ", lambda b: b / (2 * N), N),
         ("W = Phi(0)=beta (KP (★★))      ", lambda b: b, 2 * N * N)]
print(f"  {'reading':<36}{'Cbar=W/u':>10}{'Ubar@b=1/12':>15}{'beta1':>12}{'x(1/48)':>10}{'1/12 in?':>10}")
for name, W, Cb in reads:
    b1 = findroot(lambda b: E_DELTA * W(b) - 1, 0.1)
    Ub = E_DELTA * W(beta)
    print(f"  {name:<36}{mp.nstr(Cb,8):>10}{mp.nstr(Ub,10):>15}{mp.nstr(b1,10):>12}"
          f"{mp.nstr(b1/I48,8):>10}{('YES' if beta < b1 else 'no'):>10}")
print(f"  spread of window over the three readings = {(mp.nstr((2*N*N)/1,6))} x  (Cbar in [1, 2N^2])")

# ---------------------------------------------------------------- E
H("[E] exact SU(2) solves (no species truncation)")
for nm, f in [("eD*I2/I1 = 1        ", lambda b: E_DELTA * besseli(2, b) / besseli(1, b) - 1),
              ("eD*(I2/I1)/N = 1    ", lambda b: E_DELTA * besseli(2, b) / besseli(1, b) / N - 1),
              ("eD*(I2/I1)/(2N^2)=1 ", lambda b: E_DELTA * besseli(2, b) / besseli(1, b) / (2 * N * N) - 1)]:
    b1 = findroot(f, 0.05)
    print(f"  {nm} beta1 = {mp.nstr(b1,12):<16} = {mp.nstr(b1/I48,8)} x (1/48)   (1/12 in? {'YES' if beta < b1 else 'no'})")
print(f"  beta=1/12 exact:  eD*c_f          = {mp.nstr(E_DELTA*cf_ex,10)}   (>1 => fails)")
print(f"                    eD*c_f/N        = {mp.nstr(E_DELTA*cf_ex/N,10)}   (<1 => holds)")
print(f"                    eD*b            = {mp.nstr(E_DELTA*beta,10)}   (>1 => fails)")

# ---------------------------------------------------------------- F
H("[F] where the N comes from: 3-cube faces, link incidence, orientation product")
# unit 3-cube, sites (x,y,z) in {0,1}^3 ; faces as (axis, value) with oriented boundary
def cube_faces():
    F = []
    for axis in range(3):
        for val in (0, 1):
            i, j = [k for k in range(3) if k != axis]
            loop = []
            for (a, b) in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                v = [0, 0, 0]; v[axis] = val; v[i] = a; v[j] = b
                loop.append(tuple(v))
            if val == 0: loop = loop[::-1]
            F.append(loop)
    return F
def sign_incidence(loop):
    L = []
    for k in range(4):
        a, b = loop[k], loop[(k + 1) % 4]
        key = tuple(sorted([a, b]))
        L.append((key, 1 if a == key[0] else -1))
    return L
faces = cube_faces()
# brute force the sign of each face so that the surface is a closed 2-cycle
best = None
for mask in range(2 ** len(faces)):
    inc = {}
    for fi, f in enumerate(faces):
        s0 = 1 if (mask >> fi) & 1 == 0 else -1
        for key, s in sign_incidence(f):
            inc.setdefault(key, []).append(s0 * s)
    bad = sum(1 for v in inc.values() if len(v) != 2 or sum(v) != 0)
    if best is None or bad < best[0]: best = (bad, mask, inc)
bad, mask, inc = best
print(f"  faces = {len(faces)}   distinct links = {len(inc)}   incidences = {sum(len(v) for v in inc.values())}")
print(f"  min #bad links over all 64 orientation assignments: {bad}")
print(f"  oriented consistently? {'YES (=closed 2-cycle)' if bad == 0 else 'no'}")
print("  => prod_links mult(1 in f (x) fbar) = 1  => |w(S)| = c_f^|S| * 1")
for n in (1, 6):
    for nm, w in [("tr-read  c_f", cf_tr), ("Z06 read u  ", cf_tr / N)]:
        W = w**n
        print(f"   n={n}  {nm}: |w|^(1/n)/u = {mp.nstr(W**(1/n)/ (cf_tr/N),10)}")
print("  *** the ONLY difference between 'N' and '1' is the (1/N) per-tile normalisation ***")

# ---------------------------------------------------------------- G
H("[G] accounting: what a 2x window requires, and what is provably available")
Ubar_now = E_DELTA * cf_ex                       # f07 value
print(f"  f07 Ubar = eD*c_f         = {mp.nstr(Ubar_now,10)}  >1  -> beta=1/12 outside")
print(f"  to get beta=1/12 inside you need W <= {mp.nstr(1/E_DELTA,10)}  (i.e. eD*W<1)")
print(f"  exact single-plaquette W at beta=1/12:  <(1/N)trU>={mp.nstr(u_z06,10)}   <trU>={mp.nstr(cf_ex,10)}   Phi(0)={mp.nstr(beta,10)}")
print(f"  => no exact single-plaquette value below u = c_f/N exists: W_min = u = {mp.nstr(u_z06,10)}")
print("  => Cbar<1 impossible (Cbar >= |w({p})|/u);  Cbar<N <=> dropping the per-tile 1/N normalisation")
print("  done.")
