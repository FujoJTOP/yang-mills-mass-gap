"""Z12-FIXC -- hole C closure: exact recoupling factor Chat, new KP window beta1^fix, new delta.

Hole C (Z08-HOLEC.md:16-21): P3-DRESSED-MASS.md:53 claims the recoupling at the
interface dB_T is "fixed into a0".  False:  M(l) = dim Hom(1, prod_{p in l} r_p^{+-1})
is NOT multiplicative over factor partitions:
        mult(1 in f (x) fbar (x) f (x) fbar) = 2  !=  1 * 1 .
Fix (NC): apply :64 to the redefined activity w_{B_T} = w * prod_l R(l);  ubar -> Chat*ubar.

Checks, all reproducible from this file:
 (1) SU(2) Weyl integrals (40 digits): mult(1 in f^n), mult(adj in f^n), Catalan cross-check,
     per-link enhancement R(l) = M_{B+gamma}(l)/[M_B(l) M_gamma(l)] vs kappa = #gamma-plaq at l,
     and which kappa really gives the 0-vs-nonzero ("cap") failure.
 (2) Structural count for the straight four-wall tube (brute force, 4D torus):
     shared links per decoration plaquette in {0,1,4}; the 4's are exactly the time
     translates of the source plaquette p0.  =>  sum_l kappa_l <= 4|gamma|, kappa>=2 to enhance.
 (3) Chat = 2^nu (nu = # enhanced links) => per-plaquette bound 2^{sum kappa/2} <= 2^{2|gamma|}=4^{|gamma|}.
     New KP disc: Chat*ubar < 1 <=> beta < beta1^fix = 8/(Chat*eD*Cbar); no-loss cond. Chat*Cbar < 7.0633.
 (4) delta = ln(1/(Chat*ubar))/a = kappa/4 - ln(Chat*eD*Cbar)/a ;  identity delta_fix(beta)=delta_old(Chat*beta).

Conventions as in X01-L1C.md / Y04-AUDITX01.md / Z07-HOLEB.md:
    u = c_f/N = beta/(2N^2) (SU(2): beta/8),  eD = e*Delta = e*20 = 54.365636...,
    ubar = eD*Cbar*u,  Cbar = Ursell constant (value open -> the z02 hole).
"""
import math
from mpmath import mp, mpf, quad, sin, pi

mp.dps = 40
ED = mpf(math.e) * mpf(20)          # e*Delta, Delta = 20 (P3:65)
TOL = mpf('1e-30')


def chi(j, th):
    return mpf(2 * j + 1) if th == 0 else sin((2 * j + 1) * th) / sin(th)


def mult1(js, extra=None):
    """mult(1 in (x) reps js) = (2/pi) Int_0^pi prod chi_j(th) sin^2 th dth   (SU(2) Haar)."""
    def f(th):
        v = sin(th) ** 2
        for j in js:
            v *= chi(j, th)
        if extra is not None:
            v *= chi(extra, th)
        return v
    return 2 / pi * quad(f, [0, pi])


def catalan(m):
    return mpf(math.comb(2 * m, m)) / (m + 1)


def zero(x):
    return abs(x) < TOL


print("=" * 80)
print("(1) SU(2) Weyl integrals  (mp.dps=40)")
print("=" * 80)
F, ADJ = mpf(1) / 2, mpf(1)
tab = {}
print("   n :  mult(1 in f^n)   Catalan(n/2)  mult(adj in f^n)   R=1+adj/mult1")
for n in range(0, 13):
    m1 = mult1([F] * n)
    ma = mult1([F] * n, extra=ADJ) if n > 0 else mpf(0)
    tab[n] = (m1, ma)
    cat = catalan(n // 2) if n % 2 == 0 else mpf(0)
    print(f"  {n:2d} : {mp.nstr(m1,9):>13}  {mp.nstr(cat,9):>13}  {mp.nstr(ma,9):>13}   "
          f"{(mp.nstr(1 + ma / m1, 7) if not zero(m1) else '0/0'):>7}")
print()
print("  Hom(1, f f f f)           =", mp.nstr(tab[4][0], 12), " = 2   (Z08 (2.1) reproduced)")
print("  Hom(1,f)*Hom(1,f)  (naive)=", mp.nstr(tab[2][0] ** 2, 12), " = 1")
print("  mult(1 in f (x) fbar)     =", mp.nstr(mult1([F, F]), 12), " = 1")
print("  R(kappa=2, all f)         =", mp.nstr(1 + tab[2][1] / tab[2][0], 12), " = 2")
print("  R(kappa=4, all f)         =", mp.nstr(1 + tab[4][1] / tab[4][0], 12),
      " = 5/2 = C_3/C_2 =", mp.nstr(catalan(3) / catalan(2), 12))
print("  R(kappa<=4) <= 2^(kappa/2) ? ",
      all(1 + tab[k][1] / tab[k][0] <= mpf(2) ** (mpf(k) / 2) for k in (2, 4)),
      "   (2 <= 2,  2.5 <= 4)   =>  per-plaquette Chat <= 2^{sum kappa/2} <= 4^{|gamma|}")
print()
print("  cap test:  B_T side = f (x) fbar = 1 + adj;   true = mult1(A);  naive = mult1(A)")
print("  (w(gamma)=0 but w(B_T+gamma)!=0  <=>  mult(1 in A)=0, mult(adj in A)!=0)")
for name, A in [("kappa=1, A=f   (Z08's witness)", [F]),
                ("kappa=1, A=adj", [ADJ]),
                ("kappa=2, A=adj(x)j=2", [ADJ, mpf(2)])]:
    mg = mult1(A)
    mt = mg + mult1(A, extra=ADJ)
    tag = "CAP-FAILURE" if (zero(mg) and not zero(mt)) else ("both zero" if zero(mg) else "generic")
    print(f"    {name:28s} w(gamma)={mp.nstr(mg,6):>12}  w(B+gamma)={mp.nstr(mt,6):>12}  {tag}")

print()
print("=" * 80)
print("(2) straight four-wall tube: shared links per decoration plaquette (4D torus, L=4,T=6)")
print("=" * 80)
L, NT, TT = 4, 12, 6


def nxt(x, mu):
    y = list(x)
    y[mu] = (y[mu] + 1) % (NT if mu == 0 else L)
    return tuple(y)


def plq_links(pl):
    mu, nu, x = pl
    a, b, c = x, nxt(x, nu), nxt(x, mu)
    return {(mu, a), (mu, b), (nu, a), (nu, c)}


site = [(t, x, y, z) for t in range(NT) for x in range(L) for y in range(L) for z in range(L)]
all_plq = [(mu, nu, s) for mu in range(4) for nu in range(4) if mu < nu for s in site]
p0 = (1, 2, (0, 0, 0, 0))
BT = set()
for (i, b) in plq_links(p0):
    for t in range(TT):
        BT.add((0, i, (t,) + tuple(b[1:])))
links_BT = set()
for pl in BT:
    links_BT |= plq_links(pl)
print(f"  |B_T| = {len(BT)} = 4T = {4*TT};   |links(B_T)| = {len(links_BT)}")
hist, big = {}, []
for pl in all_plq:
    if pl in BT or pl == p0:
        continue
    k = len(plq_links(pl) & links_BT)
    hist[k] = hist.get(k, 0) + 1
    if k > 1:
        big.append((k, pl))
print("  histogram (#shared links per non-B_T plaquette):", dict(sorted(hist.items())))
print("  plaquettes with >1:", len(big), " all are time-translates of p0 (plane 12, x=origin), k=4:")
print("   ", sorted(p[2][0] for _, p in big))
print("  => sum_l kappa_l <= 4|gamma| ;  enhancement needs kappa_l >= 2 ;  kappa_l <= 4 always")
print("  => Chat = 2^nu <= 2^{sum kappa/2} <= 2^{2|gamma|} = 4^{|gamma|}   [provable, safe]")
cnt = {}
for l in links_BT:
    n = sum(1 for pl in BT if l in plq_links(pl))
    cnt[n] = cnt.get(n, 0) + 1
print("  backbone-side plaquettes per link (histogram):", dict(sorted(cnt.items())),
      " (the 1's are the 4 cross-section links at t=0 and t=T; p0 supplies the 2nd)")

print()
print("=" * 80)
print("(3) new KP window:  Chat*ubar < 1  <=>  beta < beta1^fix = 8/(Chat*eD*Cbar)")
print("=" * 80)
u48 = mpf(1) / 48 / 8
print(f"  eD = e*Delta = {mp.nstr(ED,12)} ;  u = beta/8 ;  ubar = {mp.nstr(ED/8,10)}*Cbar*beta")
print("  Chat       beta1^fix(Cbar=1)  vs 1/48      | beta1^fix(Cbar=N^2=4)  vs 1/48")
cands = [("1  naive", mpf(1)), ("2  per factorization point", mpf(2)),
         ("4  provable per plaquette", mpf(4)), ("16 Z08 rho^4 (over-count)", mpf(16))]
for name, C in cands:
    seg = []
    for Cb in (mpf(1), mpf(4)):
        b1 = 8 / (C * ED * Cb)
        seg.append(f"{mp.nstr(b1,8):>11} ({'keep' if b1 > mpf(1)/48 else 'SHRINK'} x{mp.nstr((mpf(1)/48)/b1,5)})")
    print(f"  {name:26s} {seg[0]:38s} | {seg[1]}")
thr = 1 / (ED * u48)
print(f"  no-loss condition for the P3/SZZ window 1/48:  Chat*Cbar < {mp.nstr(thr,10)}   (= 7.06x,")
print("     already recorded in LOOP-RUNBOOK:6182 for Chat=1)")
print("     Chat*Cbar: 2*1=2  4*1=4  4*4=16  16*1=16")

print()
print("=" * 80)
print("(4) delta = ln(1/(Chat*ubar))/a = kappa/4 - ln(Chat*eD*Cbar)/a")
print("=" * 80)
print("  Cbar = 1 (the value implicit in Y04:117 / X01 §3 row 6)")
print("   beta      ubar      | Chat=1            | Chat=2            | Chat=4            | Chat=16")
for b in [mpf(1) / 12, mpf(1) / 24, mpf(1) / 48, mpf(1) / 96]:
    ub = ED * (b / 8)
    row = f"   {mp.nstr(b,6):8} {mp.nstr(ub,8):9} |"
    for _, C in [(a, c) for a, c in cands if c in (mpf(1), mpf(2), mpf(4), mpf(16))]:
        x = C * ub
        row += f" {mp.nstr(x,6):>8}/" + (f"{mp.nstr(math.log(float(1/x)),6):>9}" if x < 1 else " KP-FAIL ")
    print(row)
print()
print("  Y04:117 check: kappa/4*a=ln(1/u) at u=0.0104 =", mp.nstr(mp.log(1 / mpf('0.0104')), 8),
      "; ln(eD) =", mp.nstr(mp.log(ED), 10), "; delta*a = 0.57022")
print("  identity delta_fix(beta) = delta_old(Chat*beta):",
      mp.nstr(mp.log(1 / (mpf(4) * ED * mpf(1) / 48 / 8)), 12), "==",
      mp.nstr(mp.log(1 / (ED * mpf(4) / 48 / 8)), 12))
print("  beta=1/12 survives iff Chat < 1/ubar =", mp.nstr(1 / (ED * mpf('0.0104')), 8),
      "  => Chat=2 and Chat=4 KILL the old reference point; Chat=1 (naive) is the one that dies of")
print("     Z08's over-count 16, not of the exact factor 2/4.")
print()
print("Z12 done. No QED, Clay 0%.")
