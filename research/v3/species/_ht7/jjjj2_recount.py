"""JJJJ2 — 独立复算 HHHH2-C1SUM（research/v3/species/HHHH2-C1SUM.md）。
不读它的脚本输出，只读它的公式：字面项 n^{n-2} w^{n-1}，无序项 /(n-1)! 与 /n!。
Decimal 60 位，零 MC。"""
from decimal import Decimal as D, getcontext

getcontext().prec = 60
G = D("0.90537044375382710143")  # _ht7/_c1sum_calc2.ps1:11
CHK = {"s4_full": D("0.15782298"), "s4_half": D("0.076744123"),
       "lit8": D("1.39896"), "lit12": D("125.61")}


def s4(a, full):
    t = (-a).exp()
    r = ((1 + t) / (1 - t)) if full else (1 / (1 - t))
    return r ** 4 - 1


def run(w, N=60):
    lit = u1 = u2 = D(0)
    f = D(1)
    keep = {}
    for n in range(2, N + 1):
        t = D(n) ** (n - 2) * w ** (n - 1)
        f *= n                              # n!
        lit += t
        u1 += t * n / f                     # /(n-1)!
        u2 += t / f                         # /n!
        if n in (2, 3, 4, 5, 6, 8, 12, N):
            keep[n] = (t, lit, t * n / f, u1, t / f, u2)
    return keep, lit, u1, u2


def cross(w, N=60):
    """字面累计首次越过 0.5 / G 的 n"""
    lit = D(0)
    c5 = cg = None
    for n in range(2, N + 1):
        lit += D(n) ** (n - 2) * w ** (n - 1)
        if c5 is None and lit > D("0.5"):
            c5 = n
        if cg is None and lit > G:
            cg = n
        if c5 and cg:
            break
    return c5, cg, lit


wf, wh = s4(D(4), True), s4(D(4), False)
print(f"S4full(4)={wf:.12f}  S4half(4)={wh:.12f}  (report 0.15782298 / 0.076744123)")
assert abs(wf - CHK["s4_full"]) < D("1e-8") and abs(wh - CHK["s4_half"]) < D("1e-9")

for name, w in (("W=0.157823 (full Z^4)", wf), ("W=0.076744 (half N^4)", wh)):
    keep, lit, u1, u2 = run(w)
    print(f"\n== {name} ==\n n      literal        lit.cum      /(n-1)!      cum/(n-1)!    /n!        cum/n!")
    for n, (t, c, a, ca, b, cb) in keep.items():
        print(f"{n:2d}  {t:13.6E} {c:13.6f} {a:13.6E} {ca:13.6f} {b:13.6E} {cb:13.6f}")
    print(f"SUM n<=60: literal={lit:.4E}  /(n-1)!={u1:.9f}  /n!={u2:.9f}")
    print(f"G/lit={G/lit:.4E}  G/u1={G/u1:.4f}  G/u2={G/u2:.4f}  n=2 share={w/u1:.4f}")
    print(f"C1<=G ? {u1<=G}   C1<0.5 ? {u1<D('0.5')}   (u1={u1:.6f})")
    print(f"literal crosses 0.5 at n={cross(w)[0]}, crosses G at n={cross(w)[1]}")

keep, lit, u1, u2 = run(wf)
print(f"\ndiff vs report: lit(n=8)={keep[8][1]-CHK['lit8']:+.2E}  "
      f"lit(n=12)={keep[12][1]-CHK['lit12']:+.2E}  sum/(n-1)! = {u1-D('0.210517'):+.2E}")
print(f"0.210517 reproducible from 0.157823? {abs(u1-D('0.210517'))<D('5e-7')}")
print(f"half-line sum = {run(wh)[2]:.6f} -> 0.210517? {abs(run(wh)[2]-D('0.210517'))<D('5e-7')}")
