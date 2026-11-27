"""gm5_ps.py -- GM5-PS witness: exact SU(2) character activity + KP/PS thresholds.

Question: what beta_A threshold does a Pirogov-Sinai / contour (KP-Ursell) criterion
actually reach, in the repo's own currency (lambda_a = e*Delta, Delta=20, d=4)?

Deterministic, stdlib + scipy.special.iv only.  Prints a log; no MC, no downloads.
"""
from fractions import Fraction
from scipy.special import iv

DELTA4 = 4 * (2 * (4 - 1) - 1)   # = 20, plaquette adjacency degree in d=4
LAM = 2.718281828459045 * DELTA4  # e*Delta = 54.365636...


def c(n, b):
    return iv(n, b) - iv(n + 2, b)


def u(b):
    return c(1, b) / c(0, b)


def phi(a, b):
    s = 0.0
    c0 = c(0, b)
    for n in range(1, 40):
        s += (c(n, b) / c0) * (n + 1) * (2.718281828459045 ** (a * n))
    return s


def rho(a, b):
    return LAM * phi(a, b)


def f0(b):
    return c(0, b) / (2.718281828459045 ** b)


print("== [0] normalisation check: sum_n (n+1) c_n = e^beta (SU(2)) ==")
for b in (1.0 / 48, 0.025, 0.625, 2.5):
    s = sum((n + 1) * c(n, b) for n in range(60))
    import math
    print(f"  beta={b:.9f}  sum={s:.10f}  e^b={math.exp(b):.10f}  relerr={abs(s-math.exp(b))/math.exp(b):.2e}")

print()
print("== [1] exact activity / diluteness / Ursell ratio at the four points ==")
print("   beta_A        u=c1/c0     f0=c0/e^b     Phi(a=1)     rho_Ur(d=4)   u*LAM")
for b in (3.777034e-3, 1.0 / 48, 0.025, 0.037, 0.1, 0.625, 2.5):
    print(f"  {b:.9f}  {u(b):.6f}   {f0(b):.6f}    {phi(1.0,b):.4e}   {rho(1.0,b):.6e}  {u(b)*LAM:.4f}")

print()
print("== [2] window W(beta) from the envelope criterion  e*Delta*Phi(a,beta) <= a/(1+a) ==")
print("   scan a in (0,3], solve for beta by bisection on  G(beta;a)=LAM*Phi(a,beta)-a/(1+a) ==")


def G(b, a):
    return LAM * phi(a, b) - a / (1.0 + a)


def beta_root(a, lo=1e-9, hi=0.30):
    glo, ghi = G(lo, a), G(hi, a)
    if glo > 0 or ghi < 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if G(mid, a) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


best = (0.0, None)
rows = []
a = 0.02
while a <= 3.0001:
    r = beta_root(a)
    if r is not None:
        rows.append((a, r))
        if r > best[0]:
            best = (r, a)
    a = round(a + 0.01, 10)
for a, r in rows[::10]:
    print(f"   a={a:.2f}   W={r:.6e}")
print(f"  MAX over a: W*={best[0]:.6e} at a*={best[1]:.3f}   (repo: 3.777034e-3 at a*=0.6169)")

print()
print("== [3] Adhikari-Cao weak-coupling PS threshold beta >= (114 + 4 ln|G|)/Delta_G ==")
print("   (arXiv:2202.10375v3, Thm 1.1, read verbatim from the HTML body: 114, not 11/4)")
import math
for nG, name in ((2, "|G|=2"), (6, "|G|=6 (S3)"), (24, "|G|=24")):
    thr = (114.0 + 4 * math.log(nG))
    print(f"  {name}: numerator = {thr:.4f}")
    print(f"     needed Delta_G at beta_A=1/48 = {thr*48:.2f} ; at beta_A=0.625 = {thr/0.625:.2f}")
print("  SU(N): Delta_G = min_{g!=1} Re[chi(1)-chi(g)] = 0  =>  threshold = +infinity, rate (beta/2)*Delta_G = 0")
print()
print("== [4] large-beta PS: correlation length xi = 2/(beta*Delta_G) (frozen phase) ==")
for b in (1.0 / 48, 0.625, 2.5, 6.0):
    print(f"  beta_A={b:.6f}: xi*Delta_G = {2.0/b:.4f}  -> xi -> 0 as beta -> inf (opposite of AF: m_lat ~ Lambda*e^(-c beta))")

print()
print("== [5] crossover-window Peierls numbers (repo, beta_W convention) ==")
print("  sigma(2.5)=0.5701 ; 0-form threshold c0 in [log3,log7] = [1.0986,1.9459] ; 1-form log5=1.6094")
print(f"  sigma-c0 : {0.5701-1.0986:+.4f} (best 0-form) ; {0.5701-1.9459:+.4f} (worst 0-form) ; {0.5701-1.6094:+.4f} (1-form)")
print(f"  W(window)/0.625 = {3.777034e-3/0.625:.4e} ; (1/48)/3.777034e-3 = {(1/48)/3.777034e-3:.4f}")
print(f"  u(0.625)=u(beta_W=2.5) = {u(2.5):.6f} ; LAM*u = {LAM*u(2.5):.4f}")
print(f"  u(0.025)*LAM = {LAM*u(0.025):.4f}  (needs <1)")
print(f"  u(0.037)*LAM = {LAM*u(0.037):.4f}  (needs <1)")
