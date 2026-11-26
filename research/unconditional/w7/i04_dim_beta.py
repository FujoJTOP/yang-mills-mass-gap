# i04_dim_beta.py — I04: f1^2 的 beta 标度 / 量纲账 / 与 m 的关系
# 只做量纲+两圈 AF 的精确算术；不得当证据用（[E] 定位为"量纲计数的一致性检查"）。
# py -3.12 -I
from mpmath import mp, mpf, exp, log, pi, power
mp.dps = 40

# --- 约定（G05-BETAOFa.md:56-61, A6-TARGET.md:43-45） ---
# beta_W = 2N/g^2 ; b0 = 11N/(48 pi^2) ; b1/(2 b0^2) = 51/121
def aL(beta_W, N):
    """a*Lambda_lat, two loops (asymptotic; [C])."""
    g2 = 2 * N / beta_W
    e = 12 * pi**2 * beta_W / (11 * N**2)          # = 1/(2 b0 g^2)
    return exp(-e) * power(g2 / (16 * pi**2), mpf(51) / 121)

def logslope(N):
    """-d ln(aL)/d beta_W  -> 12 pi^2/(11 N^2) (exponent of a)."""
    return 12 * pi**2 / (11 * N**2)

print("=" * 78)
print("I04  f1^2 = |<Omega|O|G>|^2 ,  O = (1/4) tr F^2  =>  [f1^2] = mass^8 (dim 8)")
print("     f1^2_lat(beta) = a(beta)^8 * f1^2_phys   =>  f1^2_lat ~ exp(-c' beta)")
print("=" * 78)
print()
print("[A] two-loop a*Lambda_lat and the exponent 12 pi^2/(11 N^2)")
print(f"{'N':>2} {'beta_W':>8} {'g^2':>10} {'a*L_lat':>14} {'dln(aL)/dbeta':>16} {'-12pi^2/(11N^2)':>18}")
for N in (2, 3):
    for bW in (mpf('1.0'), mpf('2.0'), mpf('2.5'), mpf('6.0'), mpf('10.0'), mpf('100.0')):
        h = mpf('1e-6')
        s = (log(aL(bW + h, N)) - log(aL(bW - h, N))) / (2 * h)
        print(f"{N:>2} {float(bW):>8.2f} {float(2*N/bW):>10.4f} {float(aL(bW,N)):>14.6e} "
              f"{float(s):>16.8f} {float(-logslope(N)):>18.8f}")
    print()

print("[B] c'(d) = d * 12 pi^2/(11 N^2)  for the three dimension conventions in the corpus")
print("    d=8 : |<Omega|O|G>|^2 with O dim 4  (ch05 reading, T02-F1SQ.md:17)")
print("    d=6 : KL-measure atom nu({m1^2})      (A2-SUMRULE-CONDENSATE.md, 'dimension-6 measure')")
print("    d=2 : ch05 'operator normalized so f1^2 has dimension 2' (A2 note)")
print(f"{'d':>3} {'N':>2} {'c\'=d*12pi^2/(11N^2)':>22} {'f1^2(2.5)/f1^2(6.0)':>22}")
for d in (2, 6, 8):
    for N in (2, 3):
        c = d * logslope(N)
        r = exp(-c * (mpf('2.5') - mpf('6.0')))
        print(f"{d:>3} {N:>2} {float(c):>22.6f} {float(r):>22.4e}")
    print()

print("[C] same driver as m, different power:  m_lat ~ a*Lambda ~ exp(-12pi^2 beta/(11N^2))")
print("    => f1^2_lat / m_lat^8  is beta-independent (dimension counting, up to logs)")
N = 3
for bW in (mpf('2.5'), mpf('6.0'), mpf('10.0')):
    a8 = aL(bW, N)**8
    m1 = aL(bW, N)                      # m_lat with m/Lambda = 1
    print(f"    N={N} beta_W={float(bW):>5.1f}: a^8={float(a8):>12.4e}  m^8={float(m1**8):>12.4e}  "
          f"a^8/m^8={float(a8/m1**8):>10.6f}")

print()
print("[D] the target as written:  exists c>0  forall beta: f1^2_lat(beta) >= c ?")
N = 3
for bW in (mpf('2.5'), mpf('6.0'), mpf('10.0'), mpf('20.0')):
    print(f"    N={N} beta_W={float(bW):>5.1f}: f1^2_lat (F=1) = {float(aL(bW,N)**8):>12.4e}   "
          f"(falls by {float(((aL(mpf('2.5'),N)/aL(bW,N))**8)):>10.3e} x w.r.t. beta=2.5)")
print("    => inf over beta is 0  [NG]; the only legitimate form is the RATIO (a^-8 f1^2 = f1^2_phys)")
print()

print("[E] window anchor (RUN3 via B04.md:204, beta=2.5, L=8, SU(3)): w1=0.3427, C(0)=19.888")
w1 = mpf('0.3427'); C0 = mpf('19.888')
f1lat = w1 * C0
print(f"    f1^2_lat(2.5) = w1*M0 = {float(f1lat):.4f}   (lattice units, O(1) -- NOT a^-8 suppressed)")
print(f"    (a*L)^8 at beta=2.5 = {float(aL(mpf('2.5'),3)**8):.3e}  =>  the two differ by "
      f"{float(f1lat/aL(mpf('2.5'),3)**8):.3e} x")
print("    => the a^8 law is an ASYMPTOTIC ([C], AF-window) statement; the accessible window")
print("       (beta<=2.5) carries contact/tadpole contamination and cannot confirm it.")
print()
print("[F] strong-coupling floor r0 = (1/4)(2N)^{-2} vs the beta=2.5 lattice value")
for N in (2, 3):
    r0 = mpf(1) / (4 * (2 * N)**2)
    print(f"    N={N}: r0 = {float(r0):.6f} ; f1^2_lat(2.5)={float(f1lat):.4f} "
          f"=> ratio {float(f1lat/r0):.1f}  (window reading: 'floor free' -- but it is a MIXING/TADPOLE number)")
print()
print("[G] verdict summary")
print("    f1^2 dim = 8 (mass^8) => lattice beta-scaling = a(beta)^8 = exp(-8*12pi^2 beta/(11N^2))")
print("    SU(3): c'=%.4f ; SU(2): c'=%.4f" % (float(8*logslope(3)), float(8*logslope(2))))
print("    => branch (2) of the task: c(beta) ~ exp(-c' beta), NOT beta-independent (1), NOT 1/beta^n (3)")
print("    => f1^2 ~ m^8 : same AF driver a(beta), power 8 vs 1")
