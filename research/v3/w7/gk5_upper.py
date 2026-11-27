"""GK5-UPPER — window f(lambda_*) = Phi(0,.)^{-1}(1/lambda_*), upper end beta=0.625.

Objects:
  (1) exact threshold lambda_* <= 1/Phi(0,beta_target) for each target beta
  (2) class ceiling over ALL lambda_a >= 1  (Phi(0,.)^{-1}(1) = 0.765528217952)
  (3) beta_A / beta_W coordinate audit (Phi's Bessel argument = beta_W)
  (4) a *rigorous* animal-counting lower bound lambda_* >= 2^(2/3) = 1.5874
      (comb construction: rod of m plaquettes + arbitrary subset decorated)
Reproduce with:  py -3.12 -I research/v3/w7/gk5_upper.py
"""
import mpmath as mp
mp.mp.dps = 30
I = mp.besseli

def Phi(b):                      # UC-COMPLETE:41  Phi(0,b) = b e^b/(2 I1(b)) - 1
    return b * mp.e**b / (2 * I(1, b)) - 1

def beta_max(lam):               # window root: lam * Phi(0,b) = 1   (Phi strictly increasing)
    return mp.findroot(lambda b: lam * Phi(b) - 1, mp.mpf("0.02") if lam > 2 else mp.mpf("0.4"))

ED = 20 * mp.e
Lt = lambda D: mp.mpf(D)**D / mp.mpf(D - 1)**(D - 1)

P = print
P("== [1] thresholds for a target beta (exact form: lam <= 1/Phi(0,beta)) ==")
targets = [("beta_A=1/48 (=beta_W 1/12)", mp.mpf(1) / 48),
           ("beta_A=0.625 (=beta_W 2.5)", mp.mpf("0.625")),
           ("beta_W=1/12", mp.mpf(1) / 12),
           ("beta_W=2.5", mp.mpf("2.5")),
           ("beta_W=0.625", mp.mpf("0.625"))]
for nm, b in targets:
    P(f"  {nm:32s} Phi(0,beta)={mp.nstr(Phi(b),12):18s}  lam_max={mp.nstr(1/Phi(b),12)}")

P("\n== [2] window f(lam) at the constants actually on the table ==")
rows = [("1 (structural floor, [T])", mp.mpf(1)),
        ("2^(2/3)=1.5874 (this round, [T])", mp.mpf(2)**(mp.mpf(2)/3)),
        ("1/Phi(0,0.625)=1.2821989 (GJ1's number)", 1/Phi(mp.mpf("0.625"))),
        ("a_6^(1/6)=13.902 ([V])", mp.mpf("13.902")),
        ("lam_* ~ 25.2 ([V])", mp.mpf("25.2")),
        ("1/Phi(0,1/48)=47.6270568", 1/Phi(mp.mpf(1)/48)),
        ("eD=20e=54.3656366", ED),
        ("Lt(14)=36.6884124", Lt(14))]
for nm, lam in rows:
    b = beta_max(lam)
    P(f"  lam={nm:44s} beta_max(W)={mp.nstr(b,10):14s} = {mp.nstr(b/4,10):14s} in beta_A"
      f"   x(1/48)_A={mp.nstr(b/4*48,8):10s} 0.625_A/beta_max_A={mp.nstr(mp.mpf('0.625')/(b/4),8)}")

P("\n== [3] class ceiling (lam_a = 1, no animal constant can do better) ==")
ceilW = beta_max(mp.mpf(1))
P(f"  Phi(0,.)^-1(1) = {mp.nstr(ceilW,12)} in beta_W = {mp.nstr(ceilW/4,12)} in beta_A")
P(f"  vs right end 0.625_A : short by {mp.nstr(mp.mpf('0.625')/(ceilW/4),10)}x")
P(f"  vs left  end 1/48_A : longer by {mp.nstr((ceilW/4)*48,10)}x")
P(f"  class coverage of the target window: ({mp.nstr(ceilW/4,10)}-0.0208333)/(0.625-0.0208333) = "
  f"{mp.nstr((ceilW/4-mp.mpf(1)/48)/(mp.mpf('0.625')-mp.mpf(1)/48),10)}")
P(f"  unreachable upper segment [beta_A_ceil, 0.625] length fraction = "
  f"{mp.nstr((mp.mpf('0.625')-ceilW/4)/(mp.mpf('0.625')-mp.mpf(1)/48),10)}")

P("\n== [4] coordinate audit: Phi's argument is beta_W (Bessel weight e^{beta_W cos th}) ==")
P(f"  u(beta)=b(I1-I3)/(2 I1):  u(2.5)={mp.nstr(mp.mpf('2.5')*(I(1,mp.mpf('2.5'))-I(3,mp.mpf('2.5')))/(2*I(1,mp.mpf('2.5'))),10)}"
  f"   (=corpus 1.0143902001)")
P(f"  same formula at 1.25 (if arg were 2*beta_A at the physical point): "
  f"{mp.nstr(mp.mpf('1.25')*(I(1,mp.mpf('1.25'))-I(3,mp.mpf('1.25')))/(2*I(1,mp.mpf('1.25'))),10)}")
P(f"  0.0183 (window, beta_W) in beta_A = {mp.nstr(mp.mpf('0.0183')/4,10)};  1/192={mp.nstr(mp.mpf(1)/192,10)}")
P(f"  4 * (1/48) = {mp.nstr(4*mp.mpf(1)/48,10)} = 1/12 ; 4*0.625 = 2.5")

P("\n== [5] rigorous animal bound lambda_* >= 2^(2/3) (comb construction) ==")
c = mp.mpf(2)**(mp.mpf(2)/3)
P(f"  2^(2/3) = {mp.nstr(c,12)}")
worst = mp.mpf(1)
for m in range(2, 4001):
    j = m // 2
    lhs = mp.binomial(m, j)
    rhs = mp.mpf(2)**m / (m + 1)
    if lhs < rhs:
        worst = min(worst, lhs / rhs)
P(f"  min over m<=4000 of C(m,floor(m/2))/(2^m/(m+1)) = {mp.nstr(worst,10)}  (>0 => C>=2^m/(m+1) holds)")
for m in [10, 100, 1000, 10**6]:
    n = mp.mpf(3) * m / 2
    lhs = (mp.mpf(2)**(mp.mpf(2)/3) / (mp.mpf(3)*m/2 + 1)**(mp.mpf(1)/n))
    P(f"  m={m:8d}: a_n^(1/n) >= {mp.nstr(lhs,10)}   (n=3m/2={mp.nstr(n,6)})")
P(f"  needed lam_max at 0.625_A (=2.5_W): {mp.nstr(1/Phi(mp.mpf('2.5')),10)}"
  f"  vs floor 1 => short {mp.nstr(mp.mpf(1)/(1/Phi(mp.mpf('2.5'))),10)}x ;"
  f" vs 2^(2/3) => short {mp.nstr(c/(1/Phi(mp.mpf('2.5'))),10)}x")
P(f"  lam_max at 1/48_A (=1/12_W): {mp.nstr(1/Phi(mp.mpf(1)/12),10)} vs 13.902 => short "
  f"{mp.nstr(mp.mpf('13.902')/(1/Phi(mp.mpf(1)/12)),8)}x ; vs eD => short {mp.nstr(ED/(1/Phi(mp.mpf(1)/12)),8)}x")

P("\n== [6] corpus a_n check (TN1:23 / LC4:12) ==")
a = [1, 20, 458, 11132, 280175, 7218768]
P("  a_n roots:", [mp.nstr(mp.mpf(a[i])**(mp.mpf(1)/(i+1)), 8) for i in range(6)])
P(f"  a_6^(1/6)={mp.nstr(mp.mpf(a[5])**(mp.mpf(1)/6),8)}; 13.902/(1/Phi(2.5))={mp.nstr(mp.mpf('13.902')/(1/Phi(mp.mpf('2.5'))),8)}")
P(f"  1/Phi(0,0.625)={mp.nstr(1/Phi(mp.mpf('0.625')),12)}; 1/Phi(0,2.5)={mp.nstr(1/Phi(mp.mpf('2.5')),12)}")
P(f"  ratio (1/Phi(0.625))/(1/Phi(2.5)) = {mp.nstr((1/Phi(mp.mpf('0.625')))/(1/Phi(mp.mpf('2.5'))),10)}")
P("\n== [7] 2D acid test: does 'm propto Lambda' hold in 2D? ==")
P("  2D: m=3g^2/8, sigma=g^2/2, m/sqrt(sigma)=0.75 sqrt(g^2) -> depends on g (no transmutation, no pure number)")
P("  (B06.md:73 / B04.md:134-139: a=2N/(beta g^2) identity, no dimensional transmutation in 2D/3D)")
