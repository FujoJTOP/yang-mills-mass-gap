import mpmath as mp
mp.mp.dps = 30
I = mp.besseli
def c0(b): return 2*I(1,b)/b
def c1(b): return I(1,b)-I(3,b)
def u(b): return c1(b)/c0(b)          # u_P3 = c1/c0 = 2 I2/I1
def ut(b): return I(2,b)/I(1,b)       # u_tile (T5 norm)
def Phi(b): return b*mp.e**b/(2*I(1,b))-1
print("== A. u(beta) values ==")
for b in ["0.0183","0.018268566","0.0183939720585721","0.0208","0.0208333333333333","2.5"]:
    bb=mp.mpf(b); print(f"  beta={float(bb):.10f}  u={mp.nstr(u(bb),10)}  u_tile={mp.nstr(ut(bb),10)}  Phi0={mp.nstr(Phi(bb),10)}")
print("== B. roots ==")
rt = mp.findroot(lambda b: u(b)-1, 2.4); print("  u(beta)=1  ->", mp.nstr(rt,12))
rt2 = mp.findroot(lambda b: u(b)-mp.mpf(1)/2, 1.0); print("  u=1/2      ->", mp.nstr(rt2,12))
rt3 = mp.findroot(lambda b: Phi(b)-1, 0.7); print("  Phi0=1     ->", mp.nstr(rt3,12))
ED = 20*mp.e
for lam,nm in [(ED,"eD=54.3656"),(mp.mpf("53.0006865328"),"LP=53.0007"),(48,"48"),(mp.mpf("25.2"),"25.2"),(1,"1")]:
    r = mp.findroot(lambda b: lam*Phi(b)-1, mp.mpf("0.018"))
    print(f"  lam*Phi0=1, lam={nm:12s} -> beta_max={mp.nstr(r,12)}  x(1/48)={mp.nstr(r*48,8)}")
print("== C. lambda needed for a target beta ==")
for b in ["0.0208","0.0208333333333333","0.625","0.5"]:
    bb=mp.mpf(b); print(f"  beta={float(bb):.8f}: lam_max=1/Phi0={mp.nstr(1/Phi(bb),10)}   (eD/lam={mp.nstr(ED/(1/Phi(bb)),8)})")
print("== D. the 1.137x gap ==")
print("  1/48 / 0.0182686      =", mp.nstr((mp.mpf(1)/48)/mp.mpf("0.018268566"),10))
print("  1/48 / 0.0183939721   =", mp.nstr((mp.mpf(1)/48)/mp.mpf("0.0183939720585721"),10), " = eD/48 ?", mp.nstr(ED/48,10))
print("  0.0208/0.0183        =", mp.nstr(mp.mpf("0.0208")/mp.mpf("0.0183"),10))
print("  0.0208/0.0182686     =", mp.nstr(mp.mpf("0.0208")/mp.mpf("0.018268566"),10))
print("  8/7.06328527049      =", mp.nstr(8/mp.mpf("7.06328527049169"),10), " 2N^2*48/eD =", mp.nstr(2*4*48/ED,12))
print("== E. 2D acid test (task normalisation m=3g^2/8, sigma=g^2/2) ==")
g=mp.mpf(1)
print("  m a = 3g^2/8 =", mp.nstr(3*g**2/8,8), " 4 sigma a^2 = 4*(g^2/2) =", mp.nstr(4*g**2/2,8), " ratio m/(4 sigma a) =", mp.nstr((3*g**2/8)/(4*g**2/2),8))
print("== F. misc ==")
print("  u(0.0183)/1 :", mp.nstr(u(mp.mpf("0.0183")),10), " => 1-u =", mp.nstr(1-u(mp.mpf("0.0183")),10))
print("  1/(2 eD)   =", mp.nstr(1/(2*ED),10), " beta(u=1/(2eD)) =", mp.nstr(mp.findroot(lambda b: u(b)-1/(2*ED), mp.mpf("0.018")),10))
print("  u(2.5)=    ", mp.nstr(u(mp.mpf("2.5")),12))
print("  4*ln(1/u) at 1/48 (u_P3) =", mp.nstr(4*mp.log(1/u(mp.mpf(1)/48)),10))
print("  4*ln(1/u_tile) at 1/48   =", mp.nstr(4*mp.log(1/ut(mp.mpf(1)/48)),10))
print("  sigma a^2 = ln(1/u_P3)(1/48) =", mp.nstr(mp.log(1/u(mp.mpf(1)/48)),10))
print("  u_tile max over [0,10]: sup =", mp.nstr(max(ut(mp.mpf(k)/100) for k in range(1,1001)),10))
