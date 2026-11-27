# GJ2-RAT1 witness: is lambda_* <= L_t(14) = 36.6884124 a theorem or an estimate?
# py -3.12 -I ; mpmath dps=40 ; all sections <= 5 s
from mpmath import mp, mpf, besseli, exp, findroot, log, nstr
from math import comb
mp.dps = 40

def Lt(D):                      # L_t(D) = D^D/(D-1)^(D-1)
    D = mpf(D)
    return D**D / (D-1)**(D-1)

def Phi0(b):                    # Phi(0,beta) = beta e^beta/(2 I1(beta)) - 1
    b = mpf(b)
    return b*exp(b)/(2*besseli(1, b)) - 1

def bmax(lam):                  # beta_max = Phi(0,.)^{-1}(1/lam)
    return findroot(lambda b: Phi0(b) - 1/mpf(lam), mpf('0.026'))

print("== A. L_t(D) table (D = 14..20), dps=40 ==")
for D in range(14, 21):
    v = Lt(D)
    print(f"L_t({D:2d}) = {nstr(v, 34)}   /eD = {nstr(v/(mpf('2.718281828459045235360287471352662497757')*D), 12)}  eD = {nstr(mpf('2.718281828459045235360287471352662497757')*D, 12)}")

L14 = Lt(14); L15 = Lt(15); L18 = Lt(18); L19 = Lt(19); L20 = Lt(20)
TGT = mpf('47.6270568122')
print("\n== B. margins against the 1/48 threshold ==")
print("L_t(14) =", nstr(L14, 30))
print("L_t(18) =", nstr(L18, 30))
print("L_t(19) =", nstr(L19, 30))
print("L_t(20) =", nstr(L20, 30))
print("L_t(14) exact  = 14^14/13^13 =", 14**14, "/", 13**13)
print("TGT/L_t(14) =", nstr(TGT/L14, 20))
print("TGT/L_t(18) =", nstr(TGT/L18, 20))
print("TGT/L_t(19) =", nstr(TGT/L19, 20))
print("L_t(20)/TGT =", nstr(L20/TGT, 20), " (= needed relative improvement -1)")
print("L_t(20)/L_t(14) =", nstr(L20/L14, 20))
print("(main target) 1/Phi(0,1/48) =", nstr(1/Phi0(mpf(1)/48), 20))
Dstar = findroot(lambda D: Lt(D) - TGT, mpf('18.05'))
print("D* : L_t(D*) = 47.6270568122 ->", nstr(Dstar, 20), " (slots to remove from 20:", nstr(20-Dstar, 12), ")")
print("D** : L_t(D**) = 48 ->", nstr(findroot(lambda D: Lt(D) - 48, mpf('18.2')), 20))

print("\n== C. window beta_max(lambda) = Phi(0,.)^{-1}(1/lambda) ==")
for name, lam in [("L_t(14)", L14), ("L_t(18)", L18), ("L_t(20)", L20), ("47.6270568122", TGT)]:
    b = bmax(lam)
    print(f"{name:>14s}: beta_max = {nstr(b, 20)}   /(1/48) = {nstr(b*48, 16)}")
print("Phi(0,1/48) =", nstr(Phi0(mpf(1)/48), 20), " 1/that =", nstr(1/Phi0(mpf(1)/48), 20))

print("\n== D. pointwise KP form (1/n)C(Dn,n-1) vs a_n  (exact integers) ==")
A = [0, 1, 20, 458, 11132, 280175, 7218768]
def kp(D, n):
    return comb(D*n, n-1)//n if comb(D*n, n-1) % n == 0 else comb(D*n, n-1)/n
for D in (14, 18, 20):
    print(f"--- D = {D} ---")
    for n in range(1, 7):
        v = kp(D, n)
        ok = "OK " if A[n] <= v else "FAIL"
        print(f" n={n}  a_n={A[n]:>8d}  KP_D(n)={str(v):>12}  {ok}  a_n/KP = {nstr(mpf(A[n])/mpf(v), 8)}")
print("\nminimal integer D with a_n <= (1/n)C(Dn,n-1) for each n:")
for n in range(1, 7):
    D = 1
    while kp(D, n) < A[n]:
        D += 1
    print(f" n={n}: D_min = {D}")

print("\n== E. monotonicity / rate in the degree argument ==")
print("D : L_t(D) : dL/dD (numeric) : L_t(D)/D")
prev = None
for D in range(12, 23):
    v = Lt(D)
    d = (Lt(D+mpf('0.001')) - Lt(D-mpf('0.001')))/mpf('0.002')
    print(f"{D:2d} : {nstr(v, 14):>18s} : {nstr(d, 10):>12s} : {nstr(v/D, 10)}")

print("\n== F. 2D acid-test analogue: Delta(d=2) = 4(2d-3) = 4 ==")
print("L_t(4) = 256/27 =", nstr(Lt(4), 20), " e*4 =", nstr(mpf('2.718281828459045235360287471352662497757')*4, 12))

print("\n== G. growth-sequence (boundary-slope) route is vacuous: prod_{i<n}(14i+6) ==")
from mpmath import gamma
for n in (2, 3, 5, 10, 20, 50):
    # prod_{i=1}^{n-1}(14i+6) = 14^{n-1} Gamma(n+3/7)/Gamma(1+3/7)
    val = mpf(14)**(n-1)*gamma(n+mpf(3)/7)/gamma(1+mpf(3)/7)
    print(f" n={n:3d}  prod^(1/n) = {nstr(val**(mpf(1)/n), 12)}   (a_n^(1/n) <= this; diverges ~ 14n/e)")
