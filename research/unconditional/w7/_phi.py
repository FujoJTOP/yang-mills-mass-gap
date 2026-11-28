import mpmath as mp
mp.mp.dps=60
from fractions import Fraction as F
ck=[F(3,2),F(3,4),F(3,16),F(-9,16),F(-2331,1280),F(-279,64),F(-153063,14336),F(-3861,128),
    F(-6708519,65536),F(-2108187,5120)]
def exact(b): return -mp.log(mp.besseli(2,b)/mp.besseli(1,b))
def ser(b,K): return sum(mp.mpf(ck[k-1].numerator)/ck[k-1].denominator/b**k for k in range(1,K+1))
print("Is the remainder  Phi_K = -ln(I2/I1) - sum_{k<=K} c_k/beta^k  just the series tail?")
print("(if  Phi_K * beta^{K+1} -> c_{K+1}, then there is NO non-perturbative piece)")
print()
K=10
print("%10s %22s %22s"%("beta","Phi_K * beta^{K+1}","c_{K+1} = c_11"))
c11=-mp.mpf(10986758433)/5767168
for b in (50,100,200,400,800):
    rem=exact(mp.mpf(b))-ser(mp.mpf(b),K)
    print("%10d %22.12e %22.12e"%(b, rem*mp.mpf(b)**(K+1), c11))
print()
print("c_11 = %.12e"%c11)
