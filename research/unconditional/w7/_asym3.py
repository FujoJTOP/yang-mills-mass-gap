import numpy as np
from scipy.special import iv
from fractions import Fraction as F
c=[F(3,2),F(3,4),F(3,16),F(-9,16),F(-2331,1280),F(-279,64),F(-153063,14336),F(-3861,128)]
def f(b): return -np.log(iv(2,b)/iv(1,b))
def approx(b,K): return sum(float(c[k-1])/b**k for k in range(1,K+1))
print("verify the exact symbolic expansion against iv() numerically")
print("%8s %20s %20s %12s"%("beta","-ln(I2/I1)","8-term","rel.diff"))
for b in (40.0,80.0,160.0,320.0,640.0,1280.0,2560.0):
    try:
        v=f(b)
    except Exception:
        v=float("nan")
    a=approx(b,8)
    print("%8.1f %20.12e %20.12e %12.2e"%(b,v,a,abs(v-a)/abs(v)))
print()
print("coefficients (exact):")
for k in range(1,9):
    print("   c_%d = %-20s = %+.10f   |c_k|/|c_{k-1}| = %s"%(k,c[k-1],float(c[k-1]),
          "%.4f"%(abs(float(c[k-1]))/abs(float(c[k-2]))) if k>1 else "-"))
