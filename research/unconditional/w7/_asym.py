import numpy as np
from scipy.special import iv
# asymptotic expansion of  -ln(I2/I1)  in powers of 1/beta
def f(b): return -np.log(iv(2,b)/iv(1,b))
bs=np.array([40.0,60.0,80.0,120.0,160.0,240.0,320.0])
vals=np.array([f(b) for b in bs])
# fit  f = c1/b + c2/b^2 + c3/b^3 + c4/b^4 + c5/b^5
A=np.vstack([1.0/bs**k for k in range(1,6)]).T
coef,res,rank,sv=np.linalg.lstsq(A,vals,rcond=None)
print("asymptotic expansion of  -ln(I2/I1)  in 1/beta:")
for k,c in enumerate(coef,1):
    print("   c_%d = %+.10e"%(k,c))
print()
print("check: f(beta) vs the 5-term fit")
for b in (40.0,80.0,160.0,320.0,640.0,1280.0):
    fit=sum(coef[k-1]/b**k for k in range(1,6))
    print("   beta=%7.1f : f=%.10e   fit=%.10e   rel.diff=%.2e"%(b,f(b),fit,abs(f(b)-fit)/f(b)))
print()
print("=> c_1 = 3/2 ?   c_2 = 15/8 ?  (the free-field-like coefficient)")
