import numpy as np
from scipy.special import iv
def d1(f,z,h=1e-6): return (f(z+h)-f(z-h))/(2*h)
def R(z):
    P=iv(1,z)*iv(2,z); I1p=d1(lambda t: iv(1,t),z)
    return P/(1+2*z*I1p/iv(1,z))
print("derivative inequality:   I1 I2 / z   >   d/dz [ I1 I2 / (1 + 2 z I1'/I1) ]")
print("%8s %20s %20s %20s"%("z","LHS=P/z","RHS'","LHS-RHS'"))
for z in (0.1,0.3,0.5,1.0,2.0,3.0,5.0,8.0,12.0):
    P=iv(1,z)*iv(2,z)
    L=P/z
    Rp=d1(R,z)
    print("%8.2f %20.10e %20.10e %20.10e"%(z,L,Rp,L-Rp))
print()
print("also test the algebraically equivalent form   (1+a)[P(1+a)-zP'] + zP a'  > 0   with a=2zI1'/I1")
for z in (0.1,0.5,1.0,3.0,8.0):
    I1=iv(1,z); I1p=d1(lambda t: iv(1,t),z); I1pp=None
    a=2*z*I1p/I1
    ap=2*(I1p+I1p*z)/I1-2*z*I1p*I1p/I1**2 if False else 2*((I1p+z*d1(lambda t: iv(1,t),z)) )  # placeholder
    print("   z=%5.2f : a = %+.6f   (a-2 = %+.3e)"%(z,a,a-2))
