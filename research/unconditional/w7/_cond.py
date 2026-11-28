import numpy as np
from scipy.special import iv
def d2(f,z,h=1e-4): return (f(z+h)-2*f(z)+f(z-h))/h**2
f=lambda z: np.log(iv(2,z)/iv(1,z))
print("condition  (ln)'' + (ln)'^2 < 0  reduces to")
print("   int_0^z I1 I2 / t  dt   >   I1 I2 / (1 + 2 z I1'/I1)")
print()
print("%8s %18s %18s %18s"%("z","LHS","RHS","LHS-RHS"))
for z in (0.1,0.3,0.5,1.0,2.0,3.0,5.0,8.0,12.0):
    from scipy.integrate import quad
    L=quad(lambda t: iv(1,t)*iv(2,t)/t,0,z,limit=200)[0]
    P=iv(1,z)*iv(2,z)
    d1=(iv(1,z+1e-5)-iv(1,z-1e-5))/2e-5
    R=P/(1+2*z*d1/iv(1,z))
    print("%8.2f %18.10e %18.10e %18.10e"%(z,L,R,L-R))
print()
print("also the direct check of (ln)''+(ln)'^2 :")
for z in (0.1,0.5,1.0,2.0,3.0,5.0,8.0):
    v=d2(f,z)
    d1=(f(z+1e-5)-f(z-1e-5))/2e-5
    print("   z=%5.2f : (ln)''+(ln)\x27^2 = %+.8e"%(z,v+d1*d1))
