import numpy as np
from scipy.special import iv
def d1(f,z,h=1e-6): return (f(z+h)-f(z-h))/(2*h)
r=lambda z: iv(1,z)/iv(2,z)
def D(z):
    P=iv(1,z)*iv(2,z); I1p=iv(1,z)+iv(2,z)/z
    a=2*z*I1p/iv(1,z); ap=d1(lambda t: 2*t*(iv(1,t)+iv(2,t)/t)/iv(1,t),z)
    Pp=d1(lambda t: iv(1,t)*iv(2,t),z)
    return P/z-(Pp*(1+a)-P*ap)/(1+a)**2
print("check  D/I2^2 = -r^2 (3 r^2 z - 12 r - 7 z)/(z (3r+2z)^2)   with r=I1/I2")
print("%8s %20s %20s %20s %18s"%("z","D (direct)","r","formula","3r^2z-12r-7z"))
for z in (0.1,0.5,1.0,2.0,5.0,10.0,20.0):
    rr=r(z); f=-rr**2*(3*rr**2*z-12*rr-7*z)/(z*(3*rr+2*z)**2)*iv(2,z)**2
    print("%8.2f %20.10e %20.10e %20.10e %18.6e"%(z,D(z),rr,f,3*rr**2*z-12*rr-7*z))
print()
print("the needed bound:   r < 2/z + sqrt(4/z^2 + 7/3)     vs the proved bound rho = 2/z + sqrt(1+4/z^2)")
rho=lambda z: 2/z+np.sqrt(1+4/z**2)
need=lambda z: 2/z+np.sqrt(4/z**2+7/3)
for z in (0.1,1.0,5.0,20.0):
    print("   z=%5.2f : r=%.6f  rho=%.6f  need=%.6f   r<rho: %s   rho<need: %s"%(
        z,r(z),rho(z),need(z), r(z)<rho(z), rho(z)<need(z)))
