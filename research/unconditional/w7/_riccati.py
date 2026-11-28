import numpy as np
from scipy.special import iv
def d1(f,z,h=1e-6): return (f(z+h)-f(z-h))/(2*h)
r=lambda z: iv(1,z)/iv(2,z)
print("Riccati check:   r' = 1 - r^2 + 3 r / z      with r = I1/I2")
print("%8s %20s %20s"%("z","r'(FD)","1-r^2+3r/z"))
for z in (0.1,0.5,1.0,2.0,5.0,10.0,20.0):
    a=d1(r,z); b=1-r(z)**2+3*r(z)/z
    print("%8.2f %20.12e %20.12e"%(z,a,b))
print()
print("barrier rho(z) = 2/z + sqrt(1+4/z^2):   compare rho' with F(rho)=1-rho^2+3rho/z")
for z in (0.1,0.5,1.0,3.0,8.0,20.0):
    rho=lambda t: 2/t+np.sqrt(1+4/t**2)
    a=d1(rho,z); b=1-rho(z)**2+3*rho(z)/z
    print("%8.2f : rho\x27=%+.8e   F(rho)=%+.8e   rho\x27-F=%+.6e   r=%.6f  rho=%.6f"%(z,a,b,a-b,r(z),rho(z)))
print()
print("small-z asymptotics:  r ~ 4/z + z/6 ,  rho ~ 4/z + z/4")
for z in (0.01,0.05,0.1):
    rho=lambda t: 2/t+np.sqrt(1+4/t**2)
    print("   z=%.2f : r-4/z = %+.8e (z/6=%.8e)   rho-4/z = %+.8e (z/4=%.8e)"%(
        z, r(z)-4/z, z/6, rho(z)-4/z, z/4))
