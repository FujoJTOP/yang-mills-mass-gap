import numpy as np
from scipy.special import iv
def d2(f,z,h=1e-4): return (f(z+h)-2*f(z)+f(z-h))/h**2
f=lambda z: iv(2,z)/iv(1,z)
print("look for the identity:  (I2/I1)''  =  - (something) * (I1 I3 - I2^2)")
print("%8s %18s %18s %18s %18s %18s"%("z","(I2/I1)''","I1I3-I2^2","ratio","ratio*z^2*I1^2","ratio*z^2"))
for z in (0.1,0.3,0.5,1.0,2.0,3.0,5.0,8.0):
    q=d2(f,z); A=iv(1,z)*iv(3,z)-iv(2,z)**2
    r=q/A
    print("%8.2f %18.8e %18.8e %18.8e %18.8e %18.8e"%(z,q,A,r,r*z*z*iv(1,z)**2,r*z*z))
