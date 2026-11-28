import numpy as np
from scipy.special import iv
print("identity check:  I3 = I1 - (4/z) I2   (from  I_{nu-1}-I_{nu+1} = (2 nu/z) I_nu  with nu=2)")
print("%8s %20s %20s"%("z","I3 exact","I1-4I2/z"))
for z in (0.1,0.5,1.0,3.0,8.0,20.0):
    print("%8.2f %20.12e %20.12e"%(z,iv(3,z),iv(1,z)-4*iv(2,z)/z))
print()
print("hence  A = I1 I3 - I2^2  =  I1^2 - I2^2 - 4 I1 I2 / z ,   and  A<0  <=>  I1/I2 - I2/I1 < 4/z")
print("%8s %20s %20s %20s"%("z","A (direct)","I1^2-I2^2-4I1I2/z","I1/I2-I2/I1"))
for z in (0.1,0.5,1.0,3.0,8.0,20.0):
    A=iv(1,z)*iv(3,z)-iv(2,z)**2
    B=iv(1,z)**2-iv(2,z)**2-4*iv(1,z)*iv(2,z)/z
    C=iv(1,z)/iv(2,z)-iv(2,z)/iv(1,z)
    print("%8.2f %20.12e %20.12e %20.12e   (4/z = %.6f)"%(z,A,B,C,4/z))
