from mpmath import mp, besseli, mpf
import numpy as np
mp.dps = 80
print("=== T1: monotonicity of D(z) and the log-derivative ===")
print()
# check monotonicity of D
zs = [mpf("1e-4") * (mpf(10) ** (mpf(k)/mpf(20))) for k in range(0, 101)]
zs = sorted(set(zs))
Ds = [besseli(2,z)**2 - besseli(1,z)*besseli(3,z) for z in zs]
mono = all(Ds[i] < Ds[i+1] for i in range(len(Ds)-1))
print("D strictly increasing over %d log-spaced points? %s" % (len(zs), mono))
print()
# log-derivative
print("log-derivative d/dz ln D(z) at selected z:")
for z in ("1e-3","1e-2","0.1","0.5","1","2","5","10","29","100","1000"):
    zz=mpf(z)
    I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz)
    I1p=(besseli(0,zz)+besseli(2,zz))/2
    I2p=(besseli(1,zz)+besseli(3,zz))/2
    I3p=(besseli(2,zz)+besseli(4,zz))/2
    D=I2*I2-I1*I3
    Dp=2*I2*I2p-I1p*I3-I1*I3p
    print("   z=%-8s  d ln D/dz = %s"%(z, mp.nstr(Dp/D,20)))
print()
print("asymptotics of d ln D/dz:  z->0: 4/z ;  z->inf: 1/z")
print()
for z in ("1e-3","1e-2","1e-1"):
    zz=mpf(z)
    I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz)
    I1p=(besseli(0,zz)+besseli(2,zz))/2
    I2p=(besseli(1,zz)+besseli(3,zz))/2
    I3p=(besseli(2,zz)+besseli(4,zz))/2
    D=I2*I2-I1*I3
    Dp=2*I2*I2p-I1p*I3-I1*I3p
    print("   z=%-8s  z*d ln D/dz = %s   (should -> 4)"%(z, mp.nstr(zz*Dp/D,15)))
for z in ("100","1000","10000"):
    zz=mpf(z)
    I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz)
    I1p=(besseli(0,zz)+besseli(2,zz))/2
    I2p=(besseli(1,zz)+besseli(3,zz))/2
    I3p=(besseli(2,zz)+besseli(4,zz))/2
    D=I2*I2-I1*I3
    Dp=2*I2*I2p-I1p*I3-I1*I3p
    print("   z=%-8s  z*d ln D/dz = %s   (should -> 1)"%(z, mp.nstr(zz*Dp/D,15)))
