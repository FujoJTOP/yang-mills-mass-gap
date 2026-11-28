from mpmath import mp, besseli, mpf
mp.dps = 80
print("=== T1: exact asymptotics of D(z) = I_2^2 - I_1 I_3 ===")
print()
print("SMALL z: I_1 ~ z/2, I_2 ~ z^2/8, I_3 ~ z^3/48")
print("  D = I_2^2 - I_1 I_3 ~ z^4/64 - z^4/96 = z^4/192")
print()
for z in ("1e-6","1e-4","1e-2"):
    zz=mpf(z)
    I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz)
    D=I2*I2-I1*I3
    print("  z=%-8s  D/z^4 = %s   (1/192 = %s)"%(z, mp.nstr(D/(zz**4),20), mp.nstr(mpf(1)/192,20)))
print()
print("LARGE z: I_nu(z) ~ e^z/sqrt(2 pi z) [1 - (4 nu^2-1)/(8z) + ...]")
print("  D = I_2^2 - I_1 I_3  ->  I_2^2 / z")
print()
for z in (100,1000,10000,100000):
    zz=mpf(z)
    I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz)
    D=I2*I2-I1*I3
    print("  z=%-8d  z*D/I_2^2 = %s"%(z, mp.nstr(zz*D/(I2*I2),20)))
print()
print("==> both asymptotics are EXACT:  D ~ z^4/192 (z->0),  D ~ I_2^2/z (z->inf)")
