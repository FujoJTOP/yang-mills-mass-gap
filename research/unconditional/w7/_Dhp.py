import mpmath as mp
mp.mp.dps=60
def I(n,z): return mp.besseli(n,z)
def Dval(z):
    I1=I(1,z); I2=I(2,z); I1p=(I(0,z)+I(2,z))/2
    I1pp=I1+I1/z**2-I1p/z
    P=I1*I2
    Pp=I1p*I2+I1*((I1+I(3,z))/2)
    a=2*z*I1p/I1
    ap=2*((I1p+z*I1pp)*I1-z*I1p*I1p)/I1**2
    return P/z-(Pp*(1+a)-P*ap)/(1+a)**2
print("60-digit check of  D  and of the symbolic formula  D/I2^2 = -r^2(3r^2 z-12r-7z)/(z(3r+2z)^2)")
print("%8s %26s %26s %26s"%("z","D (60 digits)","formula","ratio"))
for zz in ("0.1","0.5","1","2","5","10","20"):
    z=mp.mpf(zz)
    D=Dval(z)
    r=I(1,z)/I(2,z); I2=I(2,z)
    F=-r**2*(3*r**2*z-12*r-7*z)/(z*(3*r+2*z)**2)*I2**2
    print("%8s %26.16e %26.16e %26.12f"%(zz,D,F,(D/F) if F!=0 else mp.nan))
