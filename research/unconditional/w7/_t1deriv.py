from mpmath import mp, besseli, mpf
mp.dps = 60
print("=== T1 monotonicity: exact form of D'(z) ===")
print()
print("I_nu' = (I_{nu-1}+I_{nu+1})/2")
print("2 I_2 I_2' = I_2 (I_1+I_3)")
print("I_1' I_3 = (I_0+I_2) I_3 / 2")
print("I_1 I_3' = I_1 (I_2+I_4) / 2")
print()
print("D' = I_2(I_1+I_3) - (I_0+I_2)I_3/2 - I_1(I_2+I_4)/2")
print("   = (1/2)[ I_1 I_2 + I_2 I_3 - I_0 I_3 - I_1 I_4 ]")
print()
print("So monotonicity D'>0  <=>  I_1 I_2 + I_2 I_3 > I_0 I_3 + I_1 I_4")
print()
ok=True
for z in ("1e-6","1e-4","1e-2","0.1","0.5","1","2","5","10","29","100","1000","10000"):
    zz=mpf(z)
    I0=besseli(0,zz); I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz); I4=besseli(4,zz)
    lhs=I1*I2+I2*I3
    rhs=I0*I3+I1*I4
    diff=rhs-lhs
    if diff>=0: ok=False
    print("  z=%-8s  rhs-lhs = %s   (<0 ? %s)"%(z, mp.nstr(diff,12), diff<0))
print()
print("all negative? %s"%ok)
print()
print("=== asymptotics of (rhs-lhs) ===")
print("small z: I_0~1, I_1~z/2, I_2~z^2/8, I_3~z^3/48, I_4~z^4/384")
print("  lhs ~ (z/2)(z^2/8) + (z^2/8)(z^3/48) = z^3/16 + z^5/384")
print("  rhs ~ 1*(z^3/48) + (z/2)(z^4/384) = z^3/48 + z^5/768")
print("  rhs-lhs ~ z^3(1/48-1/16) + z^5(1/768-1/384) = -z^3/24 - z^5/768")
print("  => rhs-lhs ~ -z^3/24  (negative, matches D' ~ z^3/48 = 4 D/z)")
for z in ("1e-4","1e-3","1e-2"):
    zz=mpf(z)
    I0=besseli(0,zz); I1=besseli(1,zz); I2=besseli(2,zz); I3=besseli(3,zz); I4=besseli(4,zz)
    diff=(I0*I3+I1*I4)-(I1*I2+I2*I3)
    print("  z=%-8s  (rhs-lhs)/z^3 = %s   (-1/24 = %s)"%(z, mp.nstr(diff/(zz**3),15), mp.nstr(-mpf(1)/24,15)))
