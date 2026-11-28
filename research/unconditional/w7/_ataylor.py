import sympy as sp
z=sp.symbols("z")
N=12
def serI(nu):
    return sp.series(sp.besseli(nu,z),z,0,2*N+6).removeO()
A=sp.expand(serI(1)*serI(3)-serI(2)**2)
A=sp.series(A,z,0,2*N+6).removeO()
print("Taylor coefficients of A(z) = I1 I3 - I2^2:")
cs=[]
for n in range(0,N+1):
    c=sp.simplify(A.coeff(z,2*n))
    cs.append(c)
    if c!=0:
        print("   a_%-2d (coeff of z^%-2d) = %-24s = %+.10e"%(n,2*n,str(c),float(c)))
print()
print("check the recurrence for the nonzero ones:")
for n in range(3,N+1):
    if cs[n-1]!=0:
        r=sp.nsimplify(cs[n]/cs[n-1])
        pred=sp.Rational(2*n-1,1)/(2*(n+1)*(n-2)*(n+2))
        ok="OK" if sp.simplify(r-pred)==0 else "MISMATCH"
        print("   a_%d/a_%d = %-18s   pred = %-18s  %s"%(n,n-1,str(r),str(pred),ok))
