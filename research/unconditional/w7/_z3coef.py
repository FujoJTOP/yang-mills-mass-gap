import sympy as sp
z=sp.symbols("z",positive=True)
ORD=14
I1=sp.series(sp.besseli(1,z),z,0,ORD).removeO()
I2=sp.series(sp.besseli(2,z),z,0,ORD).removeO()
I3=sp.series(sp.besseli(3,z),z,0,ORD).removeO()
I1p=sp.diff(I1,z)
P=sp.expand(I1*I2)
a=sp.series(2*z*I1p/I1,z,0,ORD-2).removeO()
LHS=sp.series(P/z,z,0,ORD-2).removeO()
RHS=sp.series(sp.diff(P/(1+a),z),z,0,ORD-4).removeO()
D=sp.series(LHS-RHS,z,0,ORD-4).removeO()
A=sp.series(I1*I3-I2**2,z,0,ORD).removeO()
print("expansion of  P/z - d/dz[P/(1+a)]   and of   I1 I3 - I2^2 :")
print()
for k in range(0,8):
    cd=sp.simplify(D.coeff(z,k)); ca=sp.simplify(A.coeff(z,k))
    if cd!=0 or ca!=0:
        print("  z^%-2d :  D = %-22s   A = %-22s"%(k,str(cd),str(ca)))
print()
print("ratio D/A term by term (for the nonzero A terms):")
for k in range(2,8):
    ca=sp.simplify(A.coeff(z,k))
    cd=sp.simplify(D.coeff(z,k))
    if ca!=0:
        print("  z^%-2d : D/A = %s"%(k,str(sp.nsimplify(cd/ca))))
