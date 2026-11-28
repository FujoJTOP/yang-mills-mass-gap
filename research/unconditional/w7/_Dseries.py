import sympy as sp
z=sp.symbols("z",positive=True)
ORD=30
I1=sp.series(sp.besseli(1,z),z,0,ORD).removeO()
I2=sp.series(sp.besseli(2,z),z,0,ORD).removeO()
I1p=sp.diff(I1,z)
P=sp.expand(I1*I2)
a=sp.series(2*z*I1p/I1,z,0,ORD-4).removeO()
LHS=sp.series(P/z,z,0,ORD-4).removeO()
RHS=sp.series(sp.diff(P/(1+a),z),z,0,ORD-8).removeO()
D=sp.series(LHS-RHS,z,0,ORD-8).removeO()
print("Taylor expansion of  D(z) = P/z - d/dz[P/(1+a)] :")
nz=0; neg=0
for k in range(0,ORD-10):
    c=sp.simplify(D.coeff(z,k))
    if c!=0:
        nz+=1
        sgn="+" if c>0 else "-"
        if c<0: neg+=1
        print("   z^%-2d : %s %-24s"%(k,sgn,str(c)))
print()
print("nonzero entries: %d,  negative entries: %d"%(nz,neg))
