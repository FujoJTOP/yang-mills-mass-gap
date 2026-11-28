import numpy as np, functools
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def make_basis(jm):
    return [(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
def mult_op(B,index,J2,M2,N2,conj=False):
    d=len(B); A=np.zeros((d,d))
    mm,nn=(-M2,-N2) if conj else (M2,N2)
    ph=((-1.0)**((M2-N2)//2)) if conj else 1.0
    for i2,(j2,m2,n2) in enumerate(B):
        m3,n3=m2-mm,n2-nn
        for j3 in range(abs(J2-j2),J2+j2+1,2):
            key=(j3,m3,n3)
            if key not in index: continue
            v=((-1.0)**((m2-n2)//2))*d3(j3,m3,J2,mm,j2,-m2)*d3(j3,n3,J2,nn,j2,-n2)
            if v!=0.0: A[index[key],i2]=np.sqrt((j3+1.0)*(j2+1.0))*v
    return ph*A
jm=3
B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
# single-link Q operator (the chi_{1/2} / q operator):
Qs=np.zeros((d,d))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        Qs+=mult_op(B,index,1,M2,N2,False)
Qs=(Qs+Qs.T)/2
i0=0; ig=[index[(1,m,n)] for (m,n) in [(1,1),(1,-1),(-1,1),(-1,-1)]]
print("single-link Q on the j=0 and j=1/2 states:")
print("   <0|Q|0> = %+.10e"%Qs[i0,i0])
print("   4x4 block <g|Q|g'> :")
print(Qs[np.ix_(ig,ig)])
print("   cross block <0|Q|g> :", " ".join("%+.6f"%Qs[i0,g] for g in ig))
print()
print("=> if <0|Q|0>=0 and the 4x4 block is 0, then Q maps j=0 <-> j=1/2 (8x8 block = (0 B; B^dag 0))")
print("   B = <g|Q|0> outer product?  sigma(B):")
bvec=np.array([Qs[i0,g] for g in ig])
print("   b =",bvec,"   |b| = %.10f"%np.linalg.norm(bvec))
print("   => the 2-link Q block has sigma = |b| = %.6f  (so lambda_Q = +-%.6f)"%(np.linalg.norm(bvec),np.linalg.norm(bvec)))
