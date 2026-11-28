import numpy as np, functools
from scipy.linalg import eigh
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
Q12=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        Q12+=np.kron(mult_op(B,index,1,M2,N2,False),mult_op(B,index,1,M2,N2,True))
Q12=(Q12+Q12.T)/2
i0=0; ig=[index[(1,m,n)] for (m,n) in [(1,1),(1,-1),(-1,1),(-1,-1)]]
# G1 = link1 at j=1/2, link2 at 0 ;  G2 = link1 at 0, link2 at j=1/2
G1=[i*d+i0 for i in ig]
G2=[i0*d+i for i in ig]
M8=np.zeros((8,8))
idx8=G1+G2
for a,i in enumerate(idx8):
    for b,j in enumerate(idx8):
        M8[a,b]=Q12[i,j]
print("8x8 block of Q_12 on G1+G2 :")
print(np.round(M8,6))
Bh=M8[:4,4:]
print()
print("the 4x4 off-diagonal block B :")
print(np.round(Bh,6))
print("   singular values of B =", np.round(np.linalg.svd(Bh,compute_uv=False),8))
ev=np.linalg.eigvalsh(M8)
print("   eigenvalues of the 8x8 block =", np.round(np.sort(ev),8))
print()
print("=> if sigma(B)=0.5 then the 12x12 cyclic block (three such B's) has lambda_max = 2*0.5 = 1")
