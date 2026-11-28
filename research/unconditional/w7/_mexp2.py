import numpy as np, functools
from scipy.special import ive
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=3
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((d,d))
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
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op(J2,M2,N2),mult_op(J2,M2,N2,conj=True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
i0=0; i11=index[(1,1,1)]
idx0=i0*d+i0; idx2=i11*d+i11
print("   beta_s     |M_{0,(1/2,1/2)}|      -ln|M|      |M|/e^{-b/2}     |M|/e^{-b/4}")
for bs in (1.0,2.0,3.0,4.0,6.0,8.0,12.0,16.0,24.0):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    mv=abs(M[idx0,idx2])
    print("%10.2f %20.10e %14.8f %16.8f %16.8f"%(bs,mv,-np.log(mv),mv/np.exp(-bs/2),mv/np.exp(-bs/4)), flush=True)

