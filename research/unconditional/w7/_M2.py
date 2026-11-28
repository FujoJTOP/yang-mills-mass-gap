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
Ops={}
for J2 in range(0,jm+1):
    O=np.zeros((D,D))
    for M2 in range(-J2,J2+1,2):
        for N2 in range(-J2,J2+1,2):
            O+=np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
    Ops[J2]=(O+O.T)/2
def Mk(bs,scaled):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else:
            cJ=2*(J2+1)*(iv(J2+1,bs)*np.exp(-bs) if scaled else iv(J2+1,bs))/bs
        M+=cJ*Ops[J2]
    return M
for scaled in (False,True):
    h=1e-4
    M0=Mk(0.0,scaled); Mp=Mk(h,scaled); Mm=Mk(-h,scaled)
    M1=(Mp-Mm)/(2*h)
    M2=(Mp-2*M0+Mm)/h**2
    print("scaled=%s :  c1'(0) check  <e1,M1 e1> = %+.10e ; <e1,M2 e1> = %+.10e"%(scaled,M0[0]@(M1@M0[0].T) if False else (np.eye(D)[0]@(M1@np.eye(D)[0])), np.eye(D)[0]@(M2@np.eye(D)[0])))
    print("             <e1,M1^2 e1> = %+.10e"%(np.eye(D)[0]@(M1@(M1@np.eye(D)[0]))))
