import numpy as np, functools
from scipy.special import ive
from scipy.linalg import eigh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=4
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
bt=2.5
td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
tvec=np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B])
hvec=np.array([b1[0]+b2[0] for b1 in B for b2 in B])
def r_trunc(bs,hmax):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    keep=np.where(hvec<=hmax)[0]
    if len(keep)<2: return float('nan')
    Mk=M[np.ix_(keep,keep)]
    Kk=np.diag(tvec[keep])
    T=Mk@Kk@Mk; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[1]/ev[0]
print("L_s=1 , jm=%d , D=%d : r from the h<=hmax truncation vs the full answer"%(jm,D))
print("%8s %14s %14s %14s %14s"%("beta_s","h<=2","h<=4","h<=6","full(h<=12)"))
for bs in (0.25,0.5,1.0,2.0,4.0):
    row=[]
    for hm in (2,4,6,8):
        row.append(r_trunc(bs,hm))
    print("%8.2f %14.10f %14.10f %14.10f %14.10f"%(bs,*row), flush=True)

