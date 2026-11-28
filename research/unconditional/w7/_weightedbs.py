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
ops2={}
for J2 in range(0,jm+1):
    ops2[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
              for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
def Ops1(J2):
    A=np.zeros((d,d))
    for M2 in range(-J2,J2+1,2):
        for N2 in range(-J2,J2+1,2):
            A+=mult_op(B,index,J2,M2,N2,False)
    return (A+A.T)/2
O1={J2:Ops1(J2) for J2 in range(0,jm+1)}
bt=2.5
t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
tvec=np.array([t1d[b[0]] for b in B])
Kdiag=np.diag(tvec)
def Mf1(bp):
    A=np.zeros((d,d))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bp==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bp)/bp
        A+=cJ*O1[J2]
    return (A+A.T)/2
def model(bp):
    M=Mf1(bp); T=M@Kdiag@M; return (T+T.T)/2
def resid(Tp):
    nT=np.linalg.norm(Tp); best=(None,1e18)
    for bp in np.linspace(0.0,6.0,241):
        nE=np.linalg.norm(Tp-model(bp))
        if nE<best[1]: best=(bp,nE)
    return best[0],best[1]/nT
print("weighted block-spin variants: which one has a SMALL residual?")
for bs in (1.0,2.0):
    M2=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops2[J2]: M2+=cJ*op
    M2=(M2+M2.T)/2
    Kd=np.diag(np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B]))
    T=M2@Kd@M2; T=(T+T.T)/2
    Kn=np.diag(np.array([1.0/(t1d[b1[0]]*t1d[b2[0]]) for b1 in B for b2 in B]))
    Ks=np.diag(np.sqrt(np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])))
    Kms=np.diag(1.0/np.sqrt(np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])))
    variants={
      "Tr2(T)":T,
      "Tr2(K^1/2 T K^1/2)":Ks@T@Ks,
      "Tr2(K^-1/2 T K^-1/2)":Kms@T@Kms,
      "Tr2(M T M)":M2@T@M2,
      "Tr2(K^-1/2 M T M K^-1/2)":None,
    }
    variants["Tr2(K^-1/2 M T M K^-1/2)"]=Kms@(M2@T@M2)@Kms
    print()
    print("  beta_s = %.1f"%bs)
    for name,W in variants.items():
        Tp=np.einsum("ijik->jk",W.reshape(d,d,d,d)); Tp=(Tp+Tp.T)/2
        bp,r=resid(Tp)
        print("     %-28s : b'*=%.3f   ||E||/||T'|| = %.6f"%(name,bp,r), flush=True)
