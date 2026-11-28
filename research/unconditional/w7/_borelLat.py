import numpy as np, functools
from scipy.special import ive
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
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
def Mf(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    return (M+M.T)/2
def lnr(bs,bt):
    td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
    Kd=np.diag(np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B]))
    M=Mf(bs)
    T=M@Kd@M; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return np.log(ev[1]/ev[0])
bs=2.0
bts=np.array([32.,64.,128.,256.,512.,1024.])
fs=np.array([lnr(bs,b) for b in bts])
print("ln r(beta_s=2, beta_t) : does it converge to a positive constant (Phi>0) or decay?")
for b,f in zip(bts,fs): print("   beta_t=%7.0f : ln r = %.12f"%(b,f))
print()
print("fit  ln r = c0 + c1/b + c2/b^2 + c3/b^3 :")
A=np.vstack([np.ones_like(bts),1/bts,1/bts**2,1/bts**3]).T
coef,_,_,_=np.linalg.lstsq(A,fs,rcond=None)
for k,c in enumerate(coef): print("   c_%d = %+.10f"%(k,c))
print()
print("=> c0>0 means ln r does NOT decay to 0 => r -> exp(-c0) < 1 => a residual gap (constant in lattice units)")
