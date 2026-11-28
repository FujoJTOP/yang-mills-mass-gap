import numpy as np, functools, time
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
def Phi_of(jm,bs,Jmax=None):
    B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
    Jm=jm if Jmax is None else min(Jmax,jm)
    M=np.zeros((D,D))
    for J2 in range(0,Jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
        for M2 in range(-J2,J2+1,2):
            for N2 in range(-J2,J2+1,2):
                M+=cJ*np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
    M=(M+M.T)/2
    ev=eigh(M,subset_by_index=[D-2,D-1],eigvals_only=True)
    ev=np.sort(ev)[::-1]
    return D,-2*np.log(ev[1]/ev[0])
print("does Phi(beta_s=16) go to 0 as j_max grows?  (continuum => mu_2 = mu_1)")
print("%10s %10s %10s %20s"%("j_max","d","D","Phi(16)"))
for jm in (3,4,5,6):
    t0=time.time()
    D,phi=Phi_of(jm,16.0)
    print("%10.1f %10d %10d %20.10f   [%.0fs]"%(jm/2,D,D,phi,time.time()-t0), flush=True)
