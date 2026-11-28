import numpy as np, functools
from scipy.special import ive
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
def Mf_of(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*ive(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    return M
# n_J on the tensor basis:  n_j = (j+1)^2 - 1   (j is doubled)
nvec=np.array([((B[i1][0]+1)**2-1)+((B[i2][0]+1)**2-1) for i1 in range(d) for i2 in range(d)],float)
print("analytic prediction:  Phi = (1/2) [ <psi2| n_J |psi2> - <psi1| n_J |psi1> ]")
print("where psi_1, psi_2 are the top two eigenvectors of M_f (biggest eigenvalues)")
print()
print("%8s %20s %20s %20s"%("beta_s","<n>_1","<n>_2","Phi_pred"))
for bs in (0.5,2.0,4.0,8.0,16.0,32.0):
    M=Mf_of(bs)
    A=(M+M.T)/2
    ev,V=np.linalg.eigh(A); idx=np.argsort(ev)[::-1]
    p1=V[:,idx[0]]; p2=V[:,idx[1]]
    n1=p1@(nvec*p1); n2=p2@(nvec*p2)
    print("%8.1f %20.8f %20.8f %20.8f"%(bs,n1,n2,0.5*(n2-n1)))
