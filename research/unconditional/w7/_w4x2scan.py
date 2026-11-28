import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=1; bt=2.5
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B)
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
u={}
for M2 in (-1,1):
    for N2 in (-1,1):
        u[(M2,N2)]=sp.csr_matrix(mult_op(1,M2,N2))
        u[(M2,N2,'c')]=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
idx=[(-1,-1),(-1,1),(1,-1),(1,1)]
Ip=sp.identity(d,format='csr')
def Wilson4(L,a,b,c,dd):
    tot=None
    for (p,q) in idx:
        A=u[(p,q)]
        for (r,s) in idx:
            Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
            mats=[Ip]*L
            mats[a]=A; mats[b]=Bm; mats[c]=C; mats[dd]=Dm
            term=mats[0]
            for k in range(1,L): term=sp.kron(term,mats[k],format='csr')
            tot=term if tot is None else tot+term
    return (0.5*tot).tocsr()
L=7
W1=Wilson4(L,0,1,2,3); W2=Wilson4(L,0,4,5,6)
Wtot=(0.5*(W1+W2)).tocsr()
D=d**L
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
Kd=np.ones(D)
for ii in range(D):
    x=ii; prod=1.0
    for _ in range(L):
        prod*=td[B[x%d][0]]; x//=d
    Kd[ii]=prod
for bs in (0.5,1.0,2.0,4.0):
    def Mf_mult(v,bs=bs): return np.exp(-2*bs)*expm_multiply(bs*Wtot, v)
    e0=np.zeros(D); e0[0]=1.0
    w0=Mf_mult(e0); scale=float(np.dot(w0,Kd*w0))
    def T_mult(v,bs=bs): return Mf_mult(Kd*Mf_mult(v))/scale
    op=LinearOperator((D,D),matvec=T_mult,dtype=float)
    ev,V=eigsh(op,k=3,which="LA",tol=1e-11,maxiter=4000)
    o=np.argsort(ev)[::-1]; ev=ev[o]; V=V[:,o]
    print("2 plaquettes, beta_s=%.2f : r=%.10f  m a=%.6f"%((bs,ev[1]/ev[0],-np.log(ev[1]/ev[0]))), flush=True)
