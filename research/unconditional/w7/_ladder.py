import numpy as np
from scipy.special import ive
from scipy.linalg import expm
bt=2.5
jm=60; nJ=jm+1
Ap=np.zeros((nJ,nJ)); idx=np.arange(nJ-1); Ap[idx,idx+1]=1.0; Ap[idx+1,idx]=1.0
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(nJ)])
print("t_j and -ln t_j (beta_t=2.5):")
for j in range(8):
    print("   j=%d  t_j=%.10f   -ln t_j=%.10f"%(j,td[j],-np.log(td[j]) if td[j]>0 else np.inf), flush=True)
print()
print("=== m a vs the -ln t_j ladder ===")
print()
for bs in (0.5,1.0,2.0,4.0,8.0):
    Mf=expm(bs*Ap/4)
    w0=Mf[:,0]; sc=float(w0 @ (td*w0))
    T=(Mf @ np.diag(td) @ Mf)/sc
    ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
    r=ev[1]/ev[0]; ma=-np.log(r)
    # find where m a falls in the ladder
    lad=[-np.log(td[j]) for j in range(1,10)]
    i=0
    while i<len(lad) and lad[i]<ma: i+=1
    print("beta_s=%.1f : m a=%.8f  falls between -ln t_%d=%.6f and -ln t_%d=%.6f"%(
        bs,ma,i-1,lad[i-2] if i>=2 else 0.0,i,lad[i-1] if i<len(lad) else np.inf), flush=True)
print()
print("=== 2-link model ===")
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
import functools
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm2=5
B=[(j,m,n) for j in range(0,jm2+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; dl=len(B)
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((dl,dl))
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
tot=None
for a in (-1,1):
    for b in (-1,1):
        A=sp.csr_matrix(mult_op(1,a,b)); Bm=sp.csr_matrix(mult_op(1,a,b,conj=True))
        term=sp.kron(A,Bm,format='csr')
        tot=term if tot is None else tot+term
W2=(0.5*tot).tocsr()
D2=dl**2
js=np.array([td[b[0]] for b in B])
Kdiag=np.ones(D2)
for ii in range(D2):
    x=ii; p=1.0
    for _ in range(2):
        p*=js[x%dl]; x//=dl
    Kdiag[ii]=p
Kd=sp.diags(Kdiag)
for bs in (0.5,1.0,2.0):
    def Mf(v): return np.exp(-2*bs)*expm_multiply(bs*W2, v)
    e0=np.zeros(D2); e0[0]=1.0
    w0=Mf(e0); sc=float(np.dot(w0,Kdiag*w0))
    def T(v): return Mf(Kd*(Mf(v)))/sc
    op=LinearOperator((D2,D2),matvec=T,dtype=float)
    ev=eigsh(op,k=2,which="LA",return_eigenvectors=False,tol=1e-13)
    ev=np.sort(ev)[::-1]
    r=ev[1]/ev[0]; ma=-np.log(r)
    print("beta_s=%.1f : m a=%.8f   (2-link)"%(bs,ma), flush=True)
