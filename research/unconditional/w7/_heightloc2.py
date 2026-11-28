import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
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
Q=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        Q+=np.kron(mult_op(1,M2,N2),mult_op(1,M2,N2,conj=True))
Q=(Q+Q.T)/2
bt=2.5
td=ive(np.arange(0,jm+1)+1,bt)/ive(1,bt)
tvec=np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B])
hvec=np.array([b1[0]+b2[0] for b1 in B for b2 in B])
Qsp=sp.csr_matrix(Q)
def T_mult(v,bs):
    w=expm_multiply((bs/2)*Qsp, v)
    w=np.exp(-bs)*w
    w=tvec*w
    w=expm_multiply((bs/2)*Qsp, w)
    return np.exp(-bs)*w
bs=2.0
op=LinearOperator((D,D),matvec=lambda v: T_mult(v,bs),dtype=float)
ev,V=eigsh(op,k=1,which='LA',tol=1e-12)
psi=V[:,0]
print("L_s=1 two-link model, beta_s=%.1f ,  jm=%d (h_max=%d)"%(bs,jm,2*jm))
print("top eigenvector's weight by height h:")
hs=np.arange(0,2*jm+1,2)
wts=[np.sum(psi[hvec==h]**2)/np.sum(psi**2) for h in hs]
for h,w in zip(hs,wts):
    if w>1e-14:
        print("   h=%2d : %.6e"%(h,w))
print()
print("geometric ratio between consecutive populated rungs:")
prev=None
for h,w in zip(hs,wts):
    if w>1e-12:
        if prev is not None:
            print("   h=%2d : ratio = %.6f"%(h,w/prev))
        prev=w

