import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=2; bs=1.0; bt=2.5
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
def Wilson4(a,b,c,dd):
    tot=None
    for (p,q) in idx:
        A=u[(p,q)]
        for (r,s) in idx:
            Bm=u[(q,r,'c')]; C=u[(r,s)]; Dm=u[(s,p,'c')]
            term=sp.kron(sp.kron(sp.kron(A,Bm,format='csr'),C,format='csr'),Dm,format='csr')
            tot=term if tot is None else tot+term
    return (0.5*tot).tocsr()
# two plaquettes: P1 on links (0,1,2,3), P2 on links (0,4,5,6)
W4a=Wilson4(0,1,2,3)  # acts on the first 4 of the 7 links
W4b=Wilson4(0,4,5,6)  # acts on links 0,4,5,6
D=d**7
print("jm=%d d=%d D=%d  nnz(W4a)=%d nnz(W4b)=%d"%(jm,d,D,W4a.nnz,W4b.nnz), flush=True)
# index mapping for W4b: we need W4b acting on (l0,l4,l5,l6) = local indices (0,1,2,3) of W4b
# reshape v (d^7) -> tensor of shape (d,d,d,d,d,d,d); W4b acts on axes (0,4,5,6)
# use np.tensordot on the 4 axes
shape7=(d,)*7
def W4a_mv(v):
    T=v.reshape(d**4, d**3)
    return (W4a @ T).ravel()
def W4b_mv(v):
    T=v.reshape(shape7)
    # permute so axes (0,4,5,6) come first, then apply W4b, then permute back
    perm=(0,4,5,6,1,2,3)
    Tp=np.transpose(T,perm).reshape(d**4, d**3)
    out=(W4b @ Tp).reshape((d,)*7)
    # inverse permute
    inv=np.argsort(perm)
    return np.transpose(out,inv).ravel()
def Wtot_mv(v): return 0.5*(W4a_mv(v)+W4b_mv(v))
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
Kd=np.ones(D)
for ii in range(D):
    x=ii; prod=1.0
    for _ in range(7):
        prod*=td[B[x%d][0]]; x//=d
    Kd[ii]=prod
print("Kd ready (%.0fs)"%0, flush=True)
def Mf_mult(v): return np.exp(-2*bs)*expm_multiply(sp.linalg.LinearOperator((D,D),matvec=Wtot_mv,dtype=float) if False else Wtot_matrix(), v)
def Wtot_matrix():
    return sp.linalg.LinearOperator((D,D),matvec=Wtot_mv,dtype=float)
e0=np.zeros(D); e0[0]=1.0
w0=Mf_mult(e0); scale=float(np.dot(w0,Kd*w0))
print("scale=%.6e  (%.0fs)"%(scale,0), flush=True)
def T_mult(v): return Mf_mult(Kd*Mf_mult(v))/scale
op=sp.linalg.LinearOperator((D,D),matvec=T_mult,dtype=float)
t1=time.time()
ev,V=eigsh(op,k=2,which="LA",tol=1e-8,maxiter=2000)
o=np.argsort(ev)[::-1]; ev=ev[o]; V=V[:,o]
print("TWO 4-link plaquettes sharing a link, jm=%d, beta_s=%.2f : r=%.10f (%.0fs)"%(jm,bs,ev[1]/ev[0],time.time()-t1), flush=True)
