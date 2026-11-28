import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
from scipy.sparse.linalg import LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=int(sys.argv[1]) if len(sys.argv)>1 else 2
bs=float(sys.argv[2]) if len(sys.argv)>2 else 1.0
bt=2.5
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
W4a=Wilson4(0,1,2,3); W4b=Wilson4(0,4,5,6)
shape7=(d,)*7
def W4a_mv(v):
    T=v.reshape(d**4, d**3); return (W4a @ T).ravel()
def W4b_mv(v):
    T=v.reshape(shape7)
    perm=(0,4,5,6,1,2,3)
    Tp=np.transpose(T,perm).reshape(d**4, d**3)
    out=(W4b @ Tp).reshape((d,)*7)
    inv=np.argsort(perm)
    return np.transpose(out,inv).ravel()
def Wtot_mv(v): return 0.5*(W4a_mv(v)+W4b_mv(v))
D=d**7
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
K1d=np.array([td[b[0]] for b in B])
Kd=np.ones(D)
tmp=K1d
for _ in range(6):
    tmp=np.kron(tmp,K1d)
Kd=tmp
print("jm=%d d=%d D=%d  (%.0fs)"%(jm,d,D,0), flush=True)
# use scipy expm_multiply with normalize trick to avoid onenormest: precompute norm bound
# |Wtot| <= 0.5*(|W1|+|W2|) <= 0.5*0.5*2 = 0.5 -> Taylor with 20 terms is plenty
def Mf_mult(v):
    out=v.copy(); term=v.copy()
    for k in range(1,25):
        term=Wtot_mv(term)*(bs/k)
        out=out+term
    return np.exp(-2*bs)*out
def T_mult(v): return Mf_mult(Kd*Mf_mult(v))/scale
e0=np.zeros(D); e0[0]=1.0
t0=time.time()
w0=Mf_mult(e0); scale=float(np.dot(w0,Kd*w0))
print("scale=%.6e (%.0fs)"%(scale,time.time()-t0), flush=True)
op=LinearOperator((D,D),matvec=T_mult,dtype=float)
t1=time.time()
try:
    ev,V=eigsh(op,k=2,which="LA",tol=1e-4,maxiter=12)
    o=np.argsort(ev)[::-1]; ev=ev[o]
    print("LANCZOS(12) : t1=%.10f t2=%.10f r=%.10f (%.0fs)"%(ev[0],ev[1],ev[1]/ev[0],time.time()-t1), flush=True)
except Exception as ex:
    print("eigsh failed:",ex, flush=True)
