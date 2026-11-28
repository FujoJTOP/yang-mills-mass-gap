import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=int(sys.argv[1]) if len(sys.argv)>1 else 1
bs=1.0; bt=2.5
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B)
def u_op(M2,N2,adjoint=False):
    A=np.zeros((d,d))
    for i2,(j2,m2,n2) in enumerate(B):
        mm,nn=(-M2,-N2) if adjoint else (M2,N2)
        for j3 in range(abs(1-j2),1+j2+1,2):
            m3,n3=m2+mm,n2+nn
            key=(j3,m3,n3)
            if key not in index: continue
            v=((-1.0)**((m2-n2)//2))*d3(j3,m3,1,mm,j2,-m2)*d3(j3,n3,1,nn,j2,-n2)
            if v!=0.0: A[index[key],i2]=np.sqrt((j3+1.0)*(j2+1.0))*v
    return A
# check normalization: sum over M2,N2 of u_op should give the character chi_{1/2} = tr(u)
I=np.eye(d)
# build W = (1/2) sum_{a,b,c,d} (u1)_{ab} (u2^-1)_{bc} (u3)_{cd} (u4^-1)_{da}
# u^-1 = u^dagger (SU(2)): (u^-1)_{bc} = conj(u_{cb}) = u_op(...,adjoint=True)_{bc}
ops={}
for M2 in (-1,1):
    for N2 in (-1,1):
        ops[(M2,N2)]=sp.csr_matrix(u_op(M2,N2))
        ops[(M2,N2,'adj')]=sp.csr_matrix(u_op(M2,N2,adjoint=True))
D=d**4
print("jm=%d d=%d D=%d"%(jm,d,D), flush=True)
W=None
t0=time.time()
for (a,b) in [(-1,-1),(-1,1),(1,-1),(1,1)]:
    A=ops[(a,b)]
    for (c,d_) in [(-1,-1),(-1,1),(1,-1),(1,1)]:
        Bm=ops[(-b,-c,'adj')] if False else ops[(c,-b,'adj')] if (c,-b) in ops else None
        # (u2^-1)_{bc} = conj(u2_{cb}) -> use adjoint with M2=c? Need the pair (b,c):
        Bm=ops[(b,c,'adj')]
        C=ops[(c,d_)]
        Dm=ops[(d_,a,'adj')]
        term=sp.kron(sp.kron(A,Bm,format='csr'),sp.kron(C,Dm,format='csr'),format='csr')
        W=term if W is None else W+term
W=(W*0.5).tocsr()
print("W built: nnz=%d  (%.0fs)"%(W.nnz,time.time()-t0), flush=True)
# check hermiticity
print("W hermitian ?", np.allclose((W-W.T).toarray() if D<1000 else (W-W.T).data, 0, atol=1e-10) if D<1000 else "skip")
# K = product of t_j
import itertools
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
Kd=np.ones(D)
for idx in range(D):
    x=idx; prod=1.0
    for _ in range(4):
        prod*=td[B[x%d][0]]; x//=d
    Kd[idx]=prod
print("K range: %.3e to %.3e"%(Kd.min(),Kd.max()), flush=True)
# M_f = exp(-beta_s(1-W)/2) = e^{-beta_s/2} exp(beta_s W / 2)
def Mf_mult(v): return np.exp(-2*bs)*expm_multiply((bs/2)*W, v)
e0=np.zeros(D); e0[0]=1.0
w0=Mf_mult(e0); scale=float(np.dot(w0,Kd*w0))
print("scale=%.6e"%scale, flush=True)
def T_mult(v): return Mf_mult(Kd*Mf_mult(v))/scale
op=LinearOperator((D,D),matvec=T_mult,dtype=float)
t1=time.time()
ev,V=eigsh(op,k=3,which="LA",tol=1e-11,maxiter=5000)
o=np.argsort(ev)[::-1]; ev=ev[o]; V=V[:,o]
print("4-link Wilson plaquette, jm=%d, beta_s=%.2f : t1=%.8f t2=%.8f r=%.10f  (%.0fs)"%(jm,bs,ev[0],ev[1],ev[1]/ev[0],time.time()-t1), flush=True)
for i in range(2):
    res=np.linalg.norm(T_mult(V[:,i])-ev[i]*V[:,i])/abs(ev[i])
    print("   eig %d residual = %.2e"%(i,res))
