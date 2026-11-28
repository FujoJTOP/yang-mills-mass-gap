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
# fundamental u_{ab}: J2=1 (doubled) i.e. j=1/2
u={}
for M2 in (-1,1):
    for N2 in (-1,1):
        u[(M2,N2)]=sp.csr_matrix(mult_op(1,M2,N2))
        u[(M2,N2,'c')]=sp.csr_matrix(mult_op(1,M2,N2,conj=True))
D=d**4
print("jm=%d d=%d D=%d nnz(u_{{00}})=%d"%(jm,d,D,u[(-1,-1)].nnz), flush=True)
# W = (1/2) sum_{a,b,c,d} (u1)_{ab} (u2^{-1})_{bc} (u3)_{cd} (u4^{-1})_{da}
# u^{-1} = conj transpose. u2^{-1} acts on link 2. Its matrix element (bc) = conj(u2_{cb}).
# Using mult_op(...,conj=True) for the adjoint (which for SU(2) = inverse)
W=None
idx_list=[(-1,-1),(-1,1),(1,-1),(1,1)]
t0=time.time()
for (a,b) in idx_list:
    A=u[(a,b)]
    for (c,dd) in idx_list:
        Bm=u[(b,c,'c')]   # (u2^{-1})_{bc}
        C=u[(c,dd)]
        Dm=u[(dd,a,'c')]  # (u4^{-1})_{da}
        term=sp.kron(sp.kron(A,Bm,format='csr'),sp.kron(C,Dm,format='csr'),format='csr')
        W=term if W is None else W+term
W=((W*0.5)+(W*0.5).T*0).tocsr()
W=0.5*W
print("W nnz=%d (%.0fs)"%(W.nnz,time.time()-t0), flush=True)
# hermiticity check on a small subspace
if D<=625:
    Wa=W.toarray(); print("  max|W-W^T| = %.3e  max|W| = %.3e"%(np.abs(Wa-Wa.T).max(), np.abs(Wa).max()), flush=True)
td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
Kd=np.ones(D)
for idx in range(D):
    x=idx; prod=1.0
    for _ in range(4):
        prod*=td[B[x%d][0]]; x//=d
    Kd[idx]=prod
def Mf_mult(v): return np.exp(-2*bs)*expm_multiply((bs/2)*W, v)
e0=np.zeros(D); e0[0]=1.0
w0=Mf_mult(e0); scale=float(np.dot(w0,Kd*w0))
print("scale=%.6e"%scale, flush=True)
def T_mult(v): return Mf_mult(Kd*Mf_mult(v))/scale
op=LinearOperator((D,D),matvec=T_mult,dtype=float)
t1=time.time()
ev,V=eigsh(op,k=3,which="LA",tol=1e-11,maxiter=5000)
o=np.argsort(ev)[::-1]; ev=ev[o]; V=V[:,o]
print("4-link Wilson plaquette, jm=%d, beta_s=%.2f : t1=%.10f t2=%.10f r=%.10f  (%.0fs)"%(jm,bs,ev[0],ev[1],ev[1]/ev[0],time.time()-t1), flush=True)
for i in range(2):
    res=np.linalg.norm(T_mult(V[:,i])-ev[i]*V[:,i])/abs(ev[i])
    print("   eig %d residual = %.2e"%(i,res))
