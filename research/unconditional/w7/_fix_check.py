import numpy as np, functools, scipy.sparse as sp, time
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
# ---- CORRECT mult_op with conj, following the ORIGINAL session convention ----
def make(jm):
    B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
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
        return sp.csr_matrix(ph*A)
    return B,dl,mult_op
jm=1
B,dl,mult_op=make(jm)
U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
Ui={(r,s):mult_op(1,r,s,conj=True) for r in (-1,1) for s in (-1,1)}
I=sp.identity(dl,format='csr')
# sanity: (u^{-1})_{rs} should be the inverse of the full operator u
# test: sum_r u_{ar} u^{-1}_{rb} = delta_ab ?
tot=None
for r in (-1,1):
    pass
print("sanity checks on mult_op conventions:",flush=True)
# u^{-1} as a matrix: build Uop[i,j] and check inverse
Uop=np.zeros((dl,dl))
for (p,q),M in U.items():
    # Uop should be the full operator; check sum over a of u_ab corresponds
    pass
# better sanity: <0| V |0> and V|0> structure
def pair_op(mu,nu):
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        A=(U[(p,q)] @ Ui[(r,s)]).tocsr()
        Bm=(U[(q,r)] @ Ui[(s,p)]).tocsr()
        term=None
        for l in range(4):
            M = A if l==mu else (Bm if l==nu else I)
            term=M if term is None else sp.kron(term,M,format='csr')
        acc=term if acc is None else acc+term
    return (0.5*acc).tocsr()
P=pair_op(0,1)
D=dl**4
e0=np.zeros(D); e0[0]=1.0
Pv=P@e0
print("pair_op(0,1)|0> : nnz=%d  norm=%.6f"%(np.count_nonzero(np.abs(Pv)>1e-14),np.linalg.norm(Pv)),flush=True)
V=None
for mu in range(4):
    for nu in range(mu+1,4):
        Q=pair_op(mu,nu); V=Q if V is None else (V+Q).tocsr()
V=(0.5*(V+V.T)).tocsr()
Vv=V@e0
print("V|0> : nnz=%d  norm=%.6f"%(np.count_nonzero(np.abs(Vv)>1e-14),np.linalg.norm(Vv)),flush=True)
w=float(eigsh(V,k=1,which="LA",return_eigenvectors=False,tol=1e-10)[0])
w1=float(eigsh(P,k=1,which="LA",return_eigenvectors=False,tol=1e-10)[0])
print("CORRECTED: W_max(V)=%.10f   6*W_max(1p)=%.10f  W_max(1p)=%.10f"%(w,6*w1,w1),flush=True)
