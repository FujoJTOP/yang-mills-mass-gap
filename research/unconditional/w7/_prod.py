import numpy as np, functools, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=2
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
    return ph*A
# single-link max eigenvector
U1=np.zeros((dl,dl))
for a in (-1,1): U1+=mult_op(1,a,a)
U1=0.5*U1
ev1,V1=np.linalg.eigh(0.5*(U1+U1.T))
vs=V1[:,-1]
print("jm=%d single-link max eig=%.10f"%(jm,ev1[-1]), flush=True)
print("  nontrivial components:", [(B[i],round(float(vs[i]),6)) for i in range(dl) if abs(vs[i])>1e-8][:10], flush=True)
# W_2 correct
tot=None
for a in (-1,1):
    for b in (-1,1):
        A=sp.csr_matrix(mult_op(1,a,b)); Bm=sp.csr_matrix(mult_op(1,a,b,conj=True))
        term=sp.kron(A,Bm,format='csr')
        tot=term if tot is None else tot+term
W2=(0.5*tot).tocsr()
D2=dl**2
ev2,V2=eigsh(W2,k=1,which="LA",tol=1e-12)
v2=V2[:,0]
# test the product ansatz: v_prod = v1 tensor v1 (in the correct order: link1 then link2)
prod=np.kron(vs,vs)
# the overlap |<prod|v2>|
ov=abs(np.dot(prod,v2))
print("W_2 max eig=%.10f"%ev2[0], flush=True)
print("  |<v1 (x) v1 | v2_max>| = %.10f   (1.0 means the product ansatz is exact)"%ov, flush=True)
# also test the Rayleigh quotient of the product
rq=float(prod @ (W2 @ prod))
print("  Rayleigh(prod) = %.10f"%(rq), flush=True)
