import numpy as np, functools, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
def P2of(jm):
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
    U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
    Ui={(r,s):mult_op(1,r,s,conj=True) for r in (-1,1) for s in (-1,1)}
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        A=(U[(p,q)] @ Ui[(r,s)]).tocsr()
        Bm=(U[(q,r)] @ Ui[(s,p)]).tocsr()
        term=sp.kron(A,Bm,format='csr')
        acc=term if acc is None else acc+term
    return dl,(0.5*acc).tocsr()
print("=== standard plaquette P2 on 2 links: spectrum ===",flush=True)
print("%-4s %-5s %-8s %-16s %-16s %-16s"%("jm","d","D2","W_max(P2)","top2","top3"),flush=True)
lam={}
for jm in (1,2,3,4,5):
    dl,P2=P2of(jm)
    k=4
    ev=np.sort(eigsh(P2,k=k,which="LA",return_eigenvectors=False,tol=1e-11))[::-1]
    lam[jm]=ev
    print("%-4d %-5d %-8d %-16.10f %-16.10f %-16.10f"%(jm,dl,dl*dl,ev[0],ev[1],ev[2]),flush=True)
print()
print("=== compare lambda_1^triv with 6*W_max(P2) ===",flush=True)
triv={1:3.0000000000,2:3.5555555556,3:4.2548102767,4:4.6984325071}
print("%-4s %-16s %-16s %-16s %-12s"%("jm","lambda_1^triv","6*W_max(P2)","ratio","rel gap"),flush=True)
relg={1:0.37500000,2:0.0,3:0.01450709,4:0.00101591}
for jm in (1,2,3,4):
    w=lam[jm][0]
    print("%-4d %-16.10f %-16.10f %-16.10f %-12.8f"%(jm,triv[jm],6*w,triv[jm]/(6*w),relg[jm]),flush=True)
print()
print("=== W_max(P2) vs cos(pi/(jm+2)) (twisted formula) ===",flush=True)
for jm in (1,2,3,4,5):
    print("  jm=%d : W_max(P2)=%.10f   cos(pi/(jm+2))=%.10f   ratio=%.8f"%(
        jm,lam[jm][0],np.cos(np.pi/(jm+2)),lam[jm][0]/np.cos(np.pi/(jm+2))),flush=True)
print("[%.1fs]"%(time.time()-t0),flush=True)
