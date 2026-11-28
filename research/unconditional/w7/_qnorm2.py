import numpy as np, functools, time
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def make_basis(jm):
    return [(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
def mult_op(B,index,J2,M2,N2,norm,conj=False):
    d=len(B); A=np.zeros((d,d))
    mm,nn=(-M2,-N2) if conj else (M2,N2)
    ph=((-1.0)**((M2-N2)//2)) if conj else 1.0
    for i2,(j2,m2,n2) in enumerate(B):
        m3,n3=m2-mm,n2-nn
        for j3 in range(abs(J2-j2),J2+j2+1,2):
            key=(j3,m3,n3)
            if key not in index: continue
            v=((-1.0)**((m2-n2)//2))*d3(j3,m3,J2,mm,j2,-m2)*d3(j3,n3,J2,nn,j2,-n2)
            if v!=0.0: A[index[key],i2]=norm(j3,j2,J2)*v
    return ph*A
def nstd(j3,j2,J2): return np.sqrt((j3+1.0)/(j2+1.0))
def ntel(j3,j2,J2): return np.sqrt((j3+1.0)*(j2+1.0))
print("dense check with eigvalsh (max|lambda|, NOT SVD)")
for jm in (2,3,4):
    t0=time.time()
    B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
    for tag,nm in (("std",nstd),("tel",ntel)):
        Mq=np.zeros((D,D))
        for M2 in range(-1,2,2):
            for N2 in range(-1,2,2):
                Mq+=np.kron(mult_op(B,index,1,M2,N2,nm,False),mult_op(B,index,1,M2,N2,nm,True))
        Mq=(Mq+Mq.T)/2
        ev=np.linalg.eigvalsh(Mq)
        print("  jmax2=%d d=%d D=%d  %s : max=%.8f min=%.8f  [%.1fs]"%(jm,d,D,tag,ev[-1],ev[0],time.time()-t0), flush=True)
