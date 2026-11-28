import numpy as np, functools, time
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def make_basis(jm):
    return [(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
# <j'm'n'| D^{j0}_{MN} |jmn> = 1/(2j'+1) * <j0 M, j m | j' m'> * <j0 N, j n | j' n'>
def mult_op_correct(B,index,J2,M2,N2,conj=False):
    d=len(B); A=np.zeros((d,d))
    mm,nn=(-M2,-N2) if conj else (M2,N2)
    for i2,(j2,m2,n2) in enumerate(B):
        m3,n3=m2+mm,n2+nn
        for j3 in range(abs(J2-j2),J2+j2+1,2):
            key=(j3,m3,n3)
            if key not in index: continue
            # <j0 M, j m | j' m'>  with j0=J2, M=mm
            cg1=d3(j2,m2,J2,mm,j3,m3)   # wigner_3j(j2,m2,J2,mm,j3,m3) ~ <j2 m2 J2 mm|j3 m3>
            cg2=d3(j2,n2,J2,nn,j3,n3)
            v=cg1*cg2/(j3+1.0)
            if v!=0.0: A[index[key],i2]=v
    return A
print("check: ||M_q||_op should be 1 for q=(1/2)chi_{1/2}")
print("%8s %8s %10s %14s %14s"%("jmax2","d","D","max|lam|","min lam"))
for jm in (2,3,4):
    t0=time.time()
    B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
    Mq=np.zeros((D,D))
    for M2 in range(-1,2,2):
        for N2 in range(-1,2,2):
            Mq+=0.5*np.kron(mult_op_correct(B,index,1,M2,N2,False),mult_op_correct(B,index,1,M2,N2,True))
    Mq=(Mq+Mq.T)/2
    ev=np.linalg.eigvalsh(Mq)
    print("%8d %8d %10d %14.8f %14.8f  [%.1fs]"%(jm,d,D,ev[-1],ev[0],time.time()-t0), flush=True)
