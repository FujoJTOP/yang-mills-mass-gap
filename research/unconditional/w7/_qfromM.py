import numpy as np, functools, time
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def make_basis(jm):
    return [(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
def mult_op(B,index,J2,M2,N2,conj=False):
    d=len(B); A=np.zeros((d,d))
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
from scipy.special import iv
jm=3
B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
def Mf_of(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs)/bs)
        for op in ops[J2]: M+=cJ*op
    return M*np.exp(-bs)
bs=2.0; h=1e-4
Mp=Mf_of(bs+h); Mm=Mf_of(bs-h); M0=Mf_of(bs)
dM=(Mp-Mm)/(2*h)
Qmat=-np.linalg.solve(M0,dM)     # Q = -M_f^{-1} dM_f/dbeta_s
Qmat=(Qmat+Qmat.T)/2
ev=np.linalg.eigvalsh(Qmat)
print("Q from  Q = -M_f^{-1} dM_f/dbeta_s   (beta_s=2.0)")
print("  ||Q||_op = %.8f   min = %.8f   max = %.8f"%(max(abs(ev[0]),abs(ev[-1])),ev[0],ev[-1]))
print("  (Q = 1 - q should have spectrum in [0,2], ||Q|| <= 2)")
