import numpy as np, functools
from scipy.special import iv
from scipy.linalg import eigh
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
jm=3
B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
Q=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        Q+=np.kron(mult_op(B,index,1,M2,N2,False),mult_op(B,index,1,M2,N2,True))
Q=(Q+Q.T)/2
bt=2.5
t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
tvec=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
Kd=np.diag(tvec)
M1=0.5*Q
T1=M1@Kd+Kd@M1; T1=(T1+T1.T)/2
ev=np.sort(tvec)[::-1]
print("beta_t=2.5 :  T(0)=K is diagonal ;  t_1=%.10f  t_2=%.10f"%(ev[0],ev[1]))
print("with the CORRECT M_1 = (1/2) Q :")
i1=int(np.argmax(tvec))
deg=[i for i in range(D) if abs(tvec[i]-ev[1])<1e-10*ev[0]]
print("   top index %d ; degenerate subspace size %d"%(i1,len(deg)))
print("   <e_1,T_1 e_1> = %+.10e"%(np.eye(D)[i1]@(T1@np.eye(D)[i1])))
W=np.eye(D)[:,deg]
A=W.T@T1@W
wa=np.sort(np.linalg.eigvalsh((A+A.T)/2))[::-1]
print("   T_1 restricted to the degenerate subspace: eigenvalues =", " ".join("%.8f"%v for v in wa[:6]))
print("   dln t_2/dbeta_s = lambda_max/t_2 = %.10f"%(wa[0]/ev[1]))
print("   dln t_1/dbeta_s = <e1,T1 e1>/t_1 = %.10f"%((np.eye(D)[i1]@(T1@np.eye(D)[i1]))/ev[0]))
print()
print("   => slope dr/dbeta_s|_0 (should be ~0.5*u = %.8f)"%(0.5*iv(2,bt)/iv(1,bt)))
