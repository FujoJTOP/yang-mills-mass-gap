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
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op(B,index,J2,M2,N2,False),mult_op(B,index,J2,M2,N2,True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
def top2(bt,bs):
    td=2*iv(np.arange(0,jm+1)+1,bt)/bt
    Kd=np.diag(np.array([td[b1[0]]*td[b2[0]] for b1 in B for b2 in B]))
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    T=M@Kd@M; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[1]/ev[0]
print("fit  ln r(beta_s)/beta_s^2  ->  lambda' (should be lambda_2' - lambda_1')")
print("%8s %12s %18s %18s"%("beta_t","u","lambda_2'-lambda_1'","lambda_1' (closed form)"))
for bt in (1.0,2.0,2.5,4.0,8.0,16.0):
    u=iv(2,bt)/iv(1,bt); u2=u*u
    hs=np.array([2e-4,5e-4,1e-3,2e-3])
    ys=np.array([np.log(top2(bt,h))/h**2 for h in hs])
    A=np.vstack([hs**2]).T
    lam2m1=np.mean(ys)     # dln r/dbs^2 = lambda_2' - lambda_1'
    lam1p=0.5*(1+3*u2)/(1-u2)
    print("%8.1f %12.8f %18.10f %18.10f"%(bt,u,lam2m1,lam1p), flush=True)
