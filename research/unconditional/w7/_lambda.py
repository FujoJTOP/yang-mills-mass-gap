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
bt=2.5
t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
Kd=np.diag(np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B]))
def top2(bs):
    M=np.zeros((D,D))
    for J2 in range(0,jm+1):
        if J2==0: cJ=1.0
        elif bs==0.0: continue
        else: cJ=2*(J2+1)*iv(J2+1,bs)/bs
        for op in ops[J2]: M+=cJ*op
    T=M@Kd@M; T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[0],ev[1]
hs=np.array([5e-5,1e-4,2e-4,3e-4,5e-4])
t10,t20=top2(0.0)
print("beta_t=2.5 :  t1(0)=%.10f  t2(0)=%.10f"%(t10,t20))
l1=[];l2=[]
for h in hs:
    a,b=top2(h)
    l1.append(np.log(a/t10)); l2.append(np.log(b/t20))
l1=np.array(l1); l2=np.array(l2)
A=np.vstack([hs,hs**2]).T
c1a,_,_,_=np.linalg.lstsq(A,l1,rcond=None)
c1b,_,_,_=np.linalg.lstsq(A,l2,rcond=None)
print("   ln(t1/t1(0)) = %.8f * bs + %.8f * bs^2"%(c1a[0],c1a[1]))
print("   ln(t2/t2(0)) = %.8f * bs + %.8f * bs^2"%(c1b[0],c1b[1]))
print()
print("   so  lambda_1 = %.8f , lambda_1' = %.8f"%(c1a[0],2*c1a[1]))
print("       lambda_2 = %.8f , lambda_2' = %.8f"%(c1b[0],2*c1b[1]))
print("       lambda_2-lambda_1 = %.8f   lambda_2'-lambda_1' = %.8f"%(c1b[0]-c1a[0],2*(c1b[1]-c1a[1])))
print()
print("   and c_2 = [ r0*lambda_2' + (r0/2)*(lambda_2-lambda_1)^2 ... ]  check:")
r0=t20/t10
print("   c_1 = %.10f   (should be 0.5)"%(r0*(c1b[0]-c1a[0])))
print("   d2r/dbs2|0 = 2*r0*c_2 ; c_2 from lambda's = %.10f"%(r0*( (c1b[0]-c1a[0])**2/2 + (c1b[1]-c1a[1]) )))
