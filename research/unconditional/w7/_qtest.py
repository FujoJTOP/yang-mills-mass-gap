import numpy as np, functools, time
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def make_basis(jm):
    return [(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
def mult_op_correct(B,index,J2,M2,N2,conj=False):
    d=len(B); A=np.zeros((d,d))
    mm,nn=(-M2,-N2) if conj else (M2,N2)
    for i2,(j2,m2,n2) in enumerate(B):
        m3,n3=m2+mm,n2+nn
        for j3 in range(abs(J2-j2),J2+j2+1,2):
            key=(j3,m3,n3)
            if key not in index: continue
            cg1=d3(j2,m2,J2,mm,j3,m3)
            cg2=d3(j2,n2,J2,nn,j3,n3)
            v=cg1*cg2/(j3+1.0)
            if v!=0.0: A[index[key],i2]=v
    return A
jm=3
B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B); D=d*d
Mq=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        Mq+=np.kron(mult_op_correct(B,index,1,M2,N2,False),mult_op_correct(B,index,1,M2,N2,True))
Mq=(Mq+Mq.T)/2
ev=np.linalg.eigvalsh(Mq)
print("raw max|lam| = %.10f ; scale to norm 1: factor = %.10f"%(max(abs(ev[0]),abs(ev[-1])),1/max(abs(ev[0]),abs(ev[-1]))))
scale=1.0/max(abs(ev[0]),abs(ev[-1]))
Mq*=scale
Qop=np.eye(D)-Mq
# now check E_i >= 0 and the derivative identity, with beta_t=2.5 and the telescoped M_f
from scipy.special import iv
beta_t=2.5
t1d=2*iv(np.arange(0,jm+1)+1,beta_t)/beta_t
Kd=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
def mult_op_tel(B,index,J2,M2,N2,conj=False):
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
ops={}
for J2 in range(0,jm+1):
    ops[J2]=[np.kron(mult_op_tel(B,index,J2,M2,N2,False),mult_op_tel(B,index,J2,M2,N2,True))
             for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
def build(bs):
    Mf=np.zeros((D,D))
    for J2 in range(0,jm+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs)/bs)
        for op in ops[J2]: Mf+=cJ*op
    Mf*=np.exp(-bs)
    T=Mf@np.diag(Kd)@Mf; T=(T+T.T)/2
    ev,V=np.linalg.eigh(T); idx=np.argsort(ev)[::-1]
    return ev,V,idx
print()
print("%8s %14s %14s %14s %14s"%("beta_s","E1","E2","-2r(E2-E1)","r\x27(FD)"))
for bs in (0.8,1.2,2.0,3.0,4.0,8.0):
    ev,V,idx=build(bs); h=0.2
    x1=V[:,idx[0]]; x2=V[:,idx[1]]; r=ev[idx[1]]/ev[idx[0]]
    evp,_,idp=build(bs+h); evm,_,idm=build(bs-h)
    rp=((evp[idp[1]]/evp[idp[0]])-(evm[idm[1]]/evm[idm[0]]))/(2*h)
    E1=x1@(Qop@x1); E2=x2@(Qop@x2)
    print("%8.2f %14.8f %14.8f %14.8f %14.8f"%(bs,E1,E2,-2*r*(E2-E1),rp))
