import numpy as np, functools
from scipy.special import iv
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
# ORTHONORMAL rescaling:  basis  v_{jmn} = sqrt(2j+1) D^j_{mn}   (unit L^2 norm)
w=np.array([np.sqrt((b[0]+1)*(c[0]+1)) for b in B for c in B])   # 2j+1 = j+1 in doubled labels
W=np.diag(1.0/w)
beta_t=2.5
t1d=2*iv(np.arange(0,jm+1)+1,beta_t)/beta_t
Kd=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
def r_of(bs, weighted):
    Mf=np.zeros((D,D))
    for J2 in range(0,jm+1):
        cJ=(1.0 if J2==0 else 2*(J2+1)*iv(J2+1,bs)/bs)
        for op in ops[J2]: Mf+=cJ*op
    Mf*=np.exp(-bs)
    if weighted:
        T=W@Mf@np.diag(Kd)@Mf@W
    else:
        T=Mf@np.diag(Kd)@Mf
    T=(T+T.T)/2
    ev=np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[1]/ev[0]
print("r(beta_s)  raw (current code)  vs  orthonormal-basis weight 1/sqrt((2j1+1)(2j2+1))")
print("%8s %18s %18s"%("beta_s","raw","weighted"))
for bs in (0.0,0.5,1.0,2.0,4.0,8.0):
    a=r_of(bs,False); b=r_of(bs,True)
    print("%8.2f %18.10f %18.10f"%(bs,a,b))
