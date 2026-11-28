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
bt=2.5
t1d=2*iv(np.arange(0,jm+1)+1,bt)/bt
tvec=np.array([t1d[b1[0]]*t1d[b2[0]] for b1 in B for b2 in B])
# M_1 = Op_1 = q operator (J2=1 sum)
M1=np.zeros((D,D))
for M2 in range(-1,2,2):
    for N2 in range(-1,2,2):
        M1+=np.kron(mult_op(B,index,1,M2,N2,False),mult_op(B,index,1,M2,N2,True))
M1=(M1+M1.T)/2
idx1=int(np.argmax(tvec))     # e_1 = the j=0 (x) j=0 state
t1v=tvec[idx1]
e1=np.zeros(D); e1[idx1]=1.0
print("e_1 is index %d with t = %.10f"%(idx1,t1v))
print("   <e1,M1 e1>          = %+.10e"%(e1@(M1@e1)))
print("   <e1,M1^2 e1>        = %+.10e"%(e1@(M1@(M1@e1))))
v=M1@e1
nz=[(j,v[j]) for j in range(D) if abs(v[j])>1e-12]
print("   M1 e1 has %d nonzero components:"%len(nz))
for j,val in nz[:8]:
    print("      component %d : %.10f   t_j = %.6f   |v|^2*t_j/(t1-t_j) = %+.10f"%(j,val,tvec[j],val*val*tvec[j]/(t1v-tvec[j])))
S=sum(val*val*tvec[j]/(t1v-tvec[j]) for j,val in nz)
print()
print("   sum_j t_j|<e_j,M1 e1>|^2/(t1-t_j) = %+.10f"%S)
print("   <e1,M2 e1> = <e1,M1^2 e1> = %.10f (since M2 = q^2)"%(e1@(M1@(M1@e1))))
lam1p=2*(e1@(M1@(M1@e1)))+2*(e1@(M1@(M1@e1)))+2*S
print("   lambda_1' = 2<M2>+2<M1^2>+2*sum = %.10f     (measured 1.19268459)"%lam1p)
