import numpy as np, functools, scipy.sparse as sp
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
print("=== verify: multiplicity-1 spectrum of (1/2)tr(u) on j<=jm equals")
print("    pm cos(k pi/(jm+2)) for k=1..floor((jm+1)/2) ===")
print()
for jm in (1,2,3,4,5,6):
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
        return ph*A
    U1=np.zeros((dl,dl))
    for a in (-1,1):
        U1 += mult_op(1,a,a)
    U1=0.5*U1
    ev=np.linalg.eigvalsh(0.5*(U1+U1.T))
    vals,counts=np.unique(np.round(ev,8),return_counts=True)
    mult1=[v for v,c in zip(vals,counts) if c==1]
    # predicted
    j=j m if False else jm
    pred=[]
    kmax=(jm+1)//2
    for k in range(1,kmax+1):
        pred.append(np.cos(k*np.pi/(jm+2)))
        pred.append(-np.cos(k*np.pi/(jm+2)))
    mult1=np.sort(np.array(mult1)); pred=np.sort(np.array(pred))
    ok = len(mult1)==len(pred) and np.allclose(mult1,pred,atol=1e-8)
    print("jm=%d (dl=%2d, jm+2=%d): mult-1 values n=%d, predicted n=%d, match=%s"%(
        jm,dl,jm+2,len(mult1),len(pred),ok), flush=True)
    if not ok:
        print("   measured:",np.round(mult1,8), flush=True)
        print("   predicted:",np.round(pred,8), flush=True)
