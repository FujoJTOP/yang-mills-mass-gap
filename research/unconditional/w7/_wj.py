import numpy as np, functools, scipy.sparse as sp
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
print("=== verify w_j = (2/(jm+2)) sin^2((j+1)pi/(jm+2)) ===")
print()
for jm in (1,2,3,4,5,6,7,8):
    B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
    dl=len(B)
    def mult_op(J2,M2,N2,conj=False):
        A=np.zeros((dl,dl))
        mm,nn=(-M2,-N2) if conj else (M2,N2)
        ph=((-1.0)**((M2-N2)//2)) if conj else 1.0
        for i2,(j2,m2,n2) in enumerate(B):
            m3,n3=m2-mm,n2-nn
            for j3 in range(abs(J2-j2),J2+j2+1,2):
                key=(j3,m3,n3)
                idx=index_b.get(key) if False else None
                # need index map
                pass
        return None
    index={b:i for i,b in enumerate(B)}
    def mult_op2(J2,M2,N2,conj=False):
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
    for a in (-1,1): U1+=mult_op2(1,a,a)
    U1=0.5*U1
    ev,V=np.linalg.eigh(0.5*(U1+U1.T))
    v=V[:,-1]
    ws=[]
    for jj in range(jm+1):
        s=sum(v[i]**2 for i in range(dl) if B[i][0]==jj)
        ws.append(s)
    ws=np.array(ws)
    pred=np.array([2.0/(jm+2)*np.sin((j+1)*np.pi/(jm+2))**2 for j in range(jm+1)])
    ok=np.allclose(ws,pred,atol=1e-8)
    print("jm=%d: measured w = %s"%(jm,np.round(ws,6)), flush=True)
    print("       predicted  = %s   match=%s"%(np.round(pred,6),ok), flush=True)
