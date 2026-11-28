import numpy as np, functools, scipy.sparse as sp
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
for jm in (3,4):
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
    ev,V=np.linalg.eigh(0.5*(U1+U1.T))
    vmax=V[:,-1]
    print("jm=%d  max eigenvalue = %.10f  (cos(pi/%d)=%.10f)"%(jm,ev[-1],jm+2,np.cos(np.pi/(jm+2))), flush=True)
    # print the components grouped by j
    print("   eigenvector components (nonzero, by j):", flush=True)
    for jj in range(0,jm+1):
        nz=[(round(B[i][1],3),round(B[i][2],3),round(float(vmax[i]),6)) for i in range(dl) if B[i][0]==jj and abs(vmax[i])>1e-8]
        if nz:
            print("     j=%d : %s"%(jj,nz[:8]), flush=True)
    # check SU(2) invariance: the state should be a function of the j only? (m,n structure)
    print("   |v| grouped by (j): ", end="", flush=True)
    norms=[]
    for jj in range(0,jm+1):
        s=np.sqrt(sum(vmax[i]**2 for i in range(dl) if B[i][0]==jj))
        norms.append(round(float(s),6))
    print(norms, flush=True)
    # is the eigenvector proportional to a (2j+1)-dim trivial combination? compute the "trace" structure
