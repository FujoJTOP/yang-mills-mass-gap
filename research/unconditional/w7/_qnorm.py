import numpy as np, functools, time
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN

@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))

def make_basis(jm):
    return [(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]

def build_norm(jm, kind):
    B=make_basis(jm); index={b:i for i,b in enumerate(B)}; d=len(B)
    rows=[];cols=[];vals=[]
    for M2 in range(-1,2,2):
        for N2 in range(-1,2,2):
            for conj in (False,True):
                mm,nn=(-M2,-N2) if conj else (M2,N2)
                ph=((-1.0)**((M2-N2)//2)) if conj else 1.0
                # single-link block
                for i2,(j2,m2,n2) in enumerate(B):
                    m3,n3=m2-mm,n2-nn
                    for j3 in range(abs(1-j2),1+j2+1,2):
                        key=(j3,m3,n3)
                        if key not in index: continue
                        v=((-1.0)**((m2-n2)//2))*d3(j3,m3,1,mm,j2,-m2)*d3(j3,n3,1,nn,j2,-n2)
                        if v==0.0: continue
                        a = np.sqrt((j3+1.0)*(j2+1.0)) if kind=="tel" else np.sqrt((j3+1.0)/(j2+1.0))
                        rows.append(index[key]); cols.append(i2); vals.append(ph*a*v)
    A=coo_matrix((vals,(rows,cols)),shape=(d,d)).toarray()
    # Mq = sum kron(A_{MN}, Ac_{MN})   -> build sparse kron
    from scipy.sparse import kron as skron
    M=0
    for M2 in range(-1,2,2):
        for N2 in range(-1,2,2):
            mm,nn=(M2,N2)
            # rebuild the two single-link blocks (non-conj and conj)
            def blk(mmm,nnn,ph):
                r=[];c=[];vv=[]
                for i2,(j2,m2,n2) in enumerate(B):
                    m3,n3=m2-mmm,n2-nnn
                    for j3 in range(abs(1-j2),1+j2+1,2):
                        key=(j3,m3,n3)
                        if key not in index: continue
                        v=((-1.0)**((m2-n2)//2))*d3(j3,m3,1,mmm,j2,-m2)*d3(j3,n3,1,nnn,j2,-n2)
                        if v==0.0: continue
                        a = np.sqrt((j3+1.0)*(j2+1.0)) if kind=="tel" else np.sqrt((j3+1.0)/(j2+1.0))
                        r.append(index[key]); c.append(i2); vv.append(ph*a*v)
                return coo_matrix((vv,(r,c)),shape=(d,d))
            A1=blk(mm,nn,1.0)
            A2=blk(-mm,-nn,((-1.0)**((M2-N2)//2)))
            M = M + skron(A1,A2,format="csr")
    return M, d

for jm in (2,3,4):
    t0=time.time()
    M,d = build_norm(jm,"std")
    D=d*d
    vmax=eigsh(M.tocsr(),k=1,which="LA",return_eigenvectors=False,maxiter=5000)
    vmin=eigsh(M.tocsr(),k=1,which="SA",return_eigenvectors=False,maxiter=5000)
    print("jmax2=%d d=%d D=%d  ||M_q||_max=%.8f  min=%.8f  [%.1fs]"%(jm,d,D,vmax[0],vmin[0],time.time()-t0), flush=True)
