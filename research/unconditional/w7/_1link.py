import numpy as np, functools, scipy.sparse as sp
from scipy.special import ive
from scipy.linalg import expm
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
print("=== 1-link model: exact T in the j-subspace (dim jm+1) ===")
print("    W|_j = (1/2) A_path ; K|_j = diag(t_j) ; T = M_f K M_f")
print()
bt=2.5
for jm in (3,5,7,10,15,20):
    # A_path on jm+1 vertices
    nJ=jm+1
    Ap=np.zeros((nJ,nJ))
    for k in range(nJ-1):
        Ap[k,k+1]=0.5; Ap[k+1,k]=0.5
    # We need W|_j = (1/2) A_path -> already Ap = (1/2) A_path (entries 0.5)
    td=np.array([ive(j+1,bt)/ive(1,bt) for j in range(jm+1)])
    Kd=np.diag(td)
    for bs in (0.5,1.0,2.0,4.0):
        Mf=expm(bs*Ap)   # exp(beta_s * (1/2) A_path) ; the full M_f = exp(beta_s W/2) = exp(beta_s*(1/2)A_path/2)?
        # careful: M_f = exp(-(beta_s/2)*4) exp(beta_s W) for 4 links; for 1 link: M_f = exp(beta_s W/2)
        # here W = W|_j = (1/2) A_path, so M_f = exp(beta_s * (1/2) A_path / 2) = exp(beta_s * A_path / 4)
        Mf=expm(bs*Ap/2)   # = exp(beta_s * A_path /4) since Ap=(1/2)A -> bs*Ap/2 = bs*A/4
        w0=Mf[:,0]
        sc=float(w0 @ (td*w0))
        T=(Mf @ Kd @ Mf)/sc
        ev=np.linalg.eigvalsh(0.5*(T+T.T))[::-1]
        r=ev[1]/ev[0]
        print("jm=%2d beta_s=%.1f : t1=%.8f t2=%.8f r=%.8f  m a=%.6f"%(
            jm,bs,ev[0],ev[1],r,-np.log(r)), flush=True)
    print()
