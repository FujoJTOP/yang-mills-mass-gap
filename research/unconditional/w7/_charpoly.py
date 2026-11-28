import numpy as np, functools, scipy.sparse as sp
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
print("=== test: char poly of (1/2)tr(u) on j<=jm vs Chebyshev U_{jm+1}(2*lambda) / 2^(jm+1) ===")
print()
for jm in (1,2,3,4,5):
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
    # char poly coefficients of the full matrix
    cf=np.poly(U1)   # polynomial with coefficients
    # Chebyshev U_{jm+1}(2 lambda): zeros at cos(k pi/(jm+2))
    # U_n(cos t) = sin((n+1)t)/sin t ; zeros at t=k pi/(n+1), k=1..n
    # So we want U_{jm+1}(x) with x=2 lambda? Actually (1/2)tr(u) is not the same as cos t.
    # Instead compare the zeros of cf with cos(k pi/(jm+2))
    roots=np.roots(cf)
    roots=np.sort(roots[np.abs(roots.imag)<1e-8].real)
    target=np.sort(np.cos(np.arange(1,jm+2)*np.pi/(jm+2)))
    # only the mult-1 part will match exactly
    print("jm=%d dl=%d"%(jm,dl), flush=True)
    print("   distinct roots (real):", np.round(np.unique(np.round(roots,6)),6), flush=True)
    print("   cos(k pi/(jm+2))     :", np.round(np.unique(np.round(target,6)),6), flush=True)
