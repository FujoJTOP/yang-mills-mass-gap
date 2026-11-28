import numpy as np, functools, scipy.sparse as sp, time, sys
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
jm=int(sys.argv[1]) if len(sys.argv)>1 else 2
bs=1.0; bt=2.5
B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
index={b:i for i,b in enumerate(B)}; d=len(B)
def mult_op(J2,M2,N2,conj=False):
    A=np.zeros((d,d))
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
# REAL 4-link plaquette: W = (1/2) tr(u1 u2^-1 u3 u4^-1)
# This is a 4-link operator. Build it as a product of 2-link q operators? No - it is a single 4-link operator.
# Use the character expansion: e^{(beta/2) W} = sum_J c_J chi_J(W). But W is a 4-link operator.
# Instead, build the operator on 4 links directly:
#   (1/2) tr(u1 u2^-1 u3 u4^-1) - build via the 2-link q operators? No.
# Practical: build the operator as a sum over M2,N2 of the 4-link kernel.
# The 4-link Wilson loop can be written as a product of two 2-link "links": A = u1 u2^-1, B = u3 u4^-1, W = (1/2) tr(A B^-1).
# Since A and B live on disjoint link sets, W is a 2-link coupling in the (A,B) variables.
# But A,B are themselves 2-link operators. The basis is 4 links. 
# Build W via: W = (1/2) sum_{a,b} A_{ab} B^{-1}_{ba}, i.e., a 2-link "contraction" of the 2-link operators A and B.
# This is a genuine 4-link operator, but can be assembled from the 2-link kernels.
print("The 4-link Wilson plaquette requires a genuine 4-link operator construction.")
print("Its dimension is d^4 =", d**4, "for jm =", jm)
print("A direct construction via the (A,B) 2-link kernels is possible but needs care.")
