import numpy as np, functools, scipy.sparse as sp, time
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
def build_ops(jm):
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
    Ip=sp.identity(d,format='csr')
    ops={}
    for J2 in range(0,jm+1):
        ops[J2]=[sp.csr_matrix(np.kron(mult_op(J2,M2,N2),mult_op(J2,M2,N2,conj=True)))
                 for M2 in range(-J2,J2+1,2) for N2 in range(-J2,J2+1,2)]
    return B,ops,d
# center at jm=2 (d0=14), leaves at jm=1 (dl=5)
Bt,opst,d0=build_ops(2)
Bl,opsl,dl=build_ops(1)
bt=2.5; bs=1.0; k=6
td=np.zeros(3)
for j in range(3): td[j]=ive(j+1,bt)/ive(1,bt)
Ip0=sp.identity(d0,format='csr'); Ipl=sp.identity(dl,format='csr')
def q_pair(i_center_op, i_leaf_link):
    # q acting on (center, leaf_{i_leaf_link}) ; other leaves = identity
    tot=None
    for J2 in range(0,2):  # J2 = 0,1 for the q operator on a 2-link pair
        for op in opsl[J2] if False else []:
            pass
    # use the 2-link q operator built from dl x dl (leaves) and d0 x d0 (center)
    # build q on a pair with center dim d0 and leaf dim dl
    pass
# simpler: build the FULL 2-link q operator for a pair where both links are at jm=2, then restrict?
# Instead: build q on the (center, leaf) pair using mult_op on the leaf jm=1 space and on the center jm=2 space
# The q operator = sum_{M2,N2} kron(mult_op_center, mult_op_leaf^conj), but the CG conventions differ.
# Use the same construction as before but with the center at jm=2 and leaf at jm=1:
# We cannot mix jm in the standard mult_op. So restrict the center to jm=1 too, then take the k=6 star at jm=1 (already done), OR:
# take the center at jm=2 and leaves at jm=2 but truncate the LEAF index range by projection after building.
# Practical compromise: build the full jm=2 star (d=14) for k=6 in the LEAF-SYMMETRIC subspace.
print("Center at jm=2 / leaves at jm=1 requires a mixed CG construction not supported by the current mult_op.")
print("Fallback: use the leaf-symmetric subspace of the full jm=2, k=6 star.")
import math
for kk in (6,):
    Nsym=math.comb(14+kk-1,kk)
    print("   k=%d, d=14: symmetric leaf dimension C(%d,%d)=%d ; total D_sym=%d"%(kk,14+kk-1,kk,Nsym,14*Nsym))
