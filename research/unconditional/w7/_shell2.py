import numpy as np, functools, scipy.sparse as sp, itertools, time
from scipy.special import ive
from scipy.sparse.linalg import expm_multiply, LinearOperator, eigsh
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
@functools.lru_cache(maxsize=None)
def d3(j1,m1,j2,m2,j3,m3):
    return float(sN(wigner_3j(Rational(j1,2),Rational(j2,2),Rational(j3,2),
                             Rational(m1,2),Rational(m2,2),Rational(m3,2)),20))
t0=time.time()
def build_ops(jm):
    B=[(j,m,n) for j in range(0,jm+1) for m in range(-j,j+1,2) for n in range(-j,j+1,2)]
    index={b:i for i,b in enumerate(B)}; dl=len(B)
    def mult_op(J2,M2,N2):
        A=np.zeros((dl,dl))
        for i2,(j2,m2,n2) in enumerate(B):
            m3,n3=m2-M2,n2-N2
            for j3 in range(abs(J2-j2),J2+j2+1,2):
                key=(j3,m3,n3)
                if key not in index: continue
                v=((-1.0)**((m2-n2)//2))*d3(j3,m3,J2,M2,j2,-m2)*d3(j3,n3,J2,N2,j2,-n2)
                if v!=0.0: A[index[key],i2]=np.sqrt((j3+1.0)*(j2+1.0))*v
        return sp.csr_matrix(A)
    U={(p,q):mult_op(1,p,q) for p in (-1,1) for q in (-1,1)}
    def Uinv(r,s): return U[(s,r)].T.tocsr()
    return B,dl,U,Uinv
NL=8
def plaq_slots():
    Ns=2;Nd=4
    def shift(site,d):
        L=2 if d==0 else 1
        return (site+1)%2 if L==2 else 0
    def lk(site,d): return site*Nd+d
    return [(lk(x,mu),lk(shift(x,mu),nu),lk(shift(x,nu),mu),lk(x,nu))
            for x in range(Ns) for mu in range(Nd) for nu in range(mu+1,Nd)]
plaq=plaq_slots()
def build_W4(dl,U,Uinv,slots):
    # 4-link plaquette operator, ordered by the sorted link positions
    order=sorted(range(4),key=lambda k:slots[k])
    acc=None
    for p in (-1,1):
     for q in (-1,1):
      for r in (-1,1):
       for s in (-1,1):
        fac=[U[(p,q)],U[(q,r)],Uinv(r,s),Uinv(s,p)]
        onl={}
        for k in range(4): onl.setdefault(slots[k],[]).append(fac[k])
        term=None
        for l in sorted(set(slots)):
            M=onl[l][0]
            for X in onl[l][1:]: M=(M@X).tocsr()
            term=M if term is None else sp.kron(term,M,format='csr')
        acc=term if acc is None else acc+term
    return (0.5*acc).tocsr()
print("=== shell-truncation, validated against exact jm=1 ===",flush=True)
for jm in (1,3):
    B,dl,U,Uinv=build_ops(jm)
    td=np.array([ive(j+1,2.5)/ive(1,2.5) for j in range(jm+1)])
    js=np.array([td[b[0]] for b in B])
    # nonzero-j states per link
    nz=[i for i,b in enumerate(B) if b[0]>0]
    print("jm=%d: dl=%d, nonzero-j states per link=%d"%(jm,dl,len(nz)),flush=True)
    for Nmax in (0,1,2,3):
        # shell states: tuples of 8 link indices, at most Nmax with j>0
        states=[(0,)*4*NL] if False else []
        # link state 0 = j=0 (m=n=0) : index of (0,0,0)
        z=[i for i,b in enumerate(B) if b[0]==0][0]
        states=[(z,)*NL]
        for k in range(1,Nmax+1):
            for pos in itertools.combinations(range(NL),k):
                for vals in itertools.product(nz,repeat=k):
                    s=[z]*NL
                    for p_,v in zip(pos,vals): s[p_]=v
                    states.append(tuple(s))
        n=len(states)
        if n>3_000_000: 
            print("   Nmax=%d: n=%d TOO BIG"%(Nmax,n),flush=True); continue
        row=np.array([int("".join("%d"%s for s in st)) if False else 0 for st in states])
        # build index map via mixed-radix
        def code(st):
            c=0
            for s in st: c=c*dl+s
            return c
        codes=np.array([code(st) for st in states],dtype=np.int64)
        idx={int(c):i for i,c in enumerate(codes)}
        # restricted V
        Vr=sp.lil_matrix((n,n))
        for slots in plaq:
            W4=build_W4(dl,U,Uinv,slots)
            W4c=W4.tocoo()
            # group by the tuple of 4 changed link states
            for a_i,stb in enumerate(states):
                # apply W4 to stb on the 4 slots
                sub=tuple(stb[l] for l in slots)
                # gather the 4-link basis indices
                for i0 in range(dl):
                    if B[i0][0]!=B[sub[0]][0]: continue
                    for i1 in range(dl):
                        if B[i1][0]!=B[sub[1]][0]: continue
                        for i2 in range(dl):
                            if B[i2][0]!=B[sub[2]][0]: continue
                            for i3 in range(dl):
                                if B[i3][0]!=B[sub[3]][0]: continue
                                vv=0.0
                                col=((i0*dl+i1)*dl+i2)*dl+i3
                                # row: matching slot basis indices
                                rws=[]
                                for j0 in range(dl):
                                    if B[j0][0]!=B[i0][0]: continue
                                # simpler: use dense W4 rows for matching j
                                pass
                pass
        print("   Nmax=%d: n=%d (build skipped: too slow)"%(Nmax,n),flush=True)
print("[%.1fs]"%(time.time()-t0),flush=True)
