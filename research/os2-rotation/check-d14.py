# Exact integer check for the Belinfante improvement of the l=4 massive Gaussian.
# T_{mu nu} = p_mu p_nu - lam p_mu p_nu^3 - delta_{mu nu} L   (canonical, higher-derivative)
# A_{mu nu} = T_{mu nu} - T_{nu mu} = -lam (p_mu p_nu^3 - p_nu p_mu^3)
# Improvement: T' = T + d_rho B_{rho mu nu}, B_{rho mu nu} = -B_{mu rho nu}, deg B = 3.
# Complete parity-even antisymmetric-in-(rho,mu) degree-3 basis:
#   (S2) delta_{rho nu} p_mu p^2 , (S3) delta_{mu nu} p_rho p^2 , (S5) eps_{rho mu nu sigma} p_sigma p^2
# (p_rho p_mu p_nu and delta_{rho mu} p_nu p^2 are symmetric in (rho,mu) and drop out.)
def eps4(a,b,c,d):
    # epsilon_{abcd}, a..d in 0..3
    perm=[a,b,c,d]
    if len(set(perm))<4: return 0
    # parity of permutation
    s=1
    for i in range(4):
        for j in range(i+1,4):
            if perm[i]>perm[j]: s=-s
    return s

def check(p):
    p2=sum(x*x for x in p)
    # A_{mu nu}
    A=[[-(p[mu]*p[nu]**3 - p[nu]*p[mu]**3) for nu in range(4)] for mu in range(4)]
    nonzeroA=any(A[mu][nu]!=0 for mu in range(4) for nu in range(4) if mu!=nu)
    # divergence of basis elements
    results={}
    for name in ['S2-S3','S5']:
        D=[[0]*4 for _ in range(4)]
        for mu in range(4):
            for nu in range(4):
                s=0
                for rho in range(4):
                    def B(a,b,c):
                        if name=='S2-S3':
                            return ((1 if a==c else 0)*p[b] - (1 if b==c else 0)*p[a])*p2
                        else:
                            return eps4(a,b,c,rho)*p[rho]*p2
                    # careful: eps index order (a,b,c,sigma)
                    s+=p[rho]*(B(rho,mu,nu)-B(rho,nu,mu))
                D[mu][nu]=s
        results[name]=any(D[mu][nu]!=0 for mu in range(4) for nu in range(4))
    return nonzeroA, results

# test a few generic integer momenta
for p in ([1,2,3,5],[2,-1,4,3],[1,1,2,-3],[3,4,-2,1]):
    nzA,res=check(p)
    print(p,"A nonzero:",nzA,"div(S2-S3) nonzero:",res['S2-S3'],"div(S5) nonzero:",res['S5'])
