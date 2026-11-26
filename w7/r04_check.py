import itertools, numpy as np
from scipy.integrate import quad
L,D=2,4
sites=list(itertools.product(range(L),repeat=D)); sid={s:i for i,s in enumerate(sites)}
links=[(x,mu) for x in sites for mu in range(D)]
rng=np.random.default_rng(1)
def rot(q):  # adjoint rotation of su(2) from unit quaternion q=(w,x,y,z)
    w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
for trial in range(3):
    J=np.zeros((3*len(links),3*len(sites)))
    for l,(pos,mu) in enumerate(links):
        posp=list(pos); posp[mu]=(pos[mu]+1)%L
        R=rot(rng.normal(size=4)); R/=np.linalg.norm(R,axis=0)[0]
        for b in range(3):
            J[3*l+b,3*sid[pos]+b]+=1.0
            J[3*l:3*l+3,3*sid[tuple(posp)]:3*sid[tuple(posp)]+3]-=R
    print("generic config gauge-Jacobian rank =",np.linalg.matrix_rank(J,tol=1e-9),"/48  => orbit dim a.e. 48 => physical dim 192-48 = 144")
print()
for n in (2,3):
    f=lambda e,n=n: quad(lambda t: (lambda d1: d1**2*(e-d1)**2*(1 if n==2 else 1))(t), 0, e)[0] if n==2 else \
        quad(lambda t: t**2*quad(lambda s: s**2*(e-t-s)**2,0,e-t)[0],0,e)[0]
    es=np.array([1e-1,1e-2,1e-3]); vals=np.array([f(e) for e in es])
    print(f"n={n}: slope d log F / d log eps = {np.polyfit(np.log(es),np.log(vals),1)[0]:.10f}  (expected 3n-1={3*n-1})")
