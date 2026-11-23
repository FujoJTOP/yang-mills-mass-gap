import sys, time, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lattice-cpb'))
import numpy as np
import su3_numba as N

L=6; beta=2.6; kappa=beta/3.0
U=np.zeros((4,L,L,L,L,3,3),dtype=np.complex128)
for mu in range(4): U[mu][...]=np.eye(3)
np.random.seed(1)
t0=time.time()
for k in range(1000): N.sweep_par(U,L,kappa)
print('therm 1000 sweeps', round(time.time()-t0,2),'s  P=',N.plaquette_mean(U,L))

def adj(A): return np.conj(np.swapaxes(A,-1,-2))
def plaq_arrays(U,b=1.0):
    beta=b
    ps=[];qs=[]
    for mu in range(4):
        for nu in range(mu+1,4):
            Up=(U[mu]@np.roll(U[nu],-1,axis=mu)@adj(np.roll(U[mu],-1,axis=nu))@adj(U[nu]))
            F=Up.reshape(-1,3,3)
            ev=np.linalg.eigvals(F)
            ang=np.angle(ev)
            q=np.sum(ang*ang,axis=-1)/6.0
            p=1.0-np.real(np.trace(F,axis1=-2,axis2=-1))/3.0
            ps.append(p); qs.append(q)
    return np.array(ps),np.array(qs)

t0=time.time(); p,q=plaq_arrays(U); print('eig 1 config',round(time.time()-t0,3),'s')
for alpha in (0.0,0.5,1.0,1.5):
    V=beta*np.sum(p-alpha*q,axis=0)
    d=V-V.mean()
    print(f'alpha={alpha}: meanV={V.mean():.3f} var={d.var():.4f} k3={np.mean(d**3):.4f} skew={np.mean(d**3)/d.var()**1.5:.3f} k4R={ (np.mean(d**4)-3*d.var()**2)/d.var()**2:.3f}')
# field control
Pf=np.mean(np.real(np.trace(np.stack([ (U[mu]@np.roll(U[nu],-1,axis=mu)@adj(np.roll(U[mu],-1,axis=nu))@adj(U[nu])).reshape(-1,3,3) for mu in range(4) for nu in range(mu+1,4)]),axis1=-2,axis2=-1))/3.0,axis=0)
print('field P: var',Pf.var(),'k3',np.mean((Pf-Pf.mean())**3))
