"""Independent cross-check on the saved lattice-cpb snapshots (40 configs per
run, separated by 400 sweeps).  Per-site and per-plaquette k3 of
V(alpha)=beta(p-alpha q) at alpha=0,0.5,1 and the field control."""
import glob, os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE,'..','lattice-cpb'))
from probe import adj, eig_angles_batch

def moments(v):
    return np.array([v.sum(),(v*v).sum(),(v**3).sum(),(v**4).sum()])

def kfrom(s,n):
    m1,m2,m3,m4=s/n
    k3=m3-3*m1*m2+2*m1**3
    k2=m2-m1*m1
    return m1,k2,k3

files=[f for f in sorted(glob.glob(os.path.join(HERE,'..','lattice-cpb','runs','L*_b*.npz')))
       if os.path.basename(f).startswith(('L6_b2.6','L6_b2.8','L6_b3.0','L6_b3.3','L8_b2.2','L6_b2.2','L6_b6.0','L8_b6.0'))
       and 'block' not in f and 's24' not in f]
print(f'{"file":26s} {"N":>3s}  alpha   per-site k3     per-plaquette k3    var_site  var_plaq  k3(field)')
for f in files:
    d=np.load(f); sn=d['snapshots']; beta=float(d['beta']); L=int(d['L'])
    N=sn.shape[0]
    acc={a:[np.zeros(4),np.zeros(4)] for a in (0.0,0.5,1.0)}
    fld=np.zeros(4)
    nsite=0; npl=0
    for c in range(N):
        U=sn[c]
        ps=[];qs=[]
        for mu in range(4):
            for nu in range(mu+1,4):
                Up=(U[mu]@np.roll(U[nu],-1,axis=mu)@adj(np.roll(U[mu],-1,axis=nu))@adj(U[nu]))
                F=Up.reshape(-1,3,3); ang=eig_angles_batch(F)
                qs.append(np.sum(ang*ang,axis=-1)/6.0)
                ps.append(1.0-np.real(np.trace(F,axis1=-2,axis2=-1))/3.0)
        ps=np.array(ps);qs=np.array(qs)
        fld+=moments(ps.reshape(-1)); npl+=ps.size
        for a in acc:
            v=beta*(ps-a*qs)
            acc[a][1]+=moments(v.reshape(-1))
            acc[a][0]+=moments(v.sum(axis=0))
        nsite+=L**4
    out=f'{os.path.basename(f):26s} {N:3d}'
    for a in (0.0,0.5,1.0):
        _,k2s,k3s=kfrom(acc[a][0],nsite)
        _,k2p,k3p=kfrom(acc[a][1],npl)
        _,_,k3f=kfrom(fld,npl)
        print(f'{out if a==0.0 else "":26s}      {a:4.1f}  {k3s:10.4f}      {k3p:10.5f}       {k2s:7.3f} {k2p:8.4f}  {k3f:+.2e}')
