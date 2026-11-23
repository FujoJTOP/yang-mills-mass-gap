import glob,os,sys,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE,'..','lattice-cpb'))
from probe import adj, eig_angles_batch
files=[f for f in sorted(glob.glob(os.path.join(HERE,'..','lattice-cpb','runs','L*_b*.npz')))
       if os.path.basename(f).startswith(('L6_b2.6','L6_b2.8','L6_b3.0','L6_b3.3','L8_b2.2','L6_b2.2','L6_b6.0','L8_b6.0'))
       and 'block' not in f and 's24' not in f]
print(f'{"file":24s} {"a":>4s} {"mean":>8s} {"maxV":>7s} {"max|V|":>7s} {"Eexp|V|/2":>11s} {"Eexp|V|":>11s} {"Eexp(|V|)":>10s}')
for f in files:
    d=np.load(f); sn=d['snapshots']; beta=float(d['beta']); N=sn.shape[0]
    vals={1.0:[],0.5:[]}
    for c in range(N):
        U=sn[c]; ps=[];qs=[]
        for mu in range(4):
            for nu in range(mu+1,4):
                Up=(U[mu]@np.roll(U[nu],-1,axis=mu)@adj(np.roll(U[mu],-1,axis=nu))@adj(U[nu]))
                F=Up.reshape(-1,3,3); ang=eig_angles_batch(F)
                qs.append(np.sum(ang*ang,axis=-1)/6.0)
                ps.append(1.0-np.real(np.trace(F,axis1=-2,axis2=-1))/3.0)
        ps=np.array(ps).reshape(-1); qs=np.array(qs).reshape(-1)
        for a in vals: vals[a].append(beta*(ps-a*qs))
    for a in (1.0,0.5):
        v=np.concatenate(vals[a]); av=np.abs(v)
        print(f'{os.path.basename(f):24s} {a:4.1f} {v.mean():8.3f} {v.max():7.2f} {av.max():7.2f} '
              f'{np.mean(np.exp(av/2)):11.3f} {np.mean(np.exp(av)):11.2f} {np.mean(np.exp(av)):10.2f}')
