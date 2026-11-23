import glob,os,sys,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE,'..','lattice-cpb'))
from probe import adj, eig_angles_batch
files=[f for f in sorted(glob.glob(os.path.join(HERE,'..','lattice-cpb','runs','L*_b*.npz')))
       if os.path.basename(f).startswith(('L6_b2.6','L6_b2.8','L6_b3.0','L6_b3.3','L8_b2.2','L6_b2.2','L6_b6.0'))
       and 'block' not in f and 's24' not in f]
print(f'{"file":24s} {"Nsite":>5s} {"mean":>8s} {"min":>7s} {"max":>6s} {"Eexp|V|/2":>11s} {"Eexp|V|":>9s} {"Eexp(2|V|)":>10s}')
for f in files:
    d=np.load(f); sn=d['snapshots']; beta=float(d['beta']); L=int(d['L'])
    acc=[]
    for c in range(sn.shape[0]):
        U=sn[c]; ps=[];qs=[]
        for mu in range(4):
            for nu in range(mu+1,4):
                Up=(U[mu]@np.roll(U[nu],-1,axis=mu)@adj(np.roll(U[mu],-1,axis=nu))@adj(U[nu]))
                F=Up.reshape(-1,3,3); ang=eig_angles_batch(F)
                qs.append(np.sum(ang*ang,axis=-1)/6.0)
                ps.append(1.0-np.real(np.trace(F,axis1=-2,axis2=-1))/3.0)
        v=beta*(np.array(ps)-np.array(qs)).sum(axis=0).reshape(-1)
        acc.append(v)
    v=np.concatenate(acc); av=np.abs(v)
    print(f'{os.path.basename(f):24s} {v.size:5d} {v.mean():8.3f} {v.min():7.2f} {v.max():6.2f} '
          f'{np.mean(np.exp(av/2)):11.2f} {np.mean(np.exp(av)):9.2f} {np.mean(np.exp(2*av)):10.2f}')
