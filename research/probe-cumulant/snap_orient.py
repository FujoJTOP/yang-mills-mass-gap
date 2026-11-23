import glob,os,sys,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE,'..','lattice-cpb'))
from probe import adj, eig_angles_batch
f=os.path.join(HERE,'..','lattice-cpb','runs','L6_b3.000_s53.npz')
d=np.load(f); sn=d['snapshots']; beta=3.0
acc={}
for c in range(sn.shape[0]):
    U=sn[c]
    for mu in range(4):
        for nu in range(mu+1,4):
            Up=(U[mu]@np.roll(U[nu],-1,axis=mu)@adj(np.roll(U[mu],-1,axis=nu))@adj(U[nu]))
            F=Up.reshape(-1,3,3); ang=eig_angles_batch(F)
            q=np.sum(ang*ang,axis=-1)/6.0; p=1.0-np.real(np.trace(F,axis1=-2,axis2=-1))/3.0
            v=beta*(p-q)
            a=acc.setdefault((mu,nu),[])
            a.append(v)
print('beta=3.0 L=6: per-orientation k3(per plaquette), mean, var')
for k in sorted(acc):
    v=np.concatenate(acc[k]); m=v.mean(); x=v-m
    print(k, 'k3=%.4f var=%.4f mean=%.4f'%(np.mean(x**3),x.var(),m))
