import numpy as np, sys, time
sys.path.insert(0,'D:/Dev/math/research/probe-center')
import mcg
d=np.load('D:/Dev/math/research/lattice-cpb/runs/L6_b6.000_s32.npz',allow_pickle=True)
U=d['snapshots'][0].copy(); L=6
g=np.zeros((L,L,L,L,3,3),np.complex128)
for i in range(3): g[...,i,i]=1.0
t0=time.time()
for sw in range(60):
    mcg.mcg_sweep(U,g,L,6)
    if sw%5==0:
        V=mcg.apply_gauge(U,g,L); R=mcg.center_functional(V,L)
        n,amb=mcg.center_project(V,L); flux=mcg.plaquette_flux(n,L)
        p,l,s=mcg.vortex_clusters(flux,L)
        print(f"sw {sw:2d} R={R:.4f} dens={len(p)/(6*L**4):.4f} nclust={len(s)} t={time.time()-t0:.1f}s",flush=True)
# random-gauge baseline density: project raw U directly
n0,_=mcg.center_project(U,L); f0=mcg.plaquette_flux(n0,L); p0,_,_=mcg.vortex_clusters(f0,L)
print("RAW (unfixed) density:", len(p0)/(6*L**4))
