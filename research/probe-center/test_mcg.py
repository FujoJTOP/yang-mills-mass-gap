import numpy as np, time, sys
sys.path.insert(0,'D:/Dev/math/research/probe-center')
import mcg
np.set_printoptions(precision=4,suppress=True,linewidth=160)

def run(fn, nsw=22):
    d=np.load(fn, allow_pickle=True)
    U=d['snapshots'][0].copy(); L=int(d['L']); beta=float(d['beta'])
    g=np.zeros((L,L,L,L,3,3),np.complex128)
    for i in range(3): g[...,i,i]=1.0
    Rs=[]
    for sw in range(nsw): mcg.mcg_sweep(U,g,L,6)
    V=mcg.apply_gauge(U,g,L)
    R=mcg.center_functional(V,L)
    n,amb=mcg.center_project(V,L)
    flux=mcg.plaquette_flux(n,L)
    plaqs,labels,sizes=mcg.vortex_clusters(flux,L)
    nv=len(plaqs)
    Wf=mcg.wilson_loop_tensor(U,L,3)
    Uz=mcg.expand_center_links(n,L)
    Wz=mcg.wilson_loop_tensor(Uz,L,3)
    # adjoint full: (|Tr|^2-1)/8 ; here Wf uses Re Tr/3 so |Tr|^2 not directly
    # compute adjoint separately below
    print(f"\n=== {fn.split('/')[-1]} L={L} beta={beta} ===")
    print(f"R={R:.4f}  amb={amb:.4f}  nvortex={nv}  density/plaq={nv/(6*L**4):.4f}  nclust={len(sizes)}")
    if len(sizes): print(f"  cluster sizes: max={sizes.max()} mean={sizes.mean():.1f} n>1={np.sum(sizes>1)}")
    k1=sum(1 for p in plaqs if p[3]==1); k2=sum(1 for p in plaqs if p[3]==2)
    print(f"  k=1:{k1} k=2:{k2}")
    print("  W_full(1,1)=%.4f  W_full(1,2)=%.4f  W_full(2,2)=%.4f"%(Wf[1,1],Wf[1,2],Wf[2,2]))
    print("  W_proj(1,1)=%.4f  W_proj(1,2)=%.4f  W_proj(2,2)=%.4f"%(Wz[1,1],Wz[1,2],Wz[2,2]))

run('D:/Dev/math/research/lattice-cpb/runs/L6_b2.600_s51.npz')
run('D:/Dev/math/research/lattice-cpb/runs/L6_b6.000_s32.npz')
run('D:/Dev/math/research/lattice-cpb/runs/L8_b6.000_s33.npz')
