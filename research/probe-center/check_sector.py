import numpy as np, sys
sys.path.insert(0,'D:/Dev/math/research/probe-center')
import mcg
d=np.load('D:/Dev/math/research/lattice-cpb/runs/L6_b6.000_s32.npz',allow_pickle=True)
U=d['snapshots'][0].copy(); L=6
# use existing probe output? just gauge fix quickly
g=np.zeros((L,L,L,L,3,3),np.complex128)
for i in range(3): g[...,i,i]=1.0
for s in range(40): mcg.mcg_sweep(U,g,L,6)
V=mcg.apply_gauge(U,g,L)
n,amb=mcg.center_project(V,L)
Uz=mcg.expand_center_links(n,L)
Tf,Tf2=mcg.wilson_loop_stats(V,L,3)      # gauge-invariant anyway
Tp,Tp2=mcg.wilson_loop_stats(Uz,L,3)
WF=Tf.real; WA=(9*Tf2-1)/8
WFp=Tp.real; WAp=(9*Tp2-1)/8
np.set_printoptions(precision=4,suppress=True,linewidth=140)
print("full fundamental W_F:\n",WF[1:,1:])
print("full adjoint   W_A:\n",WA[1:,1:])
print("proj fundamental W_F:\n",WFp[1:,1:])
print("proj adjoint max|W_A-1| =", np.abs(WAp[1:,1:]-1).max())
print("proj |z|^2 (should be 1):", Tp2[1:,1:])
