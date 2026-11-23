"""Warm-started L=8 configs at beta=2.6 then 3.0, from the equilibrated
beta=2.2 L8 config (P~0.143).  Long thermalisation, 40 snapshots each."""
import numpy as np, time, sys
sys.path.insert(0,'D:/Dev/math/research/lattice-cpb')
import su3_numba as SN

src = np.load('D:/Dev/math/research/lattice-cpb/runs/L8_b2.200_s23.npz', allow_pickle=True)
U = src['snapshots'][-1].copy(); L=8
print("warm start P =", SN.plaquette_mean(U,L), flush=True)

def run(beta, therm, sep, nsnap):
    kappa = beta/3.0; t0=time.time()
    for k in range(therm):
        SN.sweep_par(U, L, kappa)
        if (k+1)%250==0: print(f"  b={beta} therm {k+1}/{therm} P={SN.plaquette_mean(U,L):.5f} t={time.time()-t0:.0f}s",flush=True)
    snaps=[]
    for n in range(nsnap):
        for s in range(sep): SN.sweep_par(U,L,kappa)
        snaps.append(U.copy())
        if (n+1)%10==0: print(f"  b={beta} snap {n+1}/{nsnap} P={SN.plaquette_mean(U,L):.5f} t={time.time()-t0:.0f}s",flush=True)
    np.savez(f"runs_L8_b{beta:.3f}.npz", L=L, beta=beta, snapshots=np.array(snaps), therm=therm, sep=sep)
    print("saved", f"runs_L8_b{beta:.3f}.npz", "P=", SN.plaquette_mean(U,L), f"t={time.time()-t0:.0f}s",flush=True)

run(2.6, 2500, 20, 40)
run(3.0, 2500, 20, 40)
