"""Generate L=8 SU(3) configs at beta=2.6, 3.0 (primary-grid couplings).
Numba parallel heat bath; snapshots separated by `sep` sweeps."""
import numpy as np, time, sys
sys.path.insert(0,'D:/Dev/math/research/lattice-cpb')
import su3_numba as SN

def gen(L, beta, seed, therm, sep, nsnap):
    rng = np.random.default_rng(seed)
    U = np.tile(np.eye(3, dtype=np.complex128), (4, L, L, L, L, 1, 1))
    # cold start; small randomisation
    kappa = beta/6.0
    t0=time.time()
    for k in range(therm):
        SN.sweep_par(U, L, kappa)
        if (k+1) % 250 == 0:
            P = SN.plaquette_mean(U, L)
            print(f"  therm {k+1}/{therm} P={P:.5f} t={time.time()-t0:.0f}s", flush=True)
    snaps=[]
    for n in range(nsnap):
        for s in range(sep):
            SN.sweep_par(U, L, kappa)
        snaps.append(U.copy())
        if (n+1) % 10 == 0:
            print(f"  snap {n+1}/{nsnap} P={SN.plaquette_mean(U,L):.5f} t={time.time()-t0:.0f}s", flush=True)
    out=f"runs_L8_b{beta:.3f}.npz"
    np.savez(out, L=L, beta=beta, seed=seed, snapshots=np.array(snaps),
             therm=therm, sep=sep)
    print("saved", out, "nsnap", len(snaps), f"total {time.time()-t0:.0f}s", flush=True)

if __name__=='__main__':
    gen(8, 2.6, 77, 1500, 20, 40)
    gen(8, 3.0, 78, 1500, 20, 40)
