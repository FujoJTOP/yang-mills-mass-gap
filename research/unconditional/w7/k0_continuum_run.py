"""
K0-CONTINUUM run: STREAMING, memory-safe measurement of the coarse-mode
curvature G''(0) at the softest nonzero momentum k=2*pi/L on ONE large volume
(L=8), SU(3) beta=5.7 (or SU(2) beta=2.4).

THE POINT (vs. earlier COMPUTE-GPP-SU3, which used L=4,6):
  L=4,6 reach only k=pi/2, pi/3 and L changes with k, so the L=4->6 drop mixes a
  k-decrease with a volume change.  Here we fix L=8 (one job) and scan
  n=1,2  =>  k = 2*pi/L, 4*pi/L on the SAME volume, isolating k-dependence.

MEMORY DISCIPLINE (the machine OOM-crashed on parallel SU(N) L=8):
  * ONE process, ONE lattice; configs are NEVER stored (no cfg list).
  * The sampler holds exactly one lattice (~2.4 MB for SU(3) L=8); each
    measurement's per-mode results are Nmode scalars.
  * Only N x Nmode float arrays are accumulated.  Peak RSS is logged.
  * Periodic checkpoints write the accumulated scalars (never configs).

Per configuration and mode v (flat exponential chart X=log U, a=<v,X>):
  a      = <v,x>
  H      = flat d^2/ds^2 S(X+sV)|0
  G      = flat d/ds   S(X+sV)|0
  P1,P2  = flat d/ds (-logJ), d^2/ds^2 (-logJ)|0   (analytic eigenvalue pert.)
Then exactly:  G''(0) = (E[H]-Var(G)) + W ,  W = E[P2]-Var(P1)-2Cov(G,P1).
"""
import sys
import os
import time
import ctypes
import ctypes.wintypes as wt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CPB = 'D:/Dev/math/research/lattice-cpb'
sys.path.insert(0, CPB)
sys.path.insert(0, HERE)

import w4_blocked as W          # noqa: E402
import gpp_core as C            # noqa: E402
import gpp_core3 as G3          # noqa: E402
import rerun_core as RC         # noqa: E402

# 6 transverse plane-wave modes (one polarisation, one colour each)
SPEC = [(1, 2, 3), (1, 3, 3), (0, 2, 3), (2, 3, 7), (3, 1, 3), (1, 2, 7)]      # SU(3)
SPEC2 = [(1, 2, 2), (1, 3, 2), (0, 2, 2), (2, 3, 2), (3, 1, 2), (1, 2, 0)]     # SU(2)


# ----------------------------------------------------------------------
# peak RSS (Windows), MB -- no psutil needed
# ----------------------------------------------------------------------
class _PMC(ctypes.Structure):
    _fields_ = [('cb', wt.DWORD), ('PageFaultCount', wt.DWORD),
                ('PeakWorkingSetSize', ctypes.c_size_t),
                ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t),
                ('PeakPagefileUsage', ctypes.c_size_t)]


_PSAPI = ctypes.WinDLL('psapi', use_last_error=True)
_PSAPI.GetProcessMemoryInfo.argtypes = [wt.HANDLE, ctypes.POINTER(_PMC), wt.DWORD]
_PSAPI.GetProcessMemoryInfo.restype = wt.BOOL
_K32 = ctypes.WinDLL('kernel32', use_last_error=True)
_K32.GetCurrentProcess.restype = wt.HANDLE


def _pmc():
    pmc = _PMC()
    pmc.cb = ctypes.sizeof(_PMC)
    _PSAPI.GetProcessMemoryInfo(_K32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
    return pmc


def peak_rss_mb():
    return _pmc().PeakWorkingSetSize / 1024.0 / 1024.0


def rss_mb():
    return _pmc().WorkingSetSize / 1024.0 / 1024.0


# ----------------------------------------------------------------------
# per-config measurement, S0 reused across modes (bit-identical to
# rerun_run.measure, ~25% faster), streams the config out immediately.
# ----------------------------------------------------------------------
def measure_su3(U, beta, modes, h=1e-3):
    _, x = RC.su3_log_herm(U)
    Xfrom = RC.X_from_x_su3
    action = G3.S3.wilson_action
    S0 = action(RC.expi(Xfrom(x)), beta)
    a = []; H = []; G = []; P1 = []; P2 = []
    for v in modes:
        a.append(float(np.sum(v * x)))
        Sp = action(RC.expi(Xfrom(x + h * v)), beta)
        Sm = action(RC.expi(Xfrom(x - h * v)), beta)
        G.append((Sp - Sm) / (2.0 * h))
        H.append((Sp - 2.0 * S0 + Sm) / (h * h))
        p1, p2 = RC.phi_derivs_analytic(x, v, Xfrom)
        P1.append(p1); P2.append(p2)
    return dict(a=np.array(a), H=np.array(H), G=np.array(G),
                P1=np.array(P1), P2=np.array(P2))


def measure_su2(U, beta, modes, h=1e-3):
    x = C.su2_log_x(U)
    Xfrom = RC.X_from_x_su2
    action = W.action2
    S0 = action(RC.expi(Xfrom(x)), beta)
    a = []; H = []; G = []; P1 = []; P2 = []
    for v in modes:
        a.append(float(np.sum(v * x)))
        Sp = action(RC.expi(Xfrom(x + h * v)), beta)
        Sm = action(RC.expi(Xfrom(x - h * v)), beta)
        G.append((Sp - Sm) / (2.0 * h))
        H.append((Sp - 2.0 * S0 + Sm) / (h * h))
        p1, p2 = RC.phi_derivs_analytic(x, v, Xfrom)
        P1.append(p1); P2.append(p2)
    return dict(a=np.array(a), H=np.array(H), G=np.array(G),
                P1=np.array(P1), P2=np.array(P2))


def build_modes(group, L, nmax):
    modes = []
    names = []
    kvals = []
    spec = SPEC if group == 'su3' else SPEC2
    for n in range(1, nmax + 1):
        for (ax, pl, c) in spec:
            if group == 'su3':
                v = G3.plane_mode(L, ax, pl, c, n=n, kind='cos')
            else:
                v = C.plane_mode(L, ax, pl, c, n=n, kind='cos')
            modes.append(v)
            names.append(f'ax{ax}_p{pl}_c{c}_n{n}')
            kvals.append(2.0 * np.pi * n / L)
    return modes, names, np.array(kvals)


def main(group, L, beta, N, gap, therm, seed, nmax, out):
    stream = (measure_su3 if group == 'su3' else measure_su2)
    modes, names, kvals = build_modes(group, L, nmax)
    Nm = len(modes)
    print(f"# K0-CONTINUUM {group} L={L} beta={beta} N={N} gap={gap} therm={therm} "
          f"seed={seed} nmax={nmax} modes={Nm}", flush=True)
    print(f"  k = {sorted(set(np.round(kvals,5)))}  (2pi/L={2*np.pi/L:.5f})", flush=True)

    if group == 'su3':
        import su3 as S3
        rng = np.random.default_rng(seed)
        lat = S3.Lattice(L, beta, rng, cold=True)
        sweep = lat.sweep
        getU = lambda: lat.U
        import su3 as _s3
        plaq_fn = lambda: _s3.mean_plaquette(lat.U)
    else:
        rng = np.random.default_rng(seed)
        lat_U = [W.random_su2(rng, (4, L, L, L, L))]   # single-element holder

        def sweep():
            U = lat_U[0]
            Uf = U.reshape(4, -1, 2, 2)
            ns = L ** 4
            idx = np.arange(ns).reshape(L, L, L, L)
            x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
            par = (x0 + x1 + x2 + x3) % 2
            idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
            kappa = beta / 2.0
            for c in (0, 1):
                ids = idxp[c]
                for mu in range(4):
                    S = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                    A = Uf[mu, ids]
                    M = A @ S
                    v = __import__('su3').su2_heatbath(kappa, M, 0, 1, rng)
                    Uf[mu, ids] = v @ A
        getU = lambda: lat_U[0]
        # action2 = beta*(6L^4 - sum 1/2 Re tr U_p); true plaquette = 1 - that ratio
        plaq_fn = lambda: 1.0 - W.action2(lat_U[0], beta) / (beta * 6.0 * L ** 4)

    t0 = time.time()
    for _ in range(therm):
        sweep()
    print(f"  therm {therm} sweeps in {time.time()-t0:.0f}s; "
          f"plaq={plaq_fn():.4f}; RSS={rss_mb():.0f}MB peak={peak_rss_mb():.0f}MB",
          flush=True)

    acc = {k: np.empty((N, Nm)) for k in ('a', 'H', 'G', 'P1', 'P2')}
    t0 = time.time()
    tlast = t0
    for i in range(N):
        for _ in range(gap):
            sweep()
        r = stream(getU(), beta, modes)
        for k in acc:
            acc[k][i] = r[k]
        if (i + 1) % 50 == 0 or i == 0:
            dt = time.time() - tlast
            tlast = time.time()
            eta = (N - i - 1) / (i + 1) * (time.time() - t0)
            print(f"    {i+1}/{N}  ({time.time()-t0:.0f}s, ETA {eta:.0f}s)  "
                  f"RSS={rss_mb():.0f} peak={peak_rss_mb():.0f}MB", flush=True)
        if (i + 1) % 100 == 0:                      # checkpoint (scalars only)
            np.savez(out, **acc, names=np.array(names), kvals=kvals,
                     L=L, beta=beta, group=group, N_done=i + 1)
            print(f"    [chk {i+1}] plaq={plaq_fn():.4f}", flush=True)

    np.savez(out, **acc, names=np.array(names), kvals=kvals,
             L=L, beta=beta, group=group, N_done=N)
    print(f"  DONE saved {out}  total {time.time()-t0:.0f}s  "
          f"peak RSS={peak_rss_mb():.0f}MB", flush=True)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--group', default='su3', choices=['su3', 'su2'])
    ap.add_argument('--L', type=int, default=8)
    ap.add_argument('--beta', type=float, default=None)
    ap.add_argument('--N', type=int, default=400)
    ap.add_argument('--gap', type=int, default=4)
    ap.add_argument('--therm', type=int, default=350)
    ap.add_argument('--seed', type=int, default=808)
    ap.add_argument('--nmax', type=int, default=2)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    beta = a.beta if a.beta is not None else (5.7 if a.group == 'su3' else 2.4)
    out = a.out or f'k0_{a.group}_b{beta}_L{a.L}_n{a.nmax}.npz'
    main(a.group, a.L, beta, a.N, a.gap, a.therm, a.seed, a.nmax, out)
