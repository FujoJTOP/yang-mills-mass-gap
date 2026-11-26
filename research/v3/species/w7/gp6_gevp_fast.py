"""GP6 -- fast observable-basis GEVP on the Lt=20 sample. O(N) pre-aggregation, no per-block recompute."""
import numpy as np, sys, time
fn = sys.argv[1] if len(sys.argv) > 1 else 'gh3_probe.npz'
OBS = [0, 3, 4]            # W11_p0, W12_p0, W22_p0
NAM = ['W11', 'W12', 'W22']
NB  = 20
t0s = time.time()
d = np.load(fn)
S, Lt, N = d['S'], int(d['Lt']), int(d['n'])
print(f"loaded {fn}  N={N} Lt={Lt}  {time.time()-t0s:.1f}s", flush=True)
A = np.ascontiguousarray(S[:, 0][:, OBS, :]).astype(np.float64)   # (N,3,Lt)
del S
t0s = time.time()
Cs = np.empty((Lt, N, 3, 3))
for t in range(Lt):
    Cs[t] = np.einsum('nct,ndt->ncd', A, np.roll(A, -t, axis=2), optimize=True) / Lt
print(f"correlators {time.time()-t0s:.1f}s", flush=True)
blocks = np.array_split(np.arange(N), NB)
Pb = np.stack([Cs[:, b].sum(1) for b in blocks], 0)      # (NB,Lt,3,3)
tot = Cs.sum(1)
del Cs
def gevp(Ct, t0=0):
    C0 = Ct[t0]
    try:
        L = np.linalg.cholesky(C0)
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh(C0); k = w > 1e-14 * max(w.max(), 1e-300)
        L = V[:, k] * np.sqrt(w[k])
    P = np.linalg.solve(L, np.eye(len(L)))
    out = []
    for t in range(Lt):
        M = P @ Ct[t] @ P.T
        out.append(np.linalg.eigvalsh(0.5 * (M + M.T))[-1])
    return np.array(out)
lt = np.array([gevp(tot)[t] for t in range(Lt)])
jk = np.array([[gevp(tot - Pb[b])[t] for t in range(Lt)] for b in range(NB)])
err = np.sqrt((NB - 1) / NB * ((jk - jk.mean(0)) ** 2).sum(0))
print("\n  t     lam_max(t)      m_eff(t)    +- JK      S/N")
for t in range(1, min(10, Lt - 1)):
    me = np.log(lt[t] / lt[t + 1])
    mek = np.array([np.log(j[b, t] / j[b, t + 1]) for b in range(NB) for j in [jk]])
    e = np.sqrt((NB - 1) / NB * ((mek - mek.mean()) ** 2).sum())
    print(f" {t:2d}  {lt[t]:12.6f}  {me:10.4f}  {e:8.4f}  {me/e if e>0 else np.nan:7.2f}")
print(f"\ntotal {time.time()-t0s:.1f}s")
