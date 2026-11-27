"""GP6b -- connected (vacuum-subtracted) GEVP, O(N) pre-aggregation."""
import numpy as np, sys, time
fn = sys.argv[1] if len(sys.argv) > 1 else 'gh3_probe.npz'
OBS = [0, 3, 4]; NAM = ['W11','W12','W22']; NB = 20
t0s = time.time()
d = np.load(fn); S, Lt, N = d['S'], int(d['Lt']), int(d['n'])
A = np.ascontiguousarray(S[:, 0][:, OBS, :]).astype(np.float64); del S
Abar = A.mean(axis=(0, 2))                      # (3,) vacuum expectation
print("vacuum means:", " ".join(f"{NAM[i]}={Abar[i]:.6f}" for i in range(3)), flush=True)
Ac = A - Abar[None, :, None]                    # connected operator
t1 = time.time()
Cs = np.empty((Lt, N, 3, 3))
for t in range(Lt):
    Cs[t] = np.einsum('nct,ndt->ncd', Ac, np.roll(Ac, -t, axis=2), optimize=True) / Lt
print(f"connected correlators {time.time()-t1:.1f}s", flush=True)
blocks = np.array_split(np.arange(N), NB)
Pb = np.stack([Cs[:, b].sum(1) for b in blocks], 0)
tot = Cs.sum(1); del Cs
print("\n  diag ratio C(t)/C(0):  " + "  ".join(f"t={t}" for t in range(1,7)))
for i in range(3):
    print(f"  {NAM[i]}: " + "  ".join(f"{tot[t,i,i]/tot[0,i,i]:.5f}" for t in range(1,7)))
def gevp(Ct, t0=0):
    C0 = Ct[t0]
    w, V = np.linalg.eigh(0.5*(C0+C0.T))
    k = w > 1e-12*max(w.max(), 1e-300)
    if k.sum() < len(w): print(f"   [note] C(0) rank {k.sum()}/{len(w)}", flush=True)
    L = V[:, k] * np.sqrt(w[k])
    P = np.linalg.solve(L, np.eye(L.shape[1]))
    return np.array([np.linalg.eigvalsh(0.5*((P@Ct[t]@P.T)+(P@Ct[t]@P.T).T))[-1] for t in range(Lt)])
lt = gevp(tot)
jk = np.array([gevp(tot - Pb[b]) for b in range(NB)])
print("\n  t     lam_max(t)      m_eff(t)    +- JK     S/N")
for t in range(1, min(10, Lt-1)):
    me = np.log(lt[t]/lt[t+1])
    v = np.array([np.log(j[t]/j[t+1]) for j in jk])
    e = np.sqrt((NB-1)/NB*((v-v.mean())**2).sum())
    print(f" {t:2d}  {lt[t]:12.6e}  {me:10.4f}  {e:8.4f}  {me/e if e>0 else np.nan:7.2f}")
print(f"\ntotal {time.time()-t0s:.1f}s")
