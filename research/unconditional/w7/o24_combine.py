"""O24 combine: Delta(L=2) with sigma_jk (grouped jackknife) AND sigma_tau (autocorrelation),
plus heavy-tail diagnostics (mean/median/winsor/Hill), raw and L^2-normalised verdicts."""
import numpy as np
import sys

TAG = sys.argv[1] if len(sys.argv) > 1 else 'L2A'
d = np.load(f'o24_{TAG}.npz')
K, F = d['K'], d['F']
G = K.shape[1]
N = K.shape[0]
L = 2
NB = 40                       # jackknife blocks per chain


def acf(x, ml=600):
    y = x - x.mean()
    v = (y * y).sum()
    return np.array([(y[:-k] * y[k:]).sum() / v for k in range(1, min(ml, len(y) - 2))])


def geyer(x, ml=600):
    r = np.concatenate(([1.0], acf(x, ml)))
    if len(r) % 2:
        r = r[:-1]
    s = np.maximum(r[0::2] + r[1::2], 0.0)
    tot = 0.0
    for si in s:
        if si <= 0:
            break
        tot += si
    return max(2 * tot - 1.0, 1.0)


def dl(k, f):
    return k.mean() - f.var()


def group_stats(k, f, nb=NB):
    """Delta, block-jackknife sigma, tau-corrected sigma, Neff."""
    n = len(k)
    idx = np.array_split(np.arange(n), nb)
    ds = np.array([dl(np.delete(k, b), np.delete(f, b)) for b in idx])
    s_jk = np.sqrt((nb - 1) / nb * ((ds - ds.mean()) ** 2).sum())
    tK, tF = geyer(k), geyer(f)
    vK, vF = k.var(ddof=1), f.var(ddof=1)
    sK = np.sqrt(vK * (1 + 2 * tK) / n)
    sF = np.sqrt(2 * vF ** 2 * (1 + 2 * tF) / n)
    s_tau = np.sqrt(sK ** 2 + sF ** 2)
    return dict(d=dl(k, f), s_jk=s_jk, s_tau=s_tau, tK=tK, tF=tF,
                neff=n / max(tK, tF), mK=k.mean(), vF=vF, mF=f.mean())


def hill(x, frac=0.01):
    y = np.sort(np.abs(x - np.median(x)))[::-1]
    k = max(int(frac * len(y)), 5)
    y = y[:k]
    if y[-1] <= 0:
        return np.nan
    return 1.0 / np.mean(np.log(y / y[-1]))


def wins(x, q):
    lo, hi = np.quantile(x, [q, 1 - q])
    return np.clip(x, lo, hi)


print(f"=== O24  tag={TAG} mode={d['mode']} Nsweep/chain={N} chains={G} seeds={d['seeds']} "
      f"eps={d['eps']} domain-reject={float(d['domrej']):.4f} burn=800 wall={float(d['el']):.0f}s "
      f"({1e3*float(d['el'])/(N+800):.2f} ms/sweep for all {G} chains) ===")
print(f"acc_link={np.round(d['acc'],3)}")

rows = []
for g in range(G):
    r = group_stats(K[:, g], F[:, g])
    rows.append(r)
    print(f"g={g} seed={int(d['seeds'][g])} E[k]= {r['mK']:8.4f} Var(F)={r['vF']:7.4f} E[F]={r['mF']:+.4f} "
          f"D={r['d']:+.4f} s_jk={r['s_jk']:.4f} s_tau={r['s_tau']:.4f} "
          f"tauK={r['tK']:5.1f} tauF={r['tF']:5.1f} Neff={r['neff']:6.0f}")

ds = np.array([r['d'] for r in rows])
sjk = np.array([r['s_jk'] for r in rows])
stau = np.array([r['s_tau'] for r in rows])
D = ds.mean()
S_jk = np.sqrt((sjk ** 2).sum()) / G
S_tau = np.sqrt((stau ** 2).sum()) / G
# delete-one-group jackknife over groups (catches between-group scatter)
dg = G * D - (G - 1) * ds
S_grp = np.sqrt((G - 1) / G * ((dg - dg.mean()) ** 2).sum())
print(f"\n--- 4-chain combination ---")
print(f"Delta_raw = {D:+.4f}   sigma_jk={S_jk:.4f}   sigma_tau={S_tau:.4f}   "
      f"sigma_grp={S_grp:.4f}   spread={ds.std(ddof=1):.4f}")
print(f"Delta_raw/sigma_jk = {D/S_jk:+.2f}   /sigma_tau = {D/S_tau:+.2f}   /sigma_grp = {D/S_grp:+.2f}")
print(f"chi2/dof over 4 groups = {np.sum((ds-D)**2/sjk**2)/(G-1):.2f}")

kall, fall = K.ravel(), F.ravel()
print(f"\n--- heavy-tail diagnostics (pooled {kall.size} samples) ---")
med, mad = np.median(kall), np.median(np.abs(kall - np.median(kall)))
print(f"kappa: mean={kall.mean():.4f} median={med:.4f} MAD={mad:.4f} max={kall.max():.3e} "
      f"min={kall.min():.3e} std={kall.std():.4f}  max/median={kall.max()/max(med,1e-9):.1f}")
print(f"F    : mean={fall.mean():+.4f} median={np.median(fall):+.4f} max|F|={np.abs(fall).max():.3e} std={fall.std():.4f}")
print(f"Hill alpha(kappa) top1% = {hill(kall):.2f}   Hill alpha(F) = {hill(fall):.2f}")
for q in (0.001, 0.01, 0.05):
    print(f"  winsor {100*q:.1f}%: Delta={dl(wins(kall,q), wins(fall,q)):+.4f}")
print(f"  median-kappa version: Delta={med - fall.var():+.4f}")
print(f"  frac |kappa-med|>10*MAD: {np.mean(np.abs(kall-med)>10*mad):.4f}")
# share of E[kappa] carried by the top 0.1% of samples
top = np.sort(kall)[::-1]
for frac in (0.0001, 0.001):
    m = max(int(frac * kall.size), 1)
    print(f"  top {100*frac:.2f}% of kappa carries {100*top[:m].sum()/kall.sum():.1f}% of E[kappa]")

print(f"\n--- normalisation (ch05:208/:214-228: d_v unnormalised, ||v||^2 = 1/L) ---")
for nm, fac in (('Delta_raw ', 1.0), ('Delta_norm=L^2*raw', float(L ** 2))):
    print(f"{nm}: {fac*D:+.4f} +- {fac*S_jk:.4f} (jk) / {fac*S_tau:.4f} (tau)   "
          f"ratio_jk={D/S_jk:+.2f} ratio_tau={D/S_tau:+.2f}")

print(f"\n--- needed sample for 5 sigma ---")
for nm, s in (('jk', S_jk), ('tau', S_tau)):
    need = N * (s / (abs(D) / 5)) ** 2
    print(f"  using sigma_{nm}={s:.4f}: sigma<= {abs(D)/5:.4f} needs x{(s/(abs(D)/5))**2:.2f} "
          f"=> {need:.0f} sweeps/chain (x4 chains), wall ~ {need*1e3*float(d['el'])/(N+800)/G/60:.1f} min")

