import numpy as np
tot = []
for m, s in (("A", 11), ("A", 12), ("B", 21)):
    K, F = np.load(f'o11_lat_{m}_{s}.npy')
    n = len(K); idx = np.array_split(np.arange(n), 20)

    def dl(sel):
        return K[sel].mean() - F[sel].var()

    d0 = dl(slice(None))
    ds = np.array([dl(np.concatenate([idx[j] for j in range(20) if j != i])) for i in range(20)])
    sg = np.sqrt(19 / 20 * ((ds - ds.mean()) ** 2).sum())
    tot.append((m, s, d0, sg, n))
    print(f"{m}{s}: n={n} Delta={d0:+.4f} +- {sg:.4f}")
a = [t for t in tot if t[0] == 'A']
d = 0.5 * (a[0][2] + a[1][2]); s = 0.5 * np.hypot(a[0][3], a[1][3])
print(f"A 2-group: Delta={d:+.4f} +- {s:.4f}  ({d/s:+.2f} sigma)  spread={abs(a[0][2]-a[1][2]):.4f} "
      f"(={abs(a[0][2]-a[1][2])/np.hypot(a[0][3],a[1][3]):.2f} sigma_diff)")
print(f"need sigma<={abs(d)/3:.4f} => samples x{(s/(abs(d)/3))**2:.2f} => per-group sweeps "
      f"{a[0][4]*(s/(abs(d)/3))**2:.0f}")
