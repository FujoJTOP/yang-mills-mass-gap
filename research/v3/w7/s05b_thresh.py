"""s05b_thresh.py — where does 2D-YM G'' cross zero in t?  (t = g^2 A / 2)"""
import numpy as np
from s05_acid import two_d

phi = np.linspace(0.05, np.pi - 0.05, 401)
mn = lambda t: min(two_d(p, t)[0] for p in phi)
prev = None
for t in np.arange(0.4, 3.01, 0.1):
    v = mn(t)
    if prev is not None and prev[1] < 0 <= v:
        print(f"sign change of min_phi G'' between t={prev[0]:.2f} and t={t:.2f}")
    prev = (t, v)
for t in (0.5, 1.0, 1.3, 1.5, 1.7, 2.0, 5.0):
    print(f"  t={t:4.1f}  min_phi G''(2D) = {mn(t):+9.5f}")
