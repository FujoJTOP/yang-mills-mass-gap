"""P19 fit — beta-trend of G''(0) at L=2 (same-model points, ns=5e3 x4 chains each)."""
import numpy as np
B = np.array([1.0, 2.5, 4.0, 6.0])
Y = np.array([0.2835, 0.4116, 0.5108, 0.5913])          # G''(0)
S = np.array([0.0049, 0.0421, 0.0939, 0.1732])          # sigma_jk
D = np.array([0.0148, 0.1401, 0.2221, 0.2962])          # Delta_dyn
V = np.array([0.2888, 1.5232, 2.8421, 4.7329])          # Var(F)
K = np.array([0.5723, 1.9348, 3.3529, 5.3242])          # E[kappa]
W = 0.2754
print("beta   Gpp(0)      s_jk    s_tau   D/S    e_k/b  V/b^2   Delta_dyn  Wgeo+Gpp-D")
for i, b in enumerate(B):
    print(f"{b:4.1f} {Y[i]:8.4f} {S[i]:8.4f}  {[0.0105,0.1028,0.2216,0.4663][i]:7.4f} "
          f"{V[i]/K[i]:6.3f} {K[i]/b:8.3f} {V[i]/b**2:7.4f} {D[i]:+9.4f}")

for deg, nm in ((2, 'A+Bb+Cb^2'), (3, 'cubic')):
    c = np.polyfit(B, Y, deg); p = np.poly1d(c)
    grid = np.linspace(1, 12, 2201); g = p(grid)
    i0 = np.argmin(g)
    print(f"\n{nm}: coef={np.round(c,5)}  min on [1,12] at beta={grid[i0]:.3f} -> {g[i0]:+.4f}")
    print(f"   values: " + " ".join(f"{b:g}:{p(b):+.4f}" for b in (1, 2.5, 4, 4.31, 6, 8, 10)))
    r = np.roots(np.trim_zeros(c) if deg == 2 else c)
    rr = [x.real for x in r if abs(x.imag) < 1e-9 and x.real > 1]
    print(f"   real zeros >1: {np.round(rr,3) if rr else 'none'}")

a, bb = np.polyfit(B, D, 1)
print(f"\nDelta_dyn linear in beta: {a:+.4f}*beta {bb:+.4f}")
X = np.vstack([B, -B ** 2]).T
sol = np.linalg.lstsq(X, D, rcond=None)[0]
print(f"Delta_dyn = a*beta - b*beta^2 : a={sol[0]:+.4f}  b={sol[1]:+.4f}  "
      f"(b>0 would be the NS beta^2 obstruction)")
print(f"residual {np.round(D - X @ sol, 4)}")
print(f"\npower laws (log-log over beta=1..6): "
      f"E[kap]~b^{np.polyfit(np.log(B),np.log(K),1)[0]:.2f}  "
      f"Var(F)~b^{np.polyfit(np.log(B),np.log(V),1)[0]:.2f}  "
      f"Gpp~b^{np.polyfit(np.log(B),np.log(Y),1)[0]:.2f}  "
      f"D/S~b^{np.polyfit(np.log(B),np.log(V/K),1)[0]:.2f}")
print(f"fit 1-D/S = c*beta^-p: c={np.exp(np.polyfit(np.log(B),np.log(1-V/K),1)[1]):.3f} "
      f"p={-np.polyfit(np.log(B),np.log(1-V/K),1)[0]:.3f}")
