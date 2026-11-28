import numpy as np
from scipy.special import iv

# CRITICAL CHECK (round 42):  t_j = 2 I_{2j+1}(beta)/beta is NOT decreasing in j for large beta.
# I_n(z) peaks near n ~ z, so for beta = 40, t_j INCREASES up to j ~ 20.
# Hence  r(K) = (second largest of t_j)/(largest)  is NOT I_2/I_1 for large beta.
# Recompute r(K) correctly and redo the comparison.

print("=== is t_j decreasing in j ? ===")
for beta in (0.5, 1.0, 2.0, 2.5, 5.0, 10.0, 20.0, 40.0):
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(0, min(120, int(4 * beta) + 40))])
    imax = int(np.argmax(t))
    print("  beta=%5.1f : argmax t_j at j=%.1f  (t_0=%.6e, t_max=%.6e, t_max/t_0=%.6f)"
          % (beta, imax / 2.0, t[0], t[imax], t[imax] / t[0]))

print()
print("=== correct r(K): ratio of the two largest values of t_j ===")
print("%8s %14s %14s %14s" % ("beta", "I2/I1", "r(K) correct", "ratio r(K)/(I2/I1)"))
for beta in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0, 20.0, 40.0):
    N = int(6 * beta) + 60
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(N)])
    tt = np.sort(t)[::-1]
    rK_correct = tt[1] / tt[0]
    rK_mine = float(iv(2, beta) / iv(1, beta))
    print("%8.1f %14.10f %14.10f %14.10f" % (beta, rK_mine, rK_correct, rK_correct / rK_mine))
