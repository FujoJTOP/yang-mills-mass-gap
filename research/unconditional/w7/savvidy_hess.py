import numpy as np

# Savvidy point: plaquette with q=0, |a|=1 (pure transverse magnetic), a = e_3
# Exact formula (12.1): Q_P(Y) = beta*q*|sum y_i|^2 - 2*beta*sum_{i<j} (y_i x y_j).a
# Gauge directions at a plaquette with links (1,2,3,4) in order:
#   U_1(x)->g_x U_1(x) g_{x+1}^-1  => y_1 -> y_1 + v_x - v_{x+1}
# For a translation-invariant gauge transformation v(x) = eps*x, the 4 links get
#   signs (+,-,+,-) up to ordering.  We test the pure "all-shift" mode and the
#   alternating mode explicitly.

a = np.array([0.0, 0.0, 1.0])
beta = 2.5
q = 0.0

def QP(Y):
    Y = np.asarray(Y, float).reshape(4, 3)
    s = Y.sum(axis=0)
    tot = beta * q * (s @ s)
    cr = np.zeros(3)
    for i in range(4):
        for j in range(i + 1, 4):
            cr += np.cross(Y[i], Y[j])
    tot += -2.0 * beta * (cr @ a)
    return tot

print("=== Hess S at the Savvidy point (q=0, |a|=1), per single plaquette ===")
print("  Q_P(Y) = -2 beta * sum_{i<j} (y_i x y_j).a        [q=0 branch]")
print()

# candidate modes
modes = {
    "all-shift  (v,v,v,v)        ": np.tile([1.0, 0, 0], 4),
    "alternating (v,-v,v,-v)     ": np.tile([1.0, 0, 0], 4) * np.repeat([1, -1, 1, -1], 3),
    "transverse e2 pattern       ": np.tile([0, 1.0, 0], 4),
    "rotated (y_i = R_i a)       ": None,
}
for name, Y in modes.items():
    if Y is None:
        continue
    Y = np.asarray(Y, float)
    n2 = (Y.reshape(4, 3) ** 2).sum()
    print("  %s : Q_P = %+10.6f   |Y|^2 = %.6f   ratio = %+8.4f"
          % (name, QP(Y), n2, QP(Y) / n2))

# scan: maximize the negative eigenvalue of the quadratic form on R^12
H = np.zeros((12, 12))
E = np.zeros(12)
for k in range(12):
    E[k] = 1.0
    Yp = E.copy(); Ym = E.copy(); Yp[k] += 1e-6; Ym[k] -= 1e-6
    H[k, k] = (QP(Yp) - 2 * QP(np.zeros(12)) + QP(Ym)) / 1e-12
    E[k] = 0.0
for k in range(12):
    for l in range(k + 1, 12):
        Ep = np.zeros(12); Ep[k] += 1e-6; Ep[l] += 1e-6
        Em = np.zeros(12); Em[k] += 1e-6; Em[l] -= 1e-6
        En = np.zeros(12); En[k] -= 1e-6; En[l] += 1e-6
        Emm = np.zeros(12); Emm[k] -= 1e-6; Emm[l] -= 1e-6
        v = (QP(Ep) - QP(Em) - QP(En) + QP(Emm)) / (4 * 1e-12)
        H[k, l] = H[l, k] = v

ev = np.linalg.eigvalsh(H)
print()
print("  spectrum of Hep_S on R^12 (4 links x so(3)):")
print("   ", np.array2string(ev, precision=4, suppress_small=False))
print("  min eigenvalue = %+.6f    max = %+.6f" % (ev.min(), ev.max()))
print()
print("  gauge-zero modes: infinitesimal gauge transformations shift the 4 link")
print("  variables by (+v,-v,+v,-v) (up to link orientation), so the gauge subspace")
print("  is spanned by 3 vectors; check Q_P on them:")
for d in range(3):
    Yg = np.zeros(12)
    for i, s in enumerate([1, -1, 1, -1]):
        Yg[3 * i + d] = s
    print("    dir %d : Q_P = %+10.6f   |Y|^2 = %.6f   ratio = %+8.4f"
          % (d, QP(Yg), (Yg ** 2).sum(), QP(Yg) / (Yg ** 2).sum()))
print()
print("  => if the gauge directions have ratio 0 and the minimum is negative,")
print("     the negative mode is TRANSVERSE (physical).")
