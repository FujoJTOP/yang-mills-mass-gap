import numpy as np

# Candidate lemma:  c_d > 0 and log-concave  ==>  (M_f)_{j'j} = sum_{d=|j-j'|}^{j+j'} c_d  is TN2.
# Test on many random log-concave families (not just geometric).

def build(c, Jf):
    M = np.zeros((Jf + 1, Jf + 1))
    for ia in range(Jf + 1):
        for ib in range(Jf + 1):
            ja, jb = ia / 2.0, ib / 2.0
            lo, hi = abs(ja - jb), ja + jb
            n = int(round(hi - lo))
            M[ia, ib] = sum(c[int(round(lo)) + t] for t in range(n + 1))
    return M

def min_minor(M, Jf):
    mn = np.inf
    for i in range(Jf):
        for j in range(i + 1, Jf + 1):
            for k in range(Jf):
                for l in range(k + 1, Jf + 1):
                    mn = min(mn, M[i, k] * M[j, l] - M[i, l] * M[j, k])
    return mn

def logconcave(v):
    return all(v[k] ** 2 >= v[k - 1] * v[k + 1] - 1e-15 * v[k] ** 2 for k in range(1, len(v) - 1))

rng = np.random.default_rng(3)
Jf = 10

print("=== random LOG-CONCAVE c (built as ratios r_d in (0,1] decreasing) ===")
fails = 0; tot = 0
for trial in range(40):
    # c_d determined by decreasing ratios q_d = c_{d+1}/c_d in (0,1]; log-concavity == q_d decreasing
    q = np.sort(rng.uniform(0.05, 0.95, Jf))[::-1]
    c = np.ones(Jf + 1)
    for d in range(Jf):
        c[d + 1] = c[d] * q[d]
    if not logconcave(c):
        continue
    M = build(c, Jf)
    mn = min_minor(M, Jf)
    tot += 1
    ok = mn >= -1e-10 * M.max() ** 2
    if not ok: fails += 1
print("  log-concave families tested = %d ; TN2 failures = %d" % (tot, fails))
print()

print("=== random NON-log-concave c (ratios NOT decreasing) ===")
fails2 = 0; tot2 = 0
for trial in range(40):
    q = rng.uniform(0.05, 0.95, Jf)          # unsorted ratios -> generally not log-concave
    c = np.ones(Jf + 1)
    for d in range(Jf):
        c[d + 1] = c[d] * q[d]
    if logconcave(c):
        continue
    M = build(c, Jf)
    mn = min_minor(M, Jf)
    tot2 += 1
    ok = mn >= -1e-10 * M.max() ** 2
    if not ok: fails2 += 1
print("  non-log-concave families tested = %d ; TN2 failures = %d" % (tot2, fails2))
