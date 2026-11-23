import math, sys

BETAS = [0.65, 0.67, 0.69, 0.71, 0.73, 0.75, 0.7613, 0.77, 0.79, 0.81, 0.85]
if len(sys.argv) > 3:
    BETAS = [float(x) for x in sys.argv[3].split(',')]

def analyze(path, L):
    N = L**4
    NLINK = 4 * N
    TOTAL = 2000
    P, W, A = {}, {}, {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        p = line.split(',')
        if p[0] == 'P':
            c = int(p[1]); s = float(p[2]); s2 = float(p[3]); n = int(p[4])
            P[c] = (s, s2, n)
        elif p[0] == 'W':
            c, r, t = int(p[1]), int(p[2]), int(p[3])
            s, s2, n = float(p[4]), float(p[5]), int(p[6])
            W[(c, r, t)] = (s, s2, n)
        elif p[0] == 'A':
            A[int(p[1])] = int(p[2])

    def stat(s, s2, n):
        m = s / n
        var = max(0.0, s2 / n - m * m)
        se = math.sqrt(var / n) if n > 1 else 0.0
        return m, se

    print(f"L={L}  N={N}  nconf={len(P)}")
    print("beta      plaquette        acc")
    for c in sorted(P):
        s, s2, n = P[c]
        m, se = stat(s, s2, n)
        plaq = m / (6.0 * N)
        err = se / (6.0 * N)
        acc = A.get(c, 0) / (TOTAL * NLINK)
        print(f"{BETAS[c]:<8.5g}  {plaq:.5f} +/- {err:.5f}   {acc:.3f}")

    print("\nWilson loops W(R,T)")
    print("beta      " + "  ".join(f"({r},{t})" for r in (1,2,3) for t in (1,2,3)))
    for c in sorted(P):
        vals = []
        for r in (1,2,3):
            for t in (1,2,3):
                s, s2, n = W[(c,r,t)]
                m, se = stat(s, s2, n)
                vals.append(f"{m/N:.5f}")
        print(f"{BETAS[c]:<8.5g}  " + "  ".join(vals))

    print("\nCreutz ratios chi(R,T)")
    print("beta      chi(2,2)   chi(2,3)   chi(3,2)   chi(3,3)")
    for c in sorted(P):
        def w(r,t):
            s,s2,n = W[(c,r,t)]
            return s/(N*n)
        def chi(r,t):
            num = w(r,t)*w(r-1,t-1)
            den = w(r-1,t)*w(r,t-1)
            if num <= 0 or den <= 0:
                return float('nan')
            return -math.log(num/den)
        print(f"{BETAS[c]:<8.5g}  {chi(2,2):8.4f}   {chi(2,3):8.4f}   {chi(3,2):8.4f}   {chi(3,3):8.4f}")

if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'results_raw_l4.csv'
    L = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    analyze(path, L)
