import math

B4 = [0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.68, 0.72]
B6 = [0.36, 0.38, 0.40, 0.42, 0.43, 0.4407, 0.45, 0.46, 0.48, 0.50, 0.54]
TOTAL = 2000


def read(path, L):
    N = L ** 4
    P, W, A = {}, {}, {}
    for line in open(path):
        p = line.strip().split(',')
        if not p:
            continue
        if p[0] == 'P':
            P[int(p[1])] = (float(p[2]), float(p[3]), int(p[4]))
        elif p[0] == 'W':
            W[(int(p[1]), int(p[2]), int(p[3]))] = (float(p[4]), float(p[5]), int(p[6]))
        elif p[0] == 'A':
            A[int(p[1])] = int(p[2])
    return N, P, W, A


def stat(s, s2, n):
    m = s / n
    var = max(0.0, s2 / n - m * m)
    return m, (math.sqrt(var / n) if n > 1 else 0.0)


summary = open('results.csv', 'w')
summary.write('L,beta,plaquette,plaq_err,acceptance\n')
loops = open('wilson_loops.csv', 'w')
loops.write('L,beta,R,T,W,err\n')
creutz = open('creutz.csv', 'w')
creutz.write('L,beta,R,T,chi\n')

for L, betas, path in [(4, B4, 'results_l4c.csv'), (6, B6, 'results_l6c.csv')]:
    N, P, W, A = read(path, L)
    NLINK = 4 * N

    def wmean(c, r, t):
        s, s2, n = W[(c, r, t)]
        return stat(s, s2, n)[0] / N

    for c in sorted(P):
        s, s2, n = P[c]
        m, se = stat(s, s2, n)
        summary.write(f'{L},{betas[c]},{m/(6.0*N):.6f},{se/(6.0*N):.6f},{A.get(c,0)/(TOTAL*NLINK):.4f}\n')
        for r in (1, 2, 3):
            for t in (1, 2, 3):
                s, s2, n = W[(c, r, t)]
                mm, ee = stat(s, s2, n)
                loops.write(f'{L},{betas[c]},{r},{t},{mm/N:.6f},{ee/N:.6f}\n')
        for r in (2, 3):
            for t in (2, 3):
                num = wmean(c, r, t) * wmean(c, r-1, t-1)
                den = wmean(c, r-1, t) * wmean(c, r, t-1)
                chi = -math.log(num/den) if num > 0 and den > 0 else float('nan')
                creutz.write(f'{L},{betas[c]},{r},{t},{chi:.6f}\n')

summary.close(); loops.close(); creutz.close()
print('wrote results.csv, wilson_loops.csv, creutz.csv')
