import csv, math, sys

def bessel_i(n, x):
    # I_n(x) = sum_k (x/2)^(2k+n) / (k! (k+n)!)
    half = x / 2.0
    term = half ** n / math.factorial(n)
    s = term
    for k in range(1, 60):
        term *= (half * half) / (k * (k + n))
        s += term
        if abs(term) < 1e-16 * abs(s):
            break
    return s

def exact_sigma(beta):
    return -math.log(bessel_i(1, beta) / bessel_i(0, beta))

def exact_plaq(beta):
    return bessel_i(1, beta) / bessel_i(0, beta)

def main(path):
    W = {}
    plaq = {}
    acc = {}
    with open(path, newline='') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if line.startswith('# acc'):
                    pass
                continue
            parts = line.split(',')
            if len(parts) != 5:
                continue
            try:
                beta = float(parts[0]); R = int(parts[1]); T = int(parts[2])
                val = float(parts[3]); err = float(parts[4])
            except ValueError:
                continue
            if R == 0 and T == 0:
                plaq[beta] = (val, err)
            else:
                W[(beta, R, T)] = (val, err)

    print("beta  plaquette(meas)  exact I1/I0   sigma_exact   sigma_from_W(1,1)")
    for beta in sorted(plaq):
        p, pe = plaq[beta]
        se = exact_sigma(beta)
        w11 = W.get((beta, 1, 1))
        s11 = -math.log(w11[0]) if w11 and w11[0] > 0 else float('nan')
        print(f"{beta:4.1f}  {p:8.5f}±{pe:.5f}   {exact_plaq(beta):8.5f}   {se:8.5f}     {s11:8.5f}")

    print("\nCreutz ratios chi(R,T) (should approach sigma_exact at large R,T)")
    for beta in sorted(plaq):
        se = exact_sigma(beta)
        print(f"\nbeta = {beta}   sigma_exact = {se:.5f}")
        for R in range(2, 5):
            for T in range(2, 5):
                try:
                    w = W[(beta, R, T)][0]
                    w11 = W[(beta, R-1, T-1)][0]
                    w10 = W[(beta, R-1, T)][0]
                    w01 = W[(beta, R, T-1)][0]
                    if min(w, w11, w10, w01) > 0:
                        chi = -math.log(w * w11 / (w10 * w01))
                        print(f"  chi({R},{T}) = {chi:8.5f}")
                except KeyError:
                    pass

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'results.csv')
