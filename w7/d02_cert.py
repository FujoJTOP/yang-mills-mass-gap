# d02_cert.py -- D02-DIRECT: what part of "strict certificate at beta=2.5" is rigorous today
# NOTE on characters: SU(2) irreps are labelled by HALF-integer j >= 0, chi_j(th) = sin((2j+1)th)/sin(th),
# dim = 2j+1 =: m in {1,2,3,...}.  The Wilson weight e^{(beta/2)tr U} = e^{beta cos th} needs ALL m,
# including the odd ones (j = 1/2, 3/2, ...).  Using integer j only is WRONG (it spans SO(3) characters).
# Contents:
#   [1] rigorous character coefficients kappa_m (closed form + quadrature check + reconstruction)
#   [2] rigorous truncation tail bound (bare Bessel bound) -> truncation level
#   [3] EXACT 1-link (convolution) model spectrum  -> a genuine strict number, but a toy
#   [4] calibration: for a 2x2 positive matrix, lambda2/lambda1 == tanh(Delta/4) EXACTLY
#   [5] structural identity: spatial-plaquette factors cancel in the Hilbert/Dobrushin ratio
#   [6] boundedness (hence Hilbert-Schmidt, hence compact) of T at L=2
import sys
import mpmath as mp

mp.mp.dps = 30

def chi(j, th):
    """SU(2) character of spin j at angle th (U ~ diag(e^{i th}, e^{-i th}))."""
    s = mp.sin(th)
    if abs(s) < mp.mpf('1e-40'):
        return mp.mpf(2*j+1)
    return mp.sin((2*j+1)*th)/s

def kappa_quad(j, beta):
    """kappa_j = int_{SU(2)} e^{beta cos th} chi_j(th) dU, normalized Haar: (2/pi) int_0^pi sin^2 th dth."""
    f = lambda th: mp.e**(beta*mp.cos(th))*chi(j, th)*mp.sin(th)**2
    return (2/mp.pi)*mp.quad(f, [0, mp.pi])

def main():
    beta = mp.mpf('2.5')
    print("="*78)
    print(f"D02 certificate inventory at beta = {beta}  (mpmath dps={mp.mp.dps})")
    print("="*78)

    # ---------- [1] character coefficients, half-integer j ----------
    print("\n[1] character coefficients kappa_j = int e^{beta cos th} chi_j dU   (closed form 2(2j+1) I_{2j+1}(beta)/beta = 2 m I_m/beta)")
    print("    irrep label j = (m-1)/2,  m = 2j+1 in {1,2,3,...}; the ODD m are the ones an integer-j scan MISSES")
    print(f"{'m':>2} {'j':>5} {'kappa (quad)':>26} {'2*I_m/beta':>26} {'rel.diff':>10}")
    kap = {}
    for m in range(1, 15):
        j = mp.mpf(m-1)/2
        q = kappa_quad(j, beta)
        cl = 2*m*mp.besseli(m, beta)/beta
        kap[m] = cl
        print(f"{m:>2} {mp.nstr(j,3):>5} {mp.nstr(q,20):>26} {mp.nstr(cl,20):>26} {mp.nstr(abs(q-cl)/cl,3):>10}")

    # identity sum_m m*kappa_m = e^beta  (this is the character reconstruction at th=0)
    for M in (5, 9, 13):
        S = sum(m*kap[m] for m in range(1, M+1))
        print(f"  sum_(m<={M:>2}) m*kappa_m = {mp.nstr(S,20):>22}  e^beta = {mp.nstr(mp.e**beta,20)}  rel = {mp.nstr(abs(S-mp.e**beta)/(mp.e**beta),3)}")

    print("  pointwise reconstruction  sum_{m<=M} kappa_m chi_{(m-1)/2}(th)  vs  e^{beta cos th}:")
    for M in (5, 9, 13):
        for th in ['0.3', '1.0', '2.0', '2.6']:
            thv = mp.mpf(th)
            rec = sum(kap[m]*chi(mp.mpf(m-1)/2, thv) for m in range(1, M+1))
            exact = mp.e**(beta*mp.cos(thv))
            print(f"    M={M:>2} th={th:>4}  rec={mp.nstr(rec,14):>18}  exact={mp.nstr(exact,14):>18}  rel={mp.nstr(abs(rec-exact)/exact,3)}")

    # ---------- [2] rigorous truncation tail ----------
    print("\n[2] RIGOROUS tail bound.  I_m(x) <= (x/2)^m e^{x^2/4}/m!  for x>0  (from I_m = sum_k (x/2)^{2k+m}/(k!(m+k)!))")
    print(f"    truncated weight sum_M := sum_(m<=M) m*kappa_m  (exact sum = sum_m 2 m^2 I_m/beta = e^beta)")
    print(f"{'M':>2} {'tail <= 2m^2(be/2)^m e^{b^2/4}/m!':>34} {'tail/e^beta':>14} {'sum_M/e^beta':>14}")
    for M in range(1, 9):
        tail = mp.mpf(0)
        for m in range(M+1, M+60):
            tail += 2*m**2*(beta/2)**m*mp.e**(beta**2/4)/mp.factorial(m)
        part = sum(m*kap[m] for m in range(1, M+1))/mp.e**beta
        print(f"{M:>2} {mp.nstr(tail,8):>34} {mp.nstr(tail/mp.e**beta,6):>14} {mp.nstr(part,10):>14}")
    print("    1 - sum_M/e^beta must be <= tail/e^beta :")
    for M in (3, 5, 7):
        tail = sum(2*m**2*(beta/2)**m*mp.e**(beta**2/4)/mp.factorial(m) for m in range(M+1, M+60))
        part = sum(m*kap[m] for m in range(1, M+1))/mp.e**beta
        print(f"    M={M}: 1-sum_M/e^beta = {mp.nstr(1-part,6):>12}   bound = {mp.nstr(tail/mp.e**beta,6):>12}   valid={1-part<=tail/mp.e**beta}")

    # ---------- [3] exact 1-link spectrum (TOY, but genuinely strict) ----------
    print("\n[3] EXACT spectrum of the 1-link convolution operator Tf(U)=int e^{beta cos(U V^dag)} f(V)dV")
    print("    Schur: D^j_{mn} are exact eigenvectors, eigenvalue = kappa_j/(2j+1) = 2 I_m/beta, m=2j+1.")
    print("    => lambda_m/lambda_1 = I_m/I_1   (closed form; multiplicity of j is (2j+1)^2 = m^2)")
    print(f"{'m':>2} {'j':>5} {'lambda_m/lambda_1':>20} {'m_m=-log(ratio)':>18} {'mult':>5}")
    for m in range(1, 8):
        r = mp.besseli(m, beta)/mp.besseli(1, beta)
        print(f"{m:>2} {mp.nstr(mp.mpf(m-1)/2,3):>5} {mp.nstr(r,18):>20} {mp.nstr(-mp.log(r),12):>18} {m*m:>5}")
    print(f"    KEY: m=2 (j=1/2, the FUNDAMENTAL) gives ratio = I_2/I_1 = {mp.nstr(mp.besseli(2,beta)/mp.besseli(1,beta),16)}")
    print(f"         = the corpus's u_tile(2.5) = 0.5071951000  ->  toy gap m_toy = {mp.nstr(-mp.log(mp.besseli(2,beta)/mp.besseli(1,beta)),12)}")
    print("    -> a STRICT closed-form gap, but for the decoupled-link toy, NOT for 4D.")

    # ---------- [4] calibration of the Hilbert/Dobrushin bound ----------
    print("\n[4] Calibration: symmetric 2x2 positive matrix [[a,b],[b,a]] -> lambda2/lambda1 = tanh(Delta/4)")
    print("    with Delta = 2 log(a/b) = projective (Hilbert) diameter.  EQUALITY in d=2")
    print("    => the Birkhoff/Hopf oscillation bound lambda2/lambda1 <= tanh(Delta/4) is SHARP in d=2.")
    for Jv in ['0.5', '2.0', '4.0', '8.0']:
        J = mp.mpf(Jv); b = mp.mpf(1); a = b*mp.e**J
        r_exact = (a-b)/(a+b)
        Delta = 2*mp.log(a/b)
        print(f"    J={Jv:>4}  exact ratio={mp.nstr(r_exact,12):>16}  tanh(Delta/4)={mp.nstr(mp.tanh(Delta/4),12):>16}  equal={abs(r_exact-mp.tanh(Delta/4))<mp.mpf('1e-25')}")

    # ---------- [5] structural identity: spatial factors cancel ----------
    print("\n[5] Structural identity (machine check on random reals):")
    print("    T[U,U'] = e^{-Ssp(U)/2} e^{-Ssp(U')/2} W(U,U')")
    print("    => T[U,U']T[V,V'] / (T[U,V']T[V,U']) = W[U,U']W[V,V']/(W[U,V']W[V,U'])   (spatial terms cancel identically)")
    import random
    random.seed(7)
    worst = mp.mpf(0)
    for _ in range(2000):
        Ssp = [mp.mpf(random.uniform(-8, 8)) for _ in range(4)]
        lgW = [[mp.mpf(random.uniform(-30, 0)) for _ in range(2)] for _ in range(2)]
        lhs = mp.e**(-Ssp[0]/2-Ssp[1]/2)*mp.e**(-Ssp[2]/2-Ssp[3]/2) * mp.e**(lgW[0][0])*mp.e**(lgW[1][1])
        den = mp.e**(-Ssp[0]/2-Ssp[3]/2)*mp.e**(-Ssp[2]/2-Ssp[1]/2) * mp.e**(lgW[0][1])*mp.e**(lgW[1][0])
        rhs = mp.e**(lgW[0][0]+lgW[1][1]-lgW[0][1]-lgW[1][0])
        worst = max(worst, abs(lhs/den-rhs))
    print(f"    2000 random draws: max |lhs-rhs| = {mp.nstr(worst,3)}  => identity exact")
    print("    Consequence: the Hilbert/Dobrushin oscillation bound for T depends ONLY on the 3D object W.")

    # ---------- [6] boundedness / compactness ----------
    print("\n[6] Boundedness at L=2 (=> T bounded, and Hilbert-Schmidt on a finite-measure space => compact)")
    L = 2
    plaq = 6*L**4
    print(f"    #plaquettes = 6L^4 = {plaq};  0 < T[U,U'] <= e^{{beta*plaq}} = e^({beta*plaq}) = {mp.nstr(mp.e**(beta*plaq),6)}")
    print("    kernel strictly positive (integrand of W > 0 on an open set)  => Jentzsch / Perron-Frobenius applies")
    print("    => lambda_1 simple, eigenvector > 0, |lambda_j| < lambda_1 for j >= 2  [QUALITATIVE ONLY: no number]")
    return 0

if __name__ == "__main__":
    sys.exit(main())
