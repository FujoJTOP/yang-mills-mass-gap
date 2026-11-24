"""V4 pincer witness -- shell-index bookkeeping for the coverage table.

Deterministic. numpy only. Run:  py -3.12 v4_pincer.py

It does NOT prove anything; it evaluates the campaign's *one-loop tower
parametrisation* (C2 / G1-L3-TOWER-POINCARE) and the *thresholds quoted in
S6/SZZ/G1-L3-UV-FORCE-BOUND*, so the shell ranges in V4-PINCER.md are
arithmetic, not assertion.  Every constant that is a physical input is
named; the two O(1) cone constants (kappa, cL) are placeholders scanned
over a range so the conclusion is shown to be inscription-independent.
"""
import numpy as np

# ---- tower parametrisation (C2 §2.1, G1-L3-TOWER-POINCARE §2) ----
#   g_j^2 = 1/(2 beta0 j) ,  beta_j = 2N/g_j^2 = 4 N beta0 j   (standard Wilson)
#   beta0 = 11 N / (48 pi^2)
def beta0(N):    return 11.0*N/(48.0*np.pi**2)
def beta_j(N,j): return 4.0*N*beta0(N)*j
def g_j(N,j):    return np.sqrt(1.0/(2.0*beta0(N)*j))

# ---- thresholds quoted in S6 / SZZ / C1 (standard-Wilson normalisation) ----
B_DOBR = 1.0/72.0      # Dobrushin alpha<1 (S6 §2)
B_SZZ  = 0.19          # SZZ beta<1/48 mapped to standard Wilson (S6 §2 table)

def first_shell_above(N, thr):
    """smallest integer j>=1 with beta_j > thr (None if none for j>=1)."""
    j = 1
    while j < 400:
        if beta_j(N,j) > thr:
            return j
        j += 1
    return None

def cone_empty_shell(N, L=2.0, a0=1.0, kappa=1.0, cL=1.0, csav=1.0):
    """smallest j with the Savvidy negative term < the Laplacian floor.

    G1-L3-UV-FORCE-BOUND §3: |F_j| = O(g_j)=t_j; the cone empties when
    sqrt(2 g_j |F_j|) << Lambda_{j+1}; equivalently the negative Savvidy
    term (~ csav*kappa*g_j^2) is beaten by the box Laplacian floor
    (~ cL / a_j^2), a_j = a0 L^{-j}.  Return minimal integer j.
    """
    j = 1
    while j < 400:
        neg = csav*kappa*g_j(N,j)**2
        floor = cL/(a0*L**(-j))**2
        if neg < floor:
            return j, neg, floor
        j += 1
    return None, None, None

def main():
    print("="*74)
    print("V4 pincer -- shell-index bookkeeping (one-loop tower, standard Wilson)")
    print("="*74)

    for N in (2,3):
        print(f"\n----- N = {N} : beta0 = {beta0(N):.5f},  beta_j = {4*N*beta0(N):.4f} * j -----")
        # strong-coupling coverage
        jd = first_shell_above(N, B_DOBR)
        js = first_shell_above(N, B_SZZ)
        print(f"  beta_1 = {beta_j(N,1):.4f}")
        print(f"  beta_Dobr = 1/72 = {B_DOBR:.4f}  -> first shell with beta_j>beta_Dobr: j={jd}")
        print(f"  beta_SZZ  ~= 0.19 = {B_SZZ:.4f}  -> first shell with beta_j>beta_SZZ : j={js}")
        n_sc_dobr = sum(1 for j in range(1,41) if beta_j(N,j) < B_DOBR)
        n_sc_szz  = sum(1 for j in range(1,41) if beta_j(N,j) < B_SZZ)
        print(f"  integer shells j>=1 with beta_j < 1/72 : {n_sc_dobr}   (<= 40 shells scanned)")
        print(f"  integer shells j>=1 with beta_j < 0.19 : {n_sc_szz}")

        # head / intermediate / UV in the one-loop index
        jhead = 2.5/(4*N*beta0(N))
        print(f"  beta_j = 2.5 at j = {jhead:.2f}  (the 'head'); j=1 gives beta_1={beta_j(N,1):.2f}")
        print("  cone-empty shells (scan of the O(1) constants):")
        for kappa in (0.5,1.0,2.0):
            for cL in (0.5,1.0,2.0):
                juv,neg,fl = cone_empty_shell(N, kappa=kappa, cL=cL)
                print(f"    kappa={kappa:<4} cL={cL:<4} -> j_UV={juv}  "
                      f"(savvidy {neg:.3e} < floor {fl:.3e})")

        # tower with a concrete depth
        L = 2.0
        for depth in (6, 10, 14):     # Lambda/a = L^depth
            J = depth
            juv,_,_ = cone_empty_shell(N, kappa=1.0, cL=1.0)
            uv_shells = [j for j in range(1,J+1) if j >= juv]
            sc_shells = [j for j in range(1,J+1) if beta_j(N,j) < B_SZZ]
            allshells = set(range(1,J+1))
            covered = set(sc_shells) | set(uv_shells)
            gap = sorted(allshells - covered)
            print(f"  tower J={J}: |SZZ shells|={len(sc_shells)}, |UV shells|={len(uv_shells)}, "
                  f"uncovered={len(gap)} -> j in {gap if len(gap)<=12 else str(gap[:6])+'...'+str(gap[-3:])}")

    print("\n" + "="*74)
    print("READ-OFF (vindicated by the numbers, not by this script):")
    print(" * strong-coupling jaws (1/72, 1/48->0.19) are crossed at j<1: ZERO tower shells.")
    print(" * UV leg turns on at j_UV = O(1) (1..3 in the physical scans); it covers")
    print("   every shell j>=j_UV, i.e. ALL but O(1) shells of an arbitrarily deep tower.")
    print(" * the uncovered set is the head band j=1..j_UV-1: a FINITE, J-independent")
    print("   number of O(1)-coupling shells sitting at beta in [beta_1, few].")
    print(" * hence there is NO large intermediate shell gap; the whole tower above")
    print("   O(1) shells is UV-leg territory, and the head is exactly the IR bottom.")
    print("="*74)

if __name__ == "__main__":
    main()
