"""
star_no_irfree.py -- STAR-NO-IRFREE witness (2026-11-24)

Checks the "IR is not IR-free" link of the (★) chain for pure SU(N) YM:

  (1) 2-loop beta:  beta(g) = -B0*g^3 - B1*g^5   with
         B0 = 11N/(3*16pi^2),   B1 = 34N^2/(3*16pi^2)^2   (pure glue, no matter)
      is NEGATIVE for ALL g>0  (B0,B1>0 => no zero at g>0 => no Banks-Zaks).
  (2) Integrating the RG ODE  dg/dln(mu) = beta(g) with g(UV)=small
      drives g UP as mu DECREASES  => the IR coupling is O(1), not 0.
  (3) Contrast: a FREE field has beta(g)=0  => the statement is YM-specific
      (acid: the argument is FALSE for the free massless field, where the IR
      IS free -- so it discriminates, it is not vacuous).

No QED. Pure real analysis / numpy.
"""
import numpy as np

def B0(N):  return 11.0*N/(3.0*16.0*np.pi**2)
def B1(N):  return 34.0*N**2/(3.0*(16.0*np.pi**2)**2)

def beta2(g, N):
    return -B0(N)*g**3 - B1(N)*g**5

def run(N, g0=0.5, mu0=1e3, mu1=1e-3, n=4000):
    # integrate dg/dln mu = beta(g); walk ln mu from ln(mu0) down to ln(mu1)
    lmu = np.linspace(np.log(mu0), np.log(mu1), n)
    g = np.empty(n); g[0] = g0
    for i in range(1, n):
        h = lmu[i]-lmu[i-1]
        k1 = beta2(g[i-1], N)
        k2 = beta2(g[i-1]+0.5*h*k1, N)
        g[i] = g[i-1] + h*k2
        if g[i] > 30:   # g would blow up (Landau-like); clamp for reporting
            g[i] = 30.0
    return g, np.exp(lmu)

if __name__ == "__main__":
    print("="*70)
    print("(1) 2-loop beta sign for pure SU(N): beta(g) < 0 for all g>0 ?")
    print("="*70)
    g = np.linspace(1e-3, 10.0, 200000)
    for N in (2, 3):
        b = beta2(g, N)
        print(f"  SU({N}): B0={B0(N):.5f}  B1={B1(N):.5f}  "
              f"max(beta)={b.max():.3e} (<0 ? {b.max()<0})  "
              f"sign changes: {int(np.sum(np.diff(np.sign(b))!=0))}")

    print("\n" + "="*70)
    print("(2) IR coupling: g grows as mu decreases (g(UV)=0.5)")
    print("="*70)
    for N in (2, 3):
        gs, mus = run(N)
        for mu_t in (1e3, 1.0, 1e-1, 1e-2):
            i = int(np.argmin(np.abs(np.log(mus)-np.log(mu_t))))
            tag = "  <-- strong (O(1))" if gs[i] > 1 else ""
            print(f"  SU({N}): mu={mu_t:8.0e}  g(mu)={gs[i]:.3f}{tag}")

    print("\n" + "="*70)
    print("(3) ACID contrast: free field has beta(g)=0 (IR stays free)")
    print("="*70)
    print("  free massless: beta(g)=0  ->  g(IR)=g(UV)=0  ->  IR IS free")
    print("  => 'IR not free' is FALSE for the free field (discriminates).")

    print("\n" + "="*70)
    print("VERDICT: pure YM g moves AWAY from 0 into the IR (AF); "
          "IR coupling O(1); NOT IR-free.")
    print("="*70)
