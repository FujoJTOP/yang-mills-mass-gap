from mpmath import mp, besseli, mpf
mp.dps = 40
print("=== is K = prod_i t_{j_i} bounded?  t_j = I_{j+1}(b)/I_1(b) ===")
print()
print("t_j values (beta_t = 2.5):")
for j in range(6):
    print("   j=%d  t_j=%.12f  <=1? %s"%(j, float(besseli(j+1,mpf("2.5"))/besseli(1,mpf("2.5"))),
          besseli(j+1,mpf("2.5"))<=besseli(1,mpf("2.5"))))
print()
print("large-beta asymptotics of t_j:")
print("  I_nu(b) ~ e^b/sqrt(2 pi b) [1 - (4 nu^2-1)/(8b) + ...]")
print("  t_j = I_{j+1}/I_1 ~ 1 - [ (4(j+1)^2-1) - (4-1) ] / (8b)")
print("       = 1 - [ 4j^2+8j ] / (8b) = 1 - j(j+2)/(2b)")
print()
print("check vs numerics:")
for b in ("10","100","1000"):
    bb=mpf(b)
    for j in (1,2,3):
        t=besseli(j+1,bb)/besseli(1,bb)
        pred=1-mpf(j*(j+2))/(2*bb)
        print("   b=%-6s j=%d  t_j=%.10f  1-j(j+2)/(2b)=%.10f"%(b,j,float(t),float(pred)))
    print()
print("==> t_j <= 1 for all j, b  =>  K = prod_i t_{j_i} <= 1  BOUNDED")
print("    (so 'K ~ beta unbounded' is FALSE for this normalization)")
