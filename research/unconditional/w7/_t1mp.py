from mpmath import mp, besseli, mpf, log
import mpmath
mp.dps = 80
print("=== T1 independent verification with mpmath (80 digits) ===")
print("    I_1(z) I_3(z) < I_2(z)^2  for all z > 0")
print()
# dense log-spaced grid plus large-z checks
import math
zs = [mpf("1e-6") * (mpf(10) ** (mpf(k)/mpf(40))) for k in range(0, 41)]   # 1e-6 .. 1e2
zs += [mpf(z) for z in (0.001,0.01,0.1,0.5,1,2,3,4,5,10,20,29,50,100,200,500,1000,1e4,1e5,1e6)]
zs = sorted(set(zs))
bad = []
mindiff = None
minz = None
for z in zs:
    I1 = besseli(1, z)
    I2 = besseli(2, z)
    I3 = besseli(3, z)
    D = I2*I2 - I1*I3
    if D <= 0:
        bad.append((z, D))
    if mindiff is None or D < mindiff:
        mindiff = D; minz = z
print("grid points: %d" % len(zs))
print("violations: %d" % len(bad))
if bad:
    for z,D in bad[:10]:
        print("   VIOLATION at z=%s : D=%s" % (z, D))
print("min of D(z) over the grid: D=%.3e at z=%s" % (float(mindiff), minz))
print()
# asymptotic check: D(z) ~ ? at large z
print("asymptotic check (large z):")
for z in (100, 1000, 1e4, 1e5):
    zz = mpf(z)
    I1 = besseli(1, zz); I2 = besseli(2, zz); I3 = besseli(3, zz)
    D = I2*I2 - I1*I3
    # ratio D/(I2^2)
    print("   z=%-8g  D/I2^2 = %.6e   D/(I2^2 z^-1) = %.6e" % (z, float(D/(I2*I2)), float(D/(I2*I2)*zz)))
print()
# small z check: D(z) ~ ? 
print("small z check:")
for z in (1e-6, 1e-4, 1e-2):
    zz = mpf(z)
    I1 = besseli(1, zz); I2 = besseli(2, zz); I3 = besseli(3, zz)
    D = I2*I2 - I1*I3
    print("   z=%-8g  D = %.6e   D/z^2 = %.6e" % (z, float(D), float(D/(zz*zz))))
