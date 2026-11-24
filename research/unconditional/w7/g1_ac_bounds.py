"""
g1_ac_bounds.py -- (G1) attack via 4D trace-anomaly central charges.

Pure-YM UV central charges (free-field sum over n_V = N^2-1 gauge vectors):
    a_V = 62/360,  c_V = 12/60 = 1/5   (per gauge vector)
Unitarity / conformal-collider bounds for a 4D unitary CFT:
    1/3 <= a/c <= 31/18           (free real scalar at 1/3)

a-theorem (Komargodski-Schwimmer):  a_IR <= a_UV for any RG flow.

Question: do these EXCLUDE a gapless CFT IR of pure YM (i.e. do they prove (G1))?
Answer (this witness): NO -- the UV data lie strictly inside the allowed a/c band,
and a_IR in (0, a_UV] is an interval, so a hypothetical IR CFT is NOT excluded.

Run: py -3.12 research/unconditional/w7/g1_ac_bounds.py
"""
# per-field trace-anomaly coefficients (normalization a = ... E_4, c = ... W^2)
a_vec, c_vec = 62/360, 12/60          # free real gauge vector
a_scalar, c_scalar = 1/360, 1/120      # free real scalar  (a/c = 1/3)
a_lo, a_hi = 1/3, 31/18                # unitarity / collider bounds on a/c

print("free-vector a/c =", a_vec/c_vec, " (bound [%s,%s])" % (a_lo, a_hi))
print("free-vector a/c in band:", a_lo <= a_vec/c_vec <= a_hi)
print()

for N in (2, 3, 4, 10):
    nV = N*N - 1
    aUV = nV * a_vec
    cUV = nV * c_vec
    print("SU(%2d): nV=%3d  a_UV=%.4f  c_UV=%.4f  a/c=%.4f  in[1/3,31/18]=%s"
          % (N, nV, aUV, cUV, aUV/cUV, a_lo <= aUV/cUV <= a_hi))
print()

# The a-theorem permits a_IR in (0, a_UV]; the ratio bound permits a_IR/c_IR in [1/3,31/18].
# A hypothetical IR CFT therefore has a whole 2D allowed region -- NOT a single point.
N = 3
aUV, cUV = (N*N-1)*a_vec, (N*N-1)*c_vec
print("SU(3) allowed IR-CFT region:")
print("   a_IR in (0, %.4f]  (a-theorem ceiling)" % aUV)
print("   a_IR/c_IR in [%.4f, %.4f]  (unitarity band)" % (a_lo, a_hi))
print("   e.g. a free massless scalar IR: a=%.5f c=%.5f a/c=%.4f -- satisfies both"
      % (a_scalar, c_scalar, a_scalar/c_scalar))
print()
print("VERDICT: central charges do NOT exclude a gapless CFT IR of pure YM.")
print("         The UV point (a/c = 31/36) sits inside [1/3, 31/18]; a_IR is an")
print("         interval, not forced to a value.  (G1) is NOT proven by a/c.")
print("         NB the only excluded CFT IR would be the *trivial* g*=0 one")
print("         (a_IR = a_UV, i.e. no flow) -- but that is the AF/no-IR-free fact,")
print("         not a central-charge bound.")
