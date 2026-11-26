# -*- coding: utf-8 -*-
"""
M3A-C  :  L3 P3 —— 验 `C <= 2.691`
Only writes: research/v3/species/M3A-C.md + w7/m3a_c_check.py|.log

Question: what is C in  m_c = a_c*ln(1/u) + O(1)  (PROOF-IN11.md:9),
and does  C <= 1/2 * ln(1/u)|_{beta=0.018394} = 2.691  hold?

The corpus has exactly ONE explicit additive O(1) term in a mass/rate formula:
    delta = ln(1/ubar)/a = kappa/4 - ln(eDelta*Cbar)/a ,  ubar = eDelta*Cbar*u
    (A02-2PPI2.md:60 ; Z07-HOLEB.md:100 ; Z12-FIXC.md:124)
with the re-coupling factor  ubar -> Chat*ubar  (Z12-FIXC.md:96, Chat=4):
    delta_fix = kappa/4 - ln(Chat*eDelta*Cbar)/a .
kappa/4 = ln(1/u)/a  is the bare straight-tube rate.  Hence

    C := ln( Chat * eDelta * Cbar )          [the additive, channel-blind constant]

Everything below is arithmetic with mpmath (50 digits).
"""
from mpmath import mp, mpf, log, exp, pi, sqrt

mp.dps = 50

E = exp(mpf(1))
EDELTA = mpf('54.3656365692')     # e*Delta, Delta=20 (P3-FULLREP-KP:89, A01-CBAR:33)
CHAT = mpf(4)                     # Z12-FIXC:70, provable per-plaquette envelope
DELTA = mpf(20)
BETA1 = mpf('0.018394')           # window: beta < 8/(Chat*Cbar*eDelta)

def hdr(s):
    print("\n" + "=" * 74 + "\n" + s + "\n" + "=" * 74)

hdr("[0] inputs")
print("  eDelta  = 54.3656365692   (Delta = 20, eDelta/e = %s)" % mp.nstr(EDELTA / E, 12))
print("  Chat    = 4   (Z12-FIXC:70)")
print("  beta1   = 0.018394  claimed = 8/(Chat*Cbar*eDelta)")

hdr("[1] window identity: beta1 = 8/(Chat*Cbar*eDelta)  <->  ubar=Chat*Cbar*eDelta*u = 1")
for Cbar, tag in [(mpf(2), "Cbar=N=2  (M07-WINDOW2, MASTER:130; Lean)"),
                  (mpf(1), "Cbar=1    (C01-CBAR1:162 'physical', sum-level)"),
                  (mpf(20), "Cbar<=Delta=20 (A01-CBAR:74, provable UPPER bound)")]:
    b = mpf(8) / (CHAT * Cbar * EDELTA)
    print("  %-46s 8/(Chat*Cbar*eDelta) = %s   (ratio to 0.018394: %s)"
          % (tag, mp.nstr(b, 12), mp.nstr(b / BETA1, 8)))

hdr("[2] u at the window edge -- TWO conventions live in the corpus")
# (a) Z12-FIXC / P3-FULLREP bookkeeping: u = beta/8  =>  ubar = Chat*Cbar*eDelta*u
# (b) W95C-ORDER / T23 measured one:      u = I2/I1 = beta/4 at beta=1/48 (=1/192)
u_a = BETA1 / 8          # 2.29925e-3
u_b = BETA1 / 4          # 4.59850e-3  (= I2/I1 scaling of the 5.208239e-3 anchor)
print("  (a) u=beta/8 : u(beta1) = %s   ln(1/u) = %s   1/2 ln(1/u) = %s"
      % (mp.nstr(u_a, 12), mp.nstr(log(1 / u_a), 12), mp.nstr(log(1 / u_a) / 2, 12)))
print("  (b) u=beta/4 : u(beta1) = %s   ln(1/u) = %s   1/2 ln(1/u) = %s"
      % (mp.nstr(u_b, 12), mp.nstr(log(1 / u_b), 12), mp.nstr(log(1 / u_b) / 2, 12)))
print("  anchor check: 1/192 = %s ; measured u(1/48)=5.208239e-3 -> rel.dev = %s"
      % (mp.nstr(mpf(1) / 192, 12), mp.nstr(abs(mpf('5.208239e-3') - mpf(1) / 192) / mpf('5.208239e-3'), 4)))
print("  >>> 0.88291*5.208239e-3 = %s  (L3C-STARSTRONG:186 anchor)"
      % mp.nstr(mpf('0.88291') * mpf('5.208239e-3'), 12))
u_l3c = mpf('0.88291') * mpf('5.208239e-3')
print("      ln(1/that) = %s   half = %s  -> L3C's 5.382 / 2.691"
      % (mp.nstr(log(1 / u_l3c), 12), mp.nstr(log(1 / u_l3c) / 2, 12)))

hdr("[3] C = ln(Chat*Cbar*eDelta)  -- value under the three Cbar readings")
for Cbar, tag in [(mpf(1), "Cbar=1"), (mpf(2), "Cbar=2 (authoritative)"), (mpf(20), "Cbar=20 (upper bd)")]:
    C = log(CHAT * Cbar * EDELTA)
    print("  %-24s C = %s      C/2.691 = %s"
          % (tag, mp.nstr(C, 12), mp.nstr(C / log(1 / u_b) / 2, 8)))

C = log(CHAT * mpf(2) * EDELTA)
hdr("[4] structural identity at the window edge (u=beta/8 bookkeeping)")
print("  ln(1/u)|_{beta1} = ln(8/beta1) = %s" % mp.nstr(log(8 / BETA1), 12))
print("  C                = ln(Chat*Cbar*eDelta) = %s" % mp.nstr(C, 12))
print("  difference       = %s   -> C == ln(1/u) at the edge BY DEFINITION of the window"
      % mp.nstr(log(8 / BETA1) - C, 6))
print("  in the u=beta/4 bookkeeping (L3C's anchor): ln(1/u)=ln(4/beta1)= %s = C - ln2"
      % mp.nstr(log(4 / BETA1), 12))
print("  half of it = %s = C/2 - ln2/2" % mp.nstr(log(4 / BETA1) / 2, 12))

hdr("[5] the check:  C <= 2.691 ?")
thr = log(1 / u_b) / 2
for Cbar, tag in [(mpf(1), "Cbar=1"), (mpf(2), "Cbar=2"), (mpf(20), "Cbar=20")]:
    Cc = log(CHAT * Cbar * EDELTA)
    print("  %-10s C = %-16s  C <= 2.691 ?  %-5s  excess = %s x"
          % (tag, mp.nstr(Cc, 10), str(Cc <= thr), mp.nstr(Cc / thr, 8)))

hdr("[6] needed threshold / needed window")
print("  C <= 2.691  <=>  Chat*Cbar*eDelta <= exp(2.691) = %s   (now %s -> factor %s short)"
      % (mp.nstr(exp(thr), 10), mp.nstr(CHAT * 2 * EDELTA, 10),
         mp.nstr(CHAT * 2 * EDELTA / exp(thr), 8)))
print("  C <= 1/2 ln(1/u)  <=>  ln(8/beta) >= 2C  <=>  beta <= 8*exp(-2C) = %s"
      % mp.nstr(8 * exp(-2 * C), 8))
print("      = %s x (1/48) = %s x beta1" % (mp.nstr(8 * exp(-2 * C) * 48, 6),
                                            mp.nstr(8 * exp(-2 * C) / BETA1, 6)))
print("  [check] 1/2 ln(1/u) floor vs the honest floor ln(1/u)-C (both at beta1):")
print("      half-floor at beta1 = %s ;  ln(1/u)-C at beta1 = %s  -> the half-floor is"
      % (mp.nstr(log(8 / BETA1) / 2, 10), mp.nstr(log(8 / BETA1) - C, 10)))
print("      C/2 short of what is needed; the C-dependent floor is >0 for ALL beta<beta1.")

hdr("[7] digit trap: 2.691 (L3C) vs 3*pi^2/11 (one-loop SU(2) slope)")
print("  L3C 2.691  = 1/2 ln(4/beta1)   = %s" % mp.nstr(log(4 / BETA1) / 2, 12))
print("  3*pi^2/11  (C1-SLOPE-EQUIV:41, A6-TARGET:50) = %s" % mp.nstr(3 * pi**2 / 11, 12))
print("  difference = %s (rel %s) -> SAME 3 DIGITS, DIFFERENT OBJECTS"
      % (mp.nstr(3 * pi**2 / 11 - log(4 / BETA1) / 2, 6),
         mp.nstr((3 * pi**2 / 11) / (log(4 / BETA1) / 2) - 1, 4)))

hdr("[8] verdict")
print("  C := ln(Chat*Cbar*eDelta) = %s  (Cbar=2, authoritative)" % mp.nstr(C, 12))
print("  1/2 ln(1/u)|_{0.018394}  = %s  (L3C's 2.691, u=beta/4 anchor)" % mp.nstr(thr, 10))
print("  C <= 2.691 :  False,  by factor %s" % mp.nstr(C / thr, 6))
print("  threshold would need Chat*Cbar*eDelta <= %s (now %s)"
      % (mp.nstr(exp(thr), 8), mp.nstr(CHAT * 2 * EDELTA, 8)))
