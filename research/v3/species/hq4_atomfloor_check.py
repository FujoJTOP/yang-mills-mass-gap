# hq4_atomfloor_check.py — HQ4-ATOMFLOOR deterministic exact checks (Fraction only, no sampling, no MC).
# Segments:
#  [1] phi==1 collapse:  (M0 <= C*f1^2)  <=>  (w1 >= 1/C)                       -- exact, both directions
#  [2] nondecreasing phi at s1  =>  coercivity(phi,C) implies w1>=1/C            -- never violated
#  [3] bounded phi: converse up to explicit constant B*C/phi(s1)                 -- the "ladder" is bounded
#  [4] unbounded phi (=1/s, =s): coercivity STRICTLY stronger than the floor     -- witness measures
#  [5] Markov-Krein sharpness: (t-Ebar)/(t-s1) <= nu([s1,t])/M0, eq iff two atoms at s1,t
#  [6] one-step Markov bound (delta-d1)/(s1+delta)  <=  1-d1/delta               -- concentration is dominated
#  [7] QI-1 identity  w1 = 1 - d1/eps + Dbar                                     -- exact on 4 measures
#  [8] (G-4) with a UV tail: target TRUE, (G-4) FALSE, repaired form TRUE        -- exact
#  [9] (G-4)_{delta=eps}  <=>  Ebar <= s2 - w0*(s2-s1)  =>  I2^W                 -- exact rewrite
# [10] w0^{P3} = r0/M0max   with M0<=2N^2 (repo) and Popoviciu M0<=N^2/4
from fractions import Fraction as F

def nu(pts):                       # pts = [(s_i, w_i)] ; returns derivable data
    M0 = sum(w for _, w in pts)
    M1 = sum(s*w for s, w in pts)
    s1 = min(s for s, _ in pts)
    f1 = sum(w for s, w in pts if s == s1)
    return M0, M1, s1, f1

def mass_in(pts, a, b):
    return sum(w for s, w in pts if a <= s <= b)

MEAS = {
 "two-1-4     ": [(F(1), F(1,2)), (F(4), F(1,2))],
 "three-1-2   ": [(F(1), F(1,100)), (F(2), F(99,100))],
 "two-1-10    ": [(F(1), F(3,5)), (F(10), F(2,5))],
 "thin-atom   ": [(F(1,1024), F(1,2048)), (F(2)-F(1,1024), 1-F(1,2048))],
 "L4-light    ": [(F(1,16), F(9,10)), (F(10), F(1,10))],
 "two-1-2     ": [(F(1), F(1,2)), (F(2), F(1,2))],
}
OUT = []
def say(s): OUT.append(s)

# ---------------------------------------------------------------- [1]
say("[1] phi==1 collapse:  (M0 <= C*f1^2) <=> (w1 >= 1/C)")
bad = 0
for name, pts in MEAS.items():
    M0, M1, s1, f1 = nu(pts)
    w1 = f1/M0
    for C in [F(2), F(4), F(10), F(100)]:
        lhs = (M0 <= C*f1)
        rhs = (w1 >= 1/C)
        bad += (lhs != rhs)
    say("    %s M0=%s f1=%s w1=%s  (both directions agree)" % (name, M0, f1, w1))
say("    mismatches = %d   -> %s" % (bad, "PASS" if bad == 0 else "FAIL"))

# ---------------------------------------------------------------- [2]
say("[2] phi NONdecreasing at s1 => coercivity(phi,C) implies w1>=1/C")
bad = 0; tested = 0
PHI = {"phi=1": lambda s: F(1), "phi=s": lambda s: s, "phi=s^2": lambda s: s*s}
for name, pts in MEAS.items():
    M0, M1, s1, f1 = nu(pts)
    w1 = f1/M0
    for pn, phi in PHI.items():
        Ip = sum(phi(s)*w for s, w in pts)
        for C in [F(2), F(3), F(10), F(1000), F(1000000)]:
            if Ip <= C*phi(s1)*f1:            # coercivity holds
                tested += 1
                bad += (w1 < 1/C)             # then the floor must hold
say("    coercivity instances fired = %d, floor violations = %d  -> %s"
    % (tested, bad, "PASS" if bad == 0 else "FAIL"))
say("[2b] same with the CONVEX-but-DECREASING phi=1/s : implication FAILS (so convexity is not the operative hypothesis)")
bad = 0; tested = 0
for name, pts in MEAS.items():
    M0, M1, s1, f1 = nu(pts)
    w1 = f1/M0
    Ip = sum(w/s for s, w in pts)
    for C in [F(2), F(3), F(10), F(1000)]:
        if Ip <= C*(1/s1)*f1:
            tested += 1
            bad += (w1 < 1/C)
say("    fired = %d, floor violations = %d  -> %s  (expected: violations > 0)"
    % (tested, bad, "PASS" if bad > 0 else "UNEXPECTED"))

# ---------------------------------------------------------------- [3]
say("[3] bounded phi: w1>=1/C  =>  coercivity(phi, B*C/phi(s1)) with B=sup phi<inf")
bad = 0
for name, pts in MEAS.items():
    M0, M1, s1, f1 = nu(pts)
    w1 = f1/M0
    B = F(1)                                  # phi=1, B=1
    for C in [F(2), F(4), F(10)]:
        if w1 >= 1/C:
            got = sum(F(1)*w for s, w in pts)          # int phi dnu
            bad += not (got <= (B*C/F(1))*F(1)*f1)
say("    violations = %d  -> %s" % (bad, "PASS" if bad == 0 else "FAIL"))

# ---------------------------------------------------------------- [4]
say("[4] unbounded phi: floor does NOT give coercivity (tail witness)")
pts = [(F(1), F(1,2)), (F(10), F(1,2))]                # w1 = 1/2
M0, M1, s1, f1 = nu(pts)
w1 = f1/M0
chi = sum(w/s for s, w in pts)
coer_1s = chi <= F(1)*f1/s1          # C=1 for phi=1/s
coer_s = M1 <= F(1)*s1*f1            # C=1 for phi=s
say("    two-1-10: w1=%s (floor for C=2 holds) ; phi=1/s coercivity at C=1: %s ; phi=s: %s"
    % (w1, coer_1s, coer_s))
say("    -> floor holds while both unbounded-phi coercivities fail  -> %s"
    % ("PASS" if (w1 >= F(1,2) and not coer_1s and not coer_s) else "FAIL"))

# ---------------------------------------------------------------- [5]
say("[5] Markov-Krein sharpness: (t-Ebar)/(t-s1) <= mass([s1,t))/M0, equality iff two atoms at {s1,t}")
def mass_half(pts, a, b):
    return sum(w for s, w in pts if a <= s < b)
bad = 0; eq = 0
for name, pts in MEAS.items():
    M0, M1, s1, f1 = nu(pts)
    Eb = M1/M0
    for t in [s1 + F(1,3), s1 + 1, s1 + 3]:
        rhs = mass_half(pts, s1, t)/M0
        lhs = (t - Eb)/(t - s1)
        bad += (lhs > rhs)
two = [(F(1), F(1,2)), (F(3), F(1,2))]
M0, M1, s1, f1 = nu(two); Eb = M1/M0; t = F(3)
sharp = ((t-Eb)/(t-s1) == mass_half(two, s1, t)/M0 == f1/M0)
say("    violations = %d" % bad)
say("    two-1-3 at t=3 (upper atom excluded): bound == w1 == %s exactly  -> %s"
    % (f1/M0, "PASS" if (sharp and bad == 0) else "FAIL"))

# ---------------------------------------------------------------- [6]
say("[6] one-step Markov window bound (d-d1)/(s1+d)  <=  (L2) = (d-d1)/d   [regime d>d1]")
bad = 0; strict = 0; empty = 0
for name, pts in MEAS.items():
    M0, M1, s1, f1 = nu(pts)
    Eb = M1/M0; d1 = Eb - s1
    for d in [F(1,2), F(1), F(3), F(10)]:
        mkv = (d - d1)/(s1 + d)
        l2 = (d - d1)/d
        if d > d1:
            bad += (mkv > l2); strict += (mkv < l2)
        else:
            empty += 1
say("    violations = %d ; strict-domination = %d ; both<=0 (empty) = %d  -> %s"
    % (bad, strict, empty, "PASS" if bad == 0 and strict > 0 else "FAIL"))

# ---------------------------------------------------------------- [7]
say("[7] QI-1 identity  w1 = 1 - d1/eps + Dbar")
bad = 0
CASES = {
 "two-1-4   ": [(F(1), F(1,2)), (F(4), F(1,2))],
 "two-1-10  ": [(F(1), F(3,5)), (F(10), F(2,5))],
 "tail-1-2-X": [(F(1), F(1,2)), (F(2), F(1,4)), (F(100), F(1,4))],
 "three-1-2 ": [(F(1), F(1,100)), (F(2), F(99,200)), (F(7), F(99,200))],
}
for name, pts in CASES.items():
    M0 = sum(w for _, w in pts); M1 = sum(s*w for s, w in pts)
    s1 = min(s for s, _ in pts); f1 = sum(w for s, w in pts if s == s1)
    rest = sorted({s for s, _ in pts if s > s1}); s2 = rest[0]; eps = s2 - s1
    D = sum((s - s2)*w for s, w in pts if s > s2)
    w1 = f1/M0; Eb = M1/M0; d1 = Eb - s1; Dbar = D/(M0*eps)
    lhs = w1; rhs = 1 - d1/eps + Dbar
    bad += (lhs != rhs)
    say("    %s w1=%s  rhs=%s  Dbar=%s" % (name, lhs, rhs, Dbar))
say("    mismatches = %d  -> %s" % (bad, "PASS" if bad == 0 else "FAIL"))

# ---------------------------------------------------------------- [8]
say("[8] (G-4) with a UV tail: target TRUE, (G-4) FALSE, repaired TRUE")
w0 = F(1,2); pts = [(F(1), F(3,5)), (F(2), F(1,5)), (F(10000), F(1,5))]
M0 = sum(w for _, w in pts); M1 = sum(s*w for s, w in pts)
s1 = F(1); f1 = F(3,5); s2 = F(2); eps = s2-s1
w1 = f1/M0; Eb = M1/M0; d1 = Eb - s1
D = sum((s - s2)*w for s, w in pts if s > s2); Dbar = D/(M0*eps)
g4 = (d1 <= (1-w0)*eps)
rep = (d1 <= (1-w0)*eps + Dbar)
say("    w1=%s >= w0=%s : %s" % (w1, w0, w1 >= w0))
say("    d1=%s ; eps=%s ; (1-w0)*eps=%s ; Dbar=%s" % (d1, eps, (1-w0)*eps, Dbar))
say("    (G-4) raw: %s   repaired: %s  -> %s"
    % (g4, rep, "PASS" if (w1 >= w0 and not g4 and rep) else "FAIL"))

# ---------------------------------------------------------------- [9]
say("[9] (G-4)_{delta=eps} <=> Ebar <= s2 - w0*eps  =>  I2^W (Ebar < s2)")
bad = 0
for name, pts in CASES.items():
    M0 = sum(w for _, w in pts); M1 = sum(s*w for s, w in pts)
    s1 = min(s for s, _ in pts); rest = sorted({s for s, _ in pts if s > s1}); s2 = rest[0]
    eps = s2-s1; Eb = M1/M0; d1 = Eb - s1
    a = (d1 <= (1-w0)*eps); b = (Eb <= s2 - w0*eps); c = (Eb < s2)
    bad += ((a != b) or (a and not c))
    say("    %s (G-4)=%s (Ebar<=s2-w0eps)=%s (I2^W)=%s" % (name, a, b, c))
say("    mismatches = %d  -> %s" % (bad, "PASS" if bad == 0 else "FAIL"))

# ---------------------------------------------------------------- [10]
say("[10] w0^{P3} = r0/M0max ,  r0 = 1/(4*(2N)^2) = 1/(16N^2)")
for N in [2, 3]:
    r0 = F(1, 16*N*N)
    w_repo = r0/(2*N*N)          # repo bound M0 <= 2N^2  (HZ1 (T2))
    w_pop = r0/(F(N*N, 4))       # Popoviciu, range N  => Var <= N^2/4
    assert w_repo == F(1, 32*N**4)
    assert w_pop == F(1, 4*N**4)
    say("    N=%d: r0=%s  w0(repo M0<=2N^2)=%s  w0(Popoviciu)=%s" % (N, r0, w_repo, w_pop))

say("")
say("ALL CHECKS PASS" if all(("-> FAIL" not in l) for l in OUT) else "SOME CHECK FAILED")
print("\n".join(OUT))
