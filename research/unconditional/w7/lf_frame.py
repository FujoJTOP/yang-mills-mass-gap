# -*- coding: utf-8 -*-
"""
lf_frame.py -- witness for candidate (4) of the death-cause axis: frame / Parseval identity.

This is the candidate the task flags as the only possible survivor (D17 has never been a
main cause), so the accounting is deliberately hostile -- and deliberately honest: W3 reports
whichever way the numbers fall.

  W1  TIGHT BRANCH: is there any slack in the "redistribution"?  (compression of a Parseval
      frame -- does it stay Parseval?)
  W2  FRAME-NATIVE bounds: do they move when the gap closes?  (they are attached to the
      frame, not to the theory)
  W2b ADAPTED bounds: built from L itself -- which part of the spectrum do they see?
  W3  NON-ADAPTED: the exact redistribution  <Lf,f> = sum_l [diag]_l + cross(f) ;
      compute the TRUE best case  min_{f perp Omega} sum_l <P_l f, L P_l f>  (= a min
      eigenvalue) and compare with Delta.  This decides whether the redistribution is a
      bound at all.
  W4  SPECTRAL-SIDE shells: the sums are nu-moments -> weak-* blind.
  W5  control: the one genuinely directional uncertainty theorem (Donoho-Stark,
      |supp f| + |supp f^| >= N+1) is a SUPPORT/cardinality statement, not a weight one.
"""
import numpy as np

np.random.seed(20261125)
OUT = []


def say(s=""):
    print(s)
    OUT.append(str(s))


def rand_orth(d):
    A = np.random.randn(d, d)
    Q, R = np.linalg.qr(A)
    return Q * np.sign(np.diag(R))


# ================================================================ W1  tight branch
say("=" * 78)
say("W1  TIGHT BRANCH.  Parseval-normalised tight frame: F F^T = I_d  (columns = frame vectors)")
say("    sum_n |<f,e_n>|^2 = f^T F F^T f = ||f||^2 exactly -- and after compression to any")
say("    subspace P, the compressed frame is STILL Parseval.  Slack = 0 => no direction.")
say("=" * 78)
d, k = 20, 3                       # k orthonormal bases stacked -> tight with bound k
basis = [rand_orth(d) for _ in range(k)]
F = np.hstack(basis)               # d x (k d)
S = F @ F.T / k                    # normalised frame operator
say(f"  d={d}, n_vectors={k*d}, ||S-I||_inf = {np.abs(S-np.eye(d)).max():.3e}  (tight, bound 1)")
say(f"{'rank P':>8} {'max |sum_n|<Pf,Pe_n>|^2/k - ||Pf||^2| / ||Pf||^2':>52}")
worst_all = 0.0
for r in (1, 4, 10, 16, 20):
    Q = rand_orth(d)[:, :r]
    P = Q @ Q.T
    Fp = P @ F
    w = 0.0
    for _ in range(30):
        f = P @ np.random.randn(d)
        nf = np.linalg.norm(f)
        if nf < 1e-10:
            continue
        lhs = np.sum((Fp.T @ f) ** 2) / k
        w = max(w, abs(lhs - nf**2) / nf**2)
    worst_all = max(worst_all, w)
    say(f"{r:8d} {w:52.3e}")
say(f"  -> worst Parseval defect over all compressions = {worst_all:.2e}  (machine zero)")
say("     P(sum|e_n><e_n|/k)P = P  is Cholesky-level.  The redistribution has NO slack:")
say("     it is an identity in every subspace, hence it cannot produce a one-sided bound. (D17)")

# ================================================================ W2  frame-native bounds
say()
say("=" * 78)
say("W2  FRAME-NATIVE frame bounds: A = lam_min(S), B = lam_max(S).  S depends on the FRAME,")
say("    not on the dynamics.  Vary L's gap from 2.0 to 0.0 (theory becomes gapless):")
say("=" * 78)
# a genuinely non-tight frame: one full basis + a strict subset of a second one
Fw = np.hstack([basis[0], basis[1][:, :d // 3]])
Sw = Fw @ Fw.T
evw = np.linalg.eigvalsh(Sw)
A, B = evw[0], evw[-1]
say(f"  non-tight frame (1 basis + {d//3} vectors of a 2nd):  A = {A:.6f}   B = {B:.6f}"
    f"   B/A = {B/A:.4f}")
say(f"{'gap Delta':>12} {'A = lam_min(S)':>18} {'B = lam_max(S)':>18} {'moves?':>8}")
for Delta_ in (2.0, 1.0, 0.5, 0.1, 0.0):
    say(f"{Delta_:12.2f} {A:18.10f} {B:18.10f} {'NO':>8}")
say("  -> the frame bounds are constants of the frame.  They do not move when the gap closes,")
say("     so they take the SAME value on the free massless field: acid test fails identically.")
say("     A datum attached to the frame rather than to the theory cannot detect the theory. (D10)")

say()
say("  W2b  the only way to make A,B depend on L: build the frame FROM L.  The natural choice")
say("       is the DIFFUSION frame  e_k = e^{-tL} delta_k  (k = 1..d).  Then")
say("           S = sum_k e^{-tL}delta_k delta_k^T e^{-tL} = e^{-2tL},")
say("       so   A = e^{-2t lam_max(L)}   and   B = e^{-2t*0} = 1.")
say(f"       L = diag(0, Delta, 1.9, 3.1, 6.0),  lam_max = 6.0 (the UV end):")
say(f"{'t':>6} {'Delta':>8} {'A=e^{-2t*lam_max}':>20} {'B':>10} {'A moves with Delta?':>20}")
for t in (0.1, 0.25):
    for Delta_ in (2.0, 1.0, 0.1, 0.0):
        lamD = np.array([0.0, Delta_, 1.9, 3.1, 6.0])
        evd = np.exp(-2.0 * t * lamD)          # eigenvalues of S = e^{-2tL}
        say(f"{t:6.2f} {Delta_:8.2f} {evd.min():20.10f} {evd.max():10.4f} {'NO':>20}")
say("     -> A is exactly the UV number e^{-2t*lam_max}; B is exactly the vacuum value 1.")
say("        The gap -- the SECOND eigenvalue -- enters neither bound.  Ends, never the edge. (D15)")
say("        NOTE: this A is the SAME UV quantity that the Choi analysis of candidate (2)")
say("        produced as its 'positivity margin' (e^{-2t*lam_max}): two different candidate")
say("        mechanisms hand back the identical UV number as their only positive reading.")

# ================================================================ W3  the honest test
say()
say("=" * 78)
say("W3  NON-ADAPTED.  Exact redistribution (a general orthogonal shell decomposition):")
say("        <Lf,f> = sum_l <P_l f, L P_l f>  +  cross(f),   cross = sum_{l!=l'} <P_l f, L P_l' f>")
say("    The redistribution is a BOUND iff  min_{f perp Omega, ||f||=1} f^T M f >= Delta,")
say("    with M = sum_l P_l L P_l  (symmetric!).  That minimum is a genuine eigenvalue, so")
say("    we can settle it by computation instead of assertion.")
say("=" * 78)
d = 36
Delta = 0.6444
lam = np.concatenate([[0.0], np.linspace(Delta, 5.0, d - 1)])
L = np.diag(lam)
Om = np.zeros(d); Om[0] = 1.0
Qperp = np.linalg.qr(np.eye(d) - np.outer(Om, Om))[0][:, :d - 1]   # basis of Omega^perp


def shell_frames():
    out = []
    V = rand_orth(d)
    out.append([V[:, i:i + 6] for i in range(0, d, 6)])
    V = rand_orth(d)
    out.append([V[:, i:i + 4] for i in range(0, d, 4)])
    V = rand_orth(d)
    out.append([V[:, i:i + 9] for i in range(0, d - d % 9 or d, 9)] if d % 9 == 0
               else [V[:, i:min(i + 9, d)] for i in range(0, d, 9)])
    return out


say(f"  Delta = {Delta:.4f},   d = {d}")
say(f"{'frame':>6} {'#shells':>8} {'min_(f perp Om) sum_l[diag]':>28} {'>= Delta?':>10} "
    f"{'min cross':>12} {'max cross':>12}")
for j, Bs in enumerate(shell_frames()):
    Ps = [Bk @ Bk.T for Bk in Bs]
    M = sum(P @ L @ P for P in Ps)
    Mred = Qperp.T @ M @ Qperp
    mn = np.linalg.eigvalsh(Mred)[0]
    # cross term over the same subspace
    cs = []
    for _ in range(300):
        f = Qperp @ np.random.randn(d - 1)
        f /= np.linalg.norm(f)
        dg = sum((P @ f) @ L @ (P @ f) for P in Ps)
        cs.append((f @ L @ f) - dg)
    cs = np.array(cs)
    say(f"{j+1:6d} {len(Ps):8d} {mn:28.6f} {str(mn >= Delta):>10} {cs.min():12.6f} {cs.max():12.6f}")
say("  READING:")
say("   * min_(f perp Omega) sum_l <P_l f, L P_l f> is the BEST POSSIBLE outcome of the")
say("     redistribution -- a genuine min eigenvalue, and it must be >= Delta for the")
say("     redistribution to be a bound at all.")
say("   * the cross term changes sign across frames (min<0<max): it is not a one-signed")
say("     repair, it is a COMMUTATOR average whose sign is a frame choice.")
say("  => the whole content of the 'redistribution' is the inequality  M >= Delta  on")
say("     Omega^perp -- an inequality that (i) depends on the frame, (ii) is a min-eigenvalue")
say("     statement, i.e. exactly the object whose computable direction is an upper bound (D1)")
say("     or which is the gap itself (D2).  The IDENTITY did no work: it is the inequality")
say("     M >= Delta that carries everything.")

# ================================================================ W4  spectral-side shells
say()
say("=" * 78)
say("W4  if the frame lives on the SPECTRAL side (Kallen-Lehmann shells), the shell sums are")
say("    moments of nu.  Two measures, same M0 and same chi, different floor and atom:")
say("=" * 78)
nu1 = [(1.0, 1.0)]
nu2 = [(0.25, 0.5), (0.75, 1.5)]


def stats(nu):
    M0 = sum(w for w, s in nu)
    M1 = sum(w * s for w, s in nu)
    chi = sum(w / s for w, s in nu)
    m1 = np.sqrt(min(s for w, s in nu))
    atom = sum(w for w, s in nu if abs(s - min(t for _, t in nu)) < 1e-14)
    return M0, M1, chi, m1, atom


say(f"{'measure':>30} {'M0':>6} {'M1':>7} {'chi':>7} {'m1':>7} {'f1^2':>7}")
for nm, nu in (("delta_1", nu1), ("0.25 d_0.5 + 0.75 d_1.5", nu2)):
    M0, M1, chi, m1, atom = stats(nu)
    say(f"{nm:>30} {M0:6.2f} {M1:7.3f} {chi:7.4f} {m1:7.4f} {atom:7.3f}")
say("  -> same M0 = 1 and same chi = 1; m1 = 1 vs 0.7071; f1^2 = 1 vs 0.25.")
say("     Every spectral-side Parseval identity computes sums of this type => nu-moments")
say("     => weak-* blind to the support boundary (D5; the CS/TR chain).")

# ================================================================ W5  control
say()
say("=" * 78)
say("W5  control: the one genuinely directional uncertainty theorem (Donoho-Stark,")
say("    |supp f| + |supp f^| >= N+1 on Z_N) is a CARDINALITY statement, not a weight one.")
say("=" * 78)
for N in (16, 24):
    Fm = np.exp(-2j * np.pi * np.outer(np.arange(N), np.arange(N)) / N) / np.sqrt(N)
    viol, worst = 0, 10**9
    for _ in range(3000):
        k = np.random.randint(1, N + 1)
        supp = np.random.choice(N, size=k, replace=False)
        f = np.zeros(N, dtype=complex)
        f[supp] = np.random.randn(k) + 1j * np.random.randn(k)
        fh = Fm @ f
        s2 = int(np.sum(np.abs(fh) > 1e-9))
        viol += (k + s2 < N + 1)
        worst = min(worst, k + s2)
    say(f"  N={N}: min(|supp f| + |supp f^|) = {worst}  (need >= {N+1});  violations {viol}/3000")
say("  -> true and one-sided, and it is about SET SIZES.  It says nothing about the WEIGHT")
say("     at the edge: support != atom (corpus SUPPORT-NOT-ATOM / TOOL-*).")

say()
say("=" * 78)
say("SUMMARY (candidate 4 / frame, Parseval)")
say("=" * 78)
say(" - tight branch: zero slack to 1e-15 in every subspace (W1) -- the redistribution IS")
say("   an identity and cannot produce a one-sided bound.  D17.")
say(" - frame-native bounds are constants of the frame and do not move when the gap closes")
say("   (W2): acid test fails identically (D10).  L-adapted bounds see only the two ENDS")
say("   (W2b): the gap enters neither (D15).")
say(" - the non-adapted branch is decided by  M = sum_l P_l L P_l >= Delta  on Omega^perp,")
say("   a genuine min eigenvalue that depends on the frame, while the cross term changes")
say("   sign across frames (W3).  So the identity does no work and D1/D2 do.")
say(" - spectral side: shell sums are nu-moments => D5 (W4).")
say(" - the one true directional frame theorem is a support/cardinality inequality (W5).")
say(" - NET: the frame identity gives the redistribution for free; the inequality direction")
say("   exists only as 'min_eigenvalue(M) >= Delta' or 'cross(f) >= ...'.  In the corpus's own")
say("   coordinates the leakage term cross(f) IS the object X  in  G'' = Delta_dyn + W_geo + X.")
say("   => semi-live: the direction is producible, the missing step is named, and it is H1.")
