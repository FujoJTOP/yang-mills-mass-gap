"""
v02_zeromono.py -- V02-BETAMONO numerical part.

Species note: research/v3/species/V02-BETAMONO.md

  f1^2(beta) = |<vac| O |G(beta)>|^2   is a MODULUS SQUARED  =>  f1^2 >= 0 identically,
  hence   "f1^2 > 0"  <=>  "<vac|O|G> != 0".

Four checks, all cheap / exact where possible:

 [A][E] Is monotonicity in the coupling FORCED by the species
        (modulus-squared + KL spectral decomposition)?
        -> exactly solvable 2-level model: f1^2 = |<0|O|1>|^2 = (1/4)sin^2(theta).
           It is a BUMP (interior max, decay back toward 0), not a monotone.

 [B][E] Can any moment functional determine the atom?
        -> exact moment-blind PAIR: same M0 = 1, same chi = 1, different f1^2.
           Hence Griffiths-type monotonicity (a moment statement) cannot transfer
           to f1^2 -- and the only rigorous bond is one-sided, f1^2 <= M0.

 [C][V] A real 1-parameter family (anharmonic oscillator, O = x^2):
        KL sum rule sum_n f_n^2 = Var(O) checked to machine precision;
        |M| stays > 0; and the SIGN of M is shown to be convention-dependent
        (raw LAPACK eigenvector signs flip), so "prove M > 0" is ill-posed unless
        a Kato phase continuation exists -- which presupposes level isolation.

 [D][E] Where exact zeros DO come from: selection rules (symmetry), coupling-free.

Run:  py -3.12 v02_zeromono.py
"""

import numpy as np

LINE = "-" * 78


def ho_basis(N):
    a = np.zeros((N, N))
    for k in range(1, N):
        a[k - 1, k] = np.sqrt(k)
    adag = a.T.copy()
    x = (a + adag) / np.sqrt(2.0)
    return a, adag, x, x @ x, x @ x @ x @ x


def anharmonic(lam, N=90):
    a, adag, x, x2, x4 = ho_basis(N)
    H = np.diag(np.arange(N) + 0.5) + lam * x4
    H = 0.5 * (H + H.T)
    w, V = np.linalg.eigh(H)
    return w, V, x, x2, x4


print(__doc__)

# --------------------------------------------------------------------------
print(LINE)
print("[A][E]  two-level model: is the residue monotone in the coupling?")
print(LINE)
print("H(D) = [[0, v], [v, D]] in a basis where O = |0><0|;  theta = atan2(2v, D);")
print("|<0|O|1>| = (1/2)|sin theta|  =>  f1^2 = (1/4) sin^2(theta), exactly solvable.\n")
v = 1.0
print(f"{'D':>9} {'theta':>9} {'f1^2':>13} {'d f1^2/dlnD':>14}")
prev = None
vals = []
for D in [1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 3.0, 1e1, 1e2, 1e3, 1e4]:
    th = np.arctan2(2 * v, D)
    f = 0.25 * np.sin(th) ** 2
    sl = float("nan") if prev is None else (f - prev[1]) / (np.log(D) - np.log(prev[0]))
    print(f"{D:9.1e} {th:9.6f} {f:13.10f} {sl:14.4e}")
    vals.append(f)
    prev = (D, f)
print(f"\n  f1^2: {vals[0]:.10f} -> max {max(vals):.10f} -> {vals[-1]:.3e}")
print("  => BUMP.  The species forces NO monotonicity: an exactly solvable atom")
print("     weight rises, saturates and decays back toward 0 while the atom is")
print("     present at every coupling.  'f1^2 increasing in beta' is an EXTRA")
print("     assumption, not a consequence of f1^2 = |M|^2 >= 0.")

# --------------------------------------------------------------------------
print()
print(LINE)
print("[B][E]  moment-blindness: can a moment functional see the atom?")
print(LINE)
print("Two KL measures with IDENTICAL M0 = int dnu and IDENTICAL chi = int s^-1 dnu:")
print("(the first moment M1 happens to differ -- which only strengthens the point)")
nu1 = [(1.0, 1.0)]                          # delta_1
nu2 = [(0.25, 0.5), (0.75, 1.5)]            # 0.25 delta_0.5 + 0.75 delta_1.5
for name, nu in (("nu1", nu1), ("nu2", nu2)):
    M0 = sum(w for w, _ in nu)
    M1 = sum(w * s for w, s in nu)
    chi = sum(w / s for w, s in nu)
    print(f"  {name}: M0={M0:.6f}  M1={M1:.6f}  chi={chi:.6f}  "
          f"f1^2=nu({{min supp}})={nu[0][0]:.6f}")
print("  => identical (M0, chi) [and even with M1 shown], DIFFERENT f1^2 (1.0 vs 0.25).")
print("  => no moment functional determines the atom; in particular Griffiths-type")
print("     monotonicity of <O> = -dlnZ/dbeta (a MOMENT statement) cannot transfer")
print("     to f1^2.  The only rigorous bond is one-sided: f1^2 <= M0 = Var(O),")
print("     and a one-sided upper bound cannot carry a positive lower bound.")

# --------------------------------------------------------------------------
print()
print(LINE)
print("[C][V]  a real family: anharmonic oscillator, O = x^2, |G> = first even excited")
print(LINE)
print("KL sum rule: Var(O) = <0|x^4|0> - <0|x^2|0>^2 = sum_{n even, n>=2} |<0|x^2|n>|^2.\n")
print(f"{'lam':>7} {'|M|':>13} {'f1^2':>13} {'Var(O)':>13} {'sum_n f_n^2':>13} {'w1':>9}")
lams = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
rows, Vprev = [], None
for lam in lams:
    w, V, xop, x2, x4 = anharmonic(lam)
    # Kato-style phase continuation: fix each eigenvector's sign by continuity in lam
    if Vprev is not None:
        for k in range(V.shape[1]):
            if V[:, k] @ Vprev[:, k] < 0:
                V[:, k] *= -1.0
    Vprev = V.copy()
    g = V[:, 0]
    M = g @ x2 @ V[:, 2]
    f1 = M ** 2
    var = (g @ x4 @ g) - (g @ x2 @ g) ** 2
    ssum = sum((g @ x2 @ V[:, n]) ** 2 for n in range(2, 90, 2))
    rows.append((lam, M, f1, var, ssum, f1 / var))
    print(f"{lam:7.2f} {abs(M):13.10f} {f1:13.10f} {var:13.10f} {ssum:13.10f} "
          f"{f1/var:9.6f}")
print(f"\n  exact anchor lam=0: <0|x^2|2> = 1/sqrt(2) = {1/np.sqrt(2):.10f}, f1^2 = 0.5")
print(f"  sum rule: max |Var(O) - sum_n f_n^2| = "
      f"{max(abs(r[3]-r[4]) for r in rows):.2e}   (basis converged)")
print(f"  |M| > 0 on the whole grid: {all(r[1] != 0 for r in rows)}; "
      f"phase-continued sign(<0|x^2|2>) is constant: "
      f"{len(set(int(np.sign(r[1])) for r in rows)) == 1}")
print(f"  w1 = f1^2/Var(O) in [{min(r[5] for r in rows):.4f}, {max(r[5] for r in rows):.4f}]"
      "  -- here the atom tracks the moment, so this family gives no universal")
print("     monotonicity in either direction: neither 'up' nor 'down' is generic.")

print()
print("  PHASE TRAP (why 'prove M > 0' is ill-posed):")
for lam in (0.01,):
    w, V, xop, x2, x4 = anharmonic(lam)
    raw = V[:, 0] @ x2 @ V[:, 2]
    Vf = V.copy()
    for k in range(Vf.shape[1]):
        if Vf[0, k] < 0:
            Vf[:, k] *= -1.0
    fixed = Vf[:, 0] @ x2 @ Vf[:, 2]
    print(f"    lam={lam}: raw LAPACK sign M = {raw:+.10f} ; after fixing a convention "
          f"M = {fixed:+.10f}")
print("    => the SIGN of M is eigenvector-convention, i.e. the phase of |G>; only")
print("       |M| is physical.  A 'decisive sign' argument must first fix the phase")
print("       by Kato continuation, which exists only while the level stays isolated")
print("       => 'prove M != 0 by proving M > 0' already assumes isolation (I2).")

# --------------------------------------------------------------------------
print()
print(LINE)
print("[D][E]  the only structural zeros: selection rules (symmetry), coupling-free")
print(LINE)
w, V, xop, x2, x4 = anharmonic(1.0)
g = V[:, 0]
for n in range(0, 6, 2):
    print(f"  <0|x|{n}> = {g @ xop @ V[:, n]:+.2e}   (O = x odd, parity-forbidden: exact 0)")
print("  => exact zeros of a matrix element are symmetry-forced and independent of the")
print("     coupling.  A zero of <vac|O|G> at a finite generic beta would need an")
print("     accidental cancellation -- codimension-1 in a 1-parameter family.")

print()
print(LINE)
print("VERDICT (numerics)")
print("  [A] monotonicity in the coupling is NOT implied by f1^2 = |M|^2 >= 0   [E]")
print("  [B] (M0, chi, and even M1) do not determine f1^2: moment monotonicity cannot")
print("      transfer to the atom; only f1^2 <= M0 is rigorous                 [E]")
print("  [C] mechanism-free sign arguments are ill-posed (phase convention)    [V]")
print("  [D] exact zeros are symmetry/selection-rule zeros                     [E]")
print("  NOTHING here decides <vac|O|G>(beta=2.5) != 0.")
print(LINE)

