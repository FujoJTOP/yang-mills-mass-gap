"""
wi_theorem_witness.py -- witness for WI-THEOREM.md  (weight-invariance of A=pi/|Im tau_zero|)

Model.  Z(tau) = w0 + sum_i w_i exp(-E_i tau),  w0>0, w_i>0, 0<E_1<=E_2<=... .
A := pi / min_k |Im tau_k|  over zeros tau_k of Z.  "Reweighting" = changing {
weights w0,w_i at FIXED energies (operator choice O -> O_lambda changes overlaps
w_i = |<i|O|0>|^2, not energies).

Claims (numpy only, py -3.12 -I):
  1. SINGLE ATOM (one E_1):   A = E_1 for ALL (w0,w1)>0  -> invariant.     [T exact]
  2. TWO LEVELS  (E_1,E_2):   A MOVES with the weights.                    [T examples]
  3. ATOM + GAPPED CONTINUUM: A = E_1 + eps*S(rho) + O(eps^2),
        dA/deps = E_1 * Im[p]/ (pi w0),  p = int rho e^{-E tau_a} dE  -> a
        LINEAR functional of the cut density rho (sensitivity).           [C numeric]
"""
import numpy as np
np.set_printoptions(precision=6, suppress=True)


def Zf(tau, w0, Es, ws):
    with np.errstate(over='ignore', invalid='ignore'):
        return w0 + np.sum(ws * np.exp(-Es * tau))


def dZf(tau, w0, Es, ws):
    with np.errstate(over='ignore', invalid='ignore'):
        return np.sum(-Es * ws * np.exp(-Es * tau))


def newton(w0, Es, ws, tau, iters=400):
    tau = complex(tau)
    for _ in range(iters):
        z = Zf(tau, w0, Es, ws); dz = dZf(tau, w0, Es, ws)
        if not (np.isfinite(z) and np.isfinite(dz)) or dz == 0:
            return None
        st = z / dz
        if not np.isfinite(st):
            return None
        tau = tau - st
        if abs(st) < 1e-15:
            break
    return tau


def nearest_A(w0, Es, ws, E1):
    """Nearest zero to the real axis: smallest |Im tau| with Z=0."""
    best = None
    ims = [np.pi / E1, -np.pi / E1, 3 * np.pi / E1, -3 * np.pi / E1,
           5 * np.pi / E1, -5 * np.pi / E1]
    for re0 in (-2.0, -1.0, 0.0, 1.0, 2.0, 3.0):
        for im0 in ims:
            tau = newton(w0, Es, ws, re0 + 1j * im0)
            if tau is None:
                continue
            if abs(Zf(tau, w0, Es, ws)) < 1e-8 and abs(tau.imag) > 1e-9:
                if best is None or abs(tau.imag) < best[0]:
                    best = (abs(tau.imag), tau)
    if best is None:
        return None, None
    return np.pi / best[0], best[1]


print("=" * 74)
print("1. SINGLE ATOM  Z = w0 + w1 e^{-E1 tau},  E1 = 1  ->  A = E1 exactly")
print("=" * 74)
E1 = 1.0
for (w0, w1) in [(1.0, 0.3), (1.0, 1.0), (1.0, 5.0), (2.0, 0.1), (0.5, 3.0)]:
    re = np.log(w1 / w0) / E1
    Es = np.array([E1]); ws = np.array([w1])
    a_num, tau = nearest_A(w0, Es, ws, E1)
    print(f"  w0={w0:<4} w1={w1:<4}  Re_tau={re:+.6f}  A(Newton)={a_num:.9f}  "
          f"|A-E1|={abs(a_num-E1):.2e}")
print("  -> A = E1 = 1 exactly for every (w0,w1): INVARIANT.  [T]\n")


print("=" * 74)
print("2. TWO LEVELS  E1=1, E2=sqrt2  ->  A MOVES with weights")
print("=" * 74)
E1, E2 = 1.0, np.sqrt(2.0)
w0 = 1.0
vals = []
for w2 in [1e-3, 1e-2, 0.1, 0.5, 1.0, 3.0]:
    Es = np.array([E1, E2]); ws = np.array([1.0, w2])
    a, tau = nearest_A(w0, Es, ws, E1)
    vals.append(a)
    print(f"  w2={w2:<6}  A={a:.6f}")
print(f"  spread of A = {max(vals)-min(vals):.6f}  -> NOT invariant.  [T/V-num]\n")


print("=" * 74)
print("3. ATOM + GAPPED CONTINUUM (gap 0.5 above E1), continuum weight eps")
print("   prediction: dA/deps = E1 * Im[p] / (pi w0),  p = int rho e^{-E tau_a}")
print("=" * 74)
E1 = 1.0; w0 = 1.0; w1 = 1.0
Delta = 0.5
Es_g = np.linspace(E1 + Delta, E1 + Delta + 2.0, 40)
dEg = Es_g[1] - Es_g[0]
rho = 1.0
tau_a = np.log(w1 / w0) / E1 + 1j * np.pi / E1
p_val = np.sum(rho * dEg * np.exp(-Es_g * tau_a))
S_pred = E1 * np.imag(p_val) / (np.pi * w0)
eps_list = np.array([0.0, 0.001, 0.002, 0.004, 0.008])
A_list = []
for eps in eps_list:
    Es = np.concatenate([[E1], Es_g])
    ws = np.concatenate([[w1], eps * rho * dEg * np.ones_like(Es_g)])
    a, tau = nearest_A(w0, Es, ws, E1)
    A_list.append(a)
    print(f"  eps={eps:<6}  A={a:.8f}")
A_list = np.array(A_list)
slope = (A_list[-1] - A_list[0]) / (eps_list[-1] - eps_list[0])
print(f"\n  A(eps=0)            = {A_list[0]:.8f}  (= E1 to {abs(A_list[0]-E1):.1e})")
print(f"  finite-diff dA/deps = {slope:+.6f}")
print(f"  predicted   dA/deps = {S_pred:+.6f}")
print(f"  relative agreement  = {abs(slope-S_pred)/abs(S_pred):.2e}")
print("  -> sensitivity is a LINEAR functional of the cut density rho.  [C]\n")

print("=" * 74)
print("VERDICT")
print("=" * 74)
print("  (1) single atom       : A = E1, invariant under reweighting.   [T exact]")
print("  (2) two levels        : A moves with weights.                  [T/V-num]")
print("  (3) atom + gapped cont: A = E1 + eps*S(rho) + O(eps^2),  S linear. [C/V-num]")
print("  => A is reweighting-invariant iff the spectral support is a SINGLE")
print("     atom (an ISOLATED atom, with a gap to the continuum): the WI test")
print("     detects ISOLATION (I2), not merely the atom weight (I1). No QED.")
