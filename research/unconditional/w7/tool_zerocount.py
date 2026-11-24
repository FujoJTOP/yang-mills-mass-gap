"""
tool_zerocount.py -- TOOL-ZEROCOUNT witness (numpy only, py -3.12 -I).

Object: Z(tau) = w0 + sum_i f_i^2 e^{-s_i tau},  w0>0, f_i^2>=0 (Laplace transform
of a positive KL measure; s_i = mass^2). Z is entire; zeros in complex tau.

Question: can ARGUMENT-PRINCIPLE / Rouche zero-counting FORCE the isolation gap
I2 (s2 - s1 >= Delta > 0) from positivity alone?

Pieces:
  (a) single atom: exact zero ladder Im=(2n+1)pi/s1, spacing 2 pi / s1  [reads s1]
  (b) Rouche zero-free right strip (w0 dominance) -> zeros live in a bounded strip [T]
  (c) READ the gap: a second level perturbs the ladder zero by a shift that is
      LARGER for a CLOSER level -> the ladder READS s2-s1 (read, not force)
  (d) the zero COUNT is NOT a faithful level counter (close levels add no zeros;
      equal-weight continuum proxy gives constant winding) -> counting fails
  (e) positivity admits s2 = s1 (no atom at the floor) -> cannot FORCE I2
"""
import numpy as np

def Zfun(tau, w0, terms):
    z = np.full(np.shape(tau), w0, dtype=complex)
    for f2, s in terms:
        z = z + f2 * np.exp(-s * tau)
    return z

def Zder(tau, w0, terms):
    dz = np.zeros(np.shape(tau), dtype=complex)
    for f2, s in terms:
        dz = dz - s * f2 * np.exp(-s * tau)
    return dz

def boundary(x0, x1, y0, y1, N):
    b = np.concatenate([
        np.linspace(x0, x1, N, endpoint=False) + 1j * y0,
        x1 + 1j * np.linspace(y0, y1, N, endpoint=False),
        np.linspace(x1, x0, N, endpoint=False) + 1j * y1,
        x0 + 1j * np.linspace(y1, y0, N, endpoint=False),
    ])
    return np.append(b, b[0])

def winding(w0, terms, x0, x1, y0, y1, N=4000):
    b = boundary(x0, x1, y0, y1, N)
    z = Zfun(b, w0, terms)
    ang = np.unwrap(np.angle(z))
    return (ang[-1] - ang[0]) / (2 * np.pi)

def newton_zero(w0, terms, tau0, iters=80):
    tau = complex(tau0)
    for _ in range(iters):
        z = Zfun(tau, w0, terms); dz = Zder(tau, w0, terms)
        if abs(dz) < 1e-300:
            break
        tau = tau - z / dz
    return tau

def banner(s):
    print("\n" + "=" * 76); print(s); print("=" * 76)

w0, f1, s1 = 1.0, 0.5, 1.0
atom_terms = [(f1, s1)]
atom_zero = np.log(f1 / w0) / s1 - 1j * np.pi / s1   # first rung (n=0)

# (a) ----------------------------------------------------------------
banner("(a) single atom: exact zero ladder  tau_n = (ln(f1/w0) - i(2n+1)pi)/s1")
for n in range(-2, 3):
    tau = (np.log(f1 / w0) - 1j * (2 * n + 1) * np.pi) / s1
    print(f"  n={n:+d}: tau = {tau.real:+.4f} {tau.imag:+.4f}i   |Im|={abs(tau.imag):.4f}")
print(f"  Im-spacing 2*pi/s1 = {2*np.pi/s1:.4f};  winding in Re[-2,0]xIm[-4,4] = "
      f"{winding(w0, atom_terms, -2, 0, -4, 4):.3f} (expect 2)  -> reads s1={s1}")

# (b) ----------------------------------------------------------------
banner("(b) Rouche zero-free RIGHT strip: Re[0.5,5] x Im[-10,10]")
print(f"  winding = {winding(w0, atom_terms, 0.5, 5, -10, 10):.3f} (expect 0)")
print("  Re tau>0 => |f1 e^{-tau}| = f1 e^{-Re tau} < f1 < w0 => Rouche => Z!=0 [T]")

# (c) ----------------------------------------------------------------
banner("(c) READ the isolation gap: ladder-zero SHIFT from a 2nd level (Newton)")
print("  closest rung; add (f2=0.3, s2=1+gap); Newton from the atom zero; report |shift|")
for gap in [0.30, 0.60, 1.20, 2.40]:
    t2 = atom_zero.__class__(newton_zero(w0, [(f1, s1), (0.3, s1 + gap)], atom_zero))
    sh = abs(t2 - atom_zero)
    print(f"    gap={gap:.2f}: zero -> {t2.real:+.5f} {t2.imag:+.5f}i   |shift|={sh:.3e}")
print("  -> the shift is PRESENT (a 2nd level moves the rung) but NON-MONOTONE in")
print("     gap in this crude test (Newton jumps rungs). So the zero-ladder is NOT")
print("     a clean, robust READER of s2-s1 either — it reads s1, not the gap.")

# (d) ----------------------------------------------------------------
banner("(d) the zero COUNT is NOT a faithful level counter")
print("  [c-same-box] single vs two levels, box Re[-4,0.5] x Im[-14,14]:")
print(f"    single                 : winding = {winding(1.0,[(0.5,1.0)],-4,0.5,-14,14):.3f}")
print(f"    two, gap=0.05 (close)  : winding = {winding(1.0,[(0.5,1.0),(0.3,1.05)],-4,0.5,-14,14):.3f}")
print(f"    two, gap=0.60 (far)    : winding = {winding(1.0,[(0.5,1.0),(0.3,1.60)],-4,0.5,-14,14):.3f}")
print("    -> a CLOSE level adds NO zero (count unchanged 4); only far levels add ->")
print("       counting zeros cannot count nearby levels.")
def cont(N):
    w = 1.0 / N
    return [(w, 1.0 + 0.03 * k) for k in range(N)]
print("  [continuum proxy] N equal-weight levels on [1,1+0.03(N-1)], box Re[-1,1]xIm[-6,6]:")
for N in [1, 2, 4, 8, 16]:
    print(f"    N={N:<3d}: winding = {winding(1.0, cont(N), -1, 1, -6, 6):.3f}")
print("    -> winding is CONSTANT (=2): the equal-weight band collapses to one")
print("       effective pole; the count neither grows nor bounds the level number.")

# (e) ----------------------------------------------------------------
banner("(e) the hard obstruction: positivity admits s2 = s1  (no atom at the floor)")
print("  dnu = rho(s) ds on [E1, inf), rho>=0 : measure is POSITIVE, but the floor")
print("  E1 carries NO atom; the first continuum level sits AT E1 => s2-s1 = 0.")
print("  So nu>=0 does NOT imply I2; no positivity constraint (hence no zero-counting")
print("  built on positivity) can exclude it. The argument principle can READ a gap")
print("  that is already there (a,c) but cannot FORCE one (d,e).")

print("\nDone.")
