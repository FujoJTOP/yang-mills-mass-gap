"""
INVENT-MODULAR witness  --  the modular / entanglement gap is a GAP-detector, not an ATOM-certifier.
Cheap, exact (<0.1 s), numpy only.  Illustrates the no-go in INVENT-MODULAR.md.
NOT a Yang-Mills calculation; NO QED.
"""
import numpy as np

# ---------------------------------------------------------------------------
# PART A | KMS.  For a thermal state rho = e^{-beta H}/Z on B(H), Tomita-Takesaki
#          gives modular operator Delta = e^{-beta H} (up to a central constant);
#          modular Hamiltonian K = -log Delta = beta H.  Hence the modular spectrum
#          IS the Boltzmann spectrum and the modular gap = beta*(energy gap).
# ---------------------------------------------------------------------------
print("=" * 74)
print("PART A | KMS: modular Hamiltonian K = beta*H  =>  modular gap = beta*(mass gap)")
print("=" * 74)
print(f"  {'sector':44s}{'beta':>6s}{'mod gap':>10s}{'entropy':>10s}")
for Delta, label in [
    (0.0, "free massless   (gap=0,  f1^2=0)"),
    (0.6, "gapped ATOMLESS (gap>0,  f1^2=0)  <-- :phi^2: proxy"),
    (2.0, "gapped glueball (gap>0,  f1^2>0)"),
]:
    for beta in (1.0, 2.5):
        E = np.array([0.0, Delta])
        K = beta * E
        modgap = K[1] - K[0]                    # = beta*Delta
        p = np.exp(-beta * E); p /= p.sum()
        S = -np.sum(p * np.log(p))
        print(f"  {label:44s}{beta:6.1f}{modgap:10.3f}{S:10.3f}")
print()
print("  => modular gap = beta*Delta is a FUNCTION of the mass gap: circular for the gap.")
print("  => ACID FAIL: the gapped-ATOMLESS sector (f1^2=0) has modular gap > 0, so the")
print("     modular gap fires where the atom weight is zero -- a gap-detector, NOT an")
print("     atom-certifier.  The free-field acid returns 0, but the discriminating step")
print("     (:phi^2:, f1^2=0, gap>0) is violated.")

# ---------------------------------------------------------------------------
# PART B | the modular spectrum is literally the Boltzmann spectrum -- full spectrum.
# ---------------------------------------------------------------------------
print()
print("=" * 74)
print("PART B | modular spectrum = Boltzmann spectrum (3-level), beta=2.5, E=(0,0.6,2.0)")
print("=" * 74)
beta = 2.5
E = np.array([0.0, 0.6, 2.0])
K = beta * E
lam = np.exp(-K); lam /= lam.sum()
print("  modular eigenvalues K_n = beta*E_n   :", np.round(K, 4))
print("  rho eigenvalues (normalised)         :", np.round(lam, 5))
print("  gaps K_n - K_0 = beta*Delta_n        :", np.round(K - K[0], 4))
print("  => the modular spectrum is exactly beta*(Hamiltonian spectrum): no new invariant.")

# ---------------------------------------------------------------------------
# PART C | sanity: a genuinely gapless spectrum gives no modular gap (acid side).
# ---------------------------------------------------------------------------
print()
print("=" * 74)
print("PART C | acid side: as the gap Delta -> 0 the modular gap beta*Delta -> 0")
print("=" * 74)
for Delta in (0.5, 0.1, 0.01, 0.0):
    print(f"  Delta={Delta:5.3f}  =>  modular gap = beta*Delta = {2.5*Delta:.4f}")
print()
print("Verdict (see INVENT-MODULAR.md): the modular spectrum is the transfer/thermal spectrum in")
print("disguise (KMS), or trivial (pure vacuum), or the entanglement spectrum (gap-detector).")
print("It is NOT a new handle for the 0++ ATOM.  Clay 0%; no QED.")
