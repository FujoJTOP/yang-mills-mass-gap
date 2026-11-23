# IN-11 (IN-b): can "lightest bulk state ∈ 0⁺⁺/2⁺⁺" be derived?

**Restatement.** Let `Δ_c` be the lightest zero-momentum bulk (gauge-invariant) mass in channel `c`.
IN-b: `Δ_{0⁺⁺} ∨ Δ_{2⁺⁺} = min_c Δ_c`. Weak form `(b′)`: no bulk state below the certificate `δ`
lies outside the C-even, P-even scalar/tensor sectors (`0⁻⁺`, odd C, spin ≥ 3 all `≥ δ`). Needed:
the defect bootstrap bounds only `0⁺⁺` (`tr F²`) and `2⁺⁺` (`T_{μν}`); a lighter state elsewhere
voids the transfer to `Δ`.

**True status.** Physical input, not symmetry-protected. *Provable (existence only):* trace anomaly
`T^μ_μ = (β/2g)tr F² + (Euler)` ties `tr F²` to spin-0; `T_{μν}` is the universal `2⁺⁺` operator —
both channels non-empty. *Not provable:* their ordering. Lattice (quenched SU(3), Morningstar–
Peardon `hep-lat/9901004`): `0⁺⁺` 1710 < `2⁺⁺` 2390 < `0⁻⁺` 2560 MeV — empirical, no Ward identity
fixes it. T1 killed the invertible symmetry/anomaly route.

**Fresh attacks — none derive it:**
- *Spectral/variational:* min–max yields upper bounds on excited energies; proving `0⁻⁺` heavy needs
  a **lower** bound in that channel = a gap statement there — same difficulty (cf. P4). Blocked.
- *Large-N:* factorization fixes `m ~ N⁰` for every single-trace dim-4 channel but orders none;
  `0⁺⁺`/`0⁻⁺` are both dim-4. No theorem; ordering is N-independent only empirically.
- *Anomaly with matter:* light quarks make the pion (`0⁻⁺`) lightest; chiral quarks make `0⁺⁺`/`0⁻⁺`
  parity partners, so `(b′)` is **false** with matter — it is dynamical, not protected.
- *Level-rank / non-invertible:* 4D non-invertible structure (continuous 2-category, Kaidi–
  Komargodski–Ohmori–Seifnashri–Shao, JHEP 2022) and 3D level-rank act on lines/surfaces/topological
  operators, not local glueball interpolators — no `J^{PC}` ordering.
- *Bootstrap bound on `0⁻⁺` (sign check):* the glueball S-matrix bootstrap (`2312.00127`) bounds
  **couplings**, not masses. And `q = tr F F̃ = ∂·K` gives, by Källén–Lehmann,
  `χ_t = Σ f_n² m_n² ≥ f₁² m₁²`, i.e. `m_{0⁻⁺}² ≤ χ_t/f₁²` — an **upper** bound. Topology/anomaly
  cannot lower-bound the pseudoscalar; the direction is wrong.

**Sharpest true statement.** The bootstrap certifies a lower bound only in the channels it sees; IN-b
is "the gap channel is visible," reducing to a channel-blind IR bound — equivalent to the mass gap.
**IN-11 is definitively physical** — neither derivable nor conditionally derivable from any known
symmetry, anomaly, large-N, or bootstrap structure.

**Flip condition.** IN-b becomes a theorem iff one proves `m_{0⁻⁺} ≥ m_{0⁺⁺}` (or any sub-`δ` lower
bound excluding every non-`0⁺⁺`/`2⁺⁺` channel). Or the conditional theorem is freed if `hcert`
(IN-10) is replaced by a channel-blind certificate needing no selection rule. Neither exists.
