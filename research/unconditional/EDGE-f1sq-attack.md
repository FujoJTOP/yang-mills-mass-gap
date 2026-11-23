# EDGE-f₁²-attack — can a YM-specific NONLINEAR proven-sign inequality lower-bound the edge residue `f₁²=ν({m₁²})`?

> 2026-11-23 · Direct hit on the one open object (LEVER-nonconvex-edge; CONVERGENCE-redteam §(b)3 identifies this route as genuinely open — only a *class* is excluded, not every route). **Verdict: not constructed.** All five nonlinear constructions collapse. The obstruction is sharper than "no linear lever": it is a **sign** obstruction. No fabricated proof.

**Setup.** `ν≥0` (OS reflection positivity), `supp ν⊂[m₁²,∞)`, gauge-invariant interpolator; `f₁²=ν({m₁²})`, total moments `M_k=∫sᵏdν`. Available YM data: anomaly coefficients (UV/total), Ward/BRST identities (linear/cohomological), RP (PSD), strong-coupling lattice (β<β₀), GZ (gauge-dependent).

**The unifying obstruction (sign).** Every nonlinear functional descending from OS positivity is **Schur-concave in the spreading/majorization order**: `log det` is concave in the Loewner order (Hadamard/Fischer), entropies are concave, variance constraints are one-sided. An edge atom `cδ_{m₁²}+ρ` is majorization-**minimal**. Hence any concave constraint `Φ(ν)≤C` is satisfied *most easily* by concentrated measures and can never forbid or force the atom; and any `Φ(ν)≥C` with concave Φ is maximized by spreading, away from the atom. To *force* `f₁²>0` one needs a **Schur-convex lower barrier**, i.e. an identity of the shape

  `[∫φ dν] ≤ C(φ)·ν({m₁²})`  — a *total bounded above by the atom* (coercivity).

That is the exact reverse of every Cauchy–Schwarz/RP inequality, which always runs *atom-implied ⇒ total larger*. **This reversed sign is the missing identity.** No proven YM structure carries it.

**Five attempts, exact death:**
1. *Nonlinear Ward/STI at k=0.* Anomaly gives `χ_t ≥ f₁²m²` (product **upper**, CF IN-11). STIs are linear/cohomological; the transverse vertex has no anomalous slot (ghost#1 dim-5 empty, W3-marginal-BPHZ) ⇒ no protected residue. Missing: an STI bounding a protected total by `f₁²`.
2. *Brascamp–Lieb / log-Sobolev on a second variation.* Hypothesis is action-convexity — killed by the Savvidy negative mode (T5, W7: G″(+) but k→0 unreachable; `E[H]` ~99% cancelled by force variance). The OS second variation is variance positivity = Cauchy–Schwarz = **upper**. Wrong sign + false hypothesis.
3. *Lattice sum rule, proven sign.* Proven only at β<β₀≈1/48 (SZZ), far below β_c≈2.47; even there any edge-residue bound ≡ a second-level `Λ` bound (residue–level-gap lemma). "Proven sign" ≠ resonance dominance, the only thing pinning an operator's overlap; no uniform-in-a control.
4. *Determinant on physical subspace.* GZ FP operator `−∂·D≥0` on the Gribov region is again PSD ⇒ det bounded **above** (Hadamard). `D(0)=0` is gauge-dependent and does **not** imply `f₁²>0`.
5. *'t Hooft anomaly / topological charge.* Anomalies are index/total statements (moments); `q=∂·K` gives `χ_t` upper on `f₁²m²`; matching pins contact terms / θ=π vacuum degeneracy (GKKS), never a dynamical 0⁺⁺ atom. Insensitivity lemma (LEVER-(B)).

**Is `D(0)=0` (GZ) provable?** No, and it would not help. It concerns the gauge-dependent gluon propagator, not the gauge-invariant residue (no implication `D(0)=0⇒f₁²>0`). The horizon condition is a variational *scenario* built on a local renormalizable action (Zwanziger), not derived from the YM measure; exact DSE does not close; and lattice shows `D` IR-suppressed but *finite*, so `D(0)=0` is not even the empirical statement.

**Sharpest line.** The nonlinear route is genuinely open as a *class* (LEVER excluded only positivity/convex/moment data) — but the specific object demands a **Schur-convex coercivity** (total ≤ atom), the one sign absent from anomaly, Ward, RP, lattice and GZ alike. Constructing it *is* the edge theorem, hence the gap. Still true: the free field (g→0) is an admissible RP theory with `f₁²=0` (edge at 0); so any bound must inject nonperturbative IR dynamics — the gap, renamed.
