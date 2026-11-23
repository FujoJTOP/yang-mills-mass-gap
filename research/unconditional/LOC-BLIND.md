# LOC-BLIND: fixed-cell coarse-mode localization cannot deliver the gap

**Verdict.** The proposed mechanism — coarse zero mode "localized", V̄(a) varying only O(1)
across a cell of size w ∼ ξ, and a conductance/Poincaré floor giving a uniform positive gap —
is **not a proof of `inf supp ν > 0`**. As stated it is *vacuous*: the free massless scalar
satisfies every premise and has no gap. The only regime where V̄ is sensitive to the mass is
w ≳ ξ = 1/m, which is circular. The one correct nearby theorem (local Poincaré + cross-cell
Dobrushin contraction ⇒ global Poincaré) is conditional on a strong-coupling hypothesis and,
at large β, *is* the coarse spectral gap restated. Separation below.

## 0. Objects.
Neutral scalar sector: gauge-invariant scalar φ, Källén–Lehmann measure ν with
⟨φ(x)φ(0)⟩_c = ∫ Δ_m(x) dν(m²), so the mass gap is inf supp ν. For a cell Λ_w of linear
size w, coarse zero mode a = |Λ_w|^{-1}∫_{Λ_w} φ, block measure ν_w(a), effective potential
V̄(a) = −log(dν_w/da) + const.

## 1. The w ∼ ξ vise (circularity). [Rigorous logic.]
ξ = 1/m is the target. Two regimes, both fatal to the proposal:
- **w ≪ ξ.** V̄ is fixed by G_m on |x−y| ≤ w ≪ 1/m, where G_m = G_0·(1 + O((mw)² ln(mw))).
  Every *local* functional of ν_w — V̄, its curvature, the local Poincaré constant — differs
  from its massless value by O((mw)² ln(mw)) → 0. The cell carries **no leading information
  about m**: it cannot certify m > 0.
- **w ≳ ξ.** Making V̄ depend on m at order 1 requires probing |x−y| ≳ 1/m, i.e. choosing w
  from ξ = 1/m — using the answer. At a continuum critical point ξ → ∞ and "w ∼ ξ" is the
  whole system; "fixed cell" self-destructs.

The criterion is sensitive to the gap *exactly when* it is circular.

## 2. Free-massless acid test. [PROVED, elementary.]
Take the 4D free massless scalar (gap 0; neutral-scalar sector, so in-sector). G_0(x) =
1/(4π²|x|²). Scaling, Var(a) = |Λ_w|^{-2}∫∫_{Λ_w}G_0 = c/w². Hence V̄(a) = a²/(2σ_w²),
σ_w² = c/w², and in units a/σ_w, V̄ = a²/2 *exactly*: O(1) variation across O(1) fluctuations
for every w. By §1 a very light massive theory is indistinguishable at this order. So
"V̄ varies only O(1) on the cell" **holds for the gapless theory** ⇒ cannot imply a gap.
Vacuity in the charter's sense.

Reverse direction also fails: a heavy mode (mw ≫ 1) has Var(a) = 1/(m²w⁴), curvature
1/Var = m²w⁴ — larger, yet the criterion value is flat across the whole window m ∈ [0, 1/w].
Variation does not track the gap.

## 3. Proof attempt: single cell ⇏ uniform L→∞. [PROVED counterexample.]
Gluing attempt via conditional variance for a cell Λ_1:
Var_ν(f) = E_ν[Var_{ν(·|σ(out))}(f)] + Var_ν(E_ν[f | σ(out)]).
Term 1 ≤ C_1 E[∫_{Λ_1}|∇f|²] — local Poincaré, always finite on a finite cell. **Term 2 is the
failure point:** it is the variance of the cell-mean as a functional of the *outside*, controlled
only if outside → cell influence contracts summably. For the free massless field every cell has
bounded C_1, yet the global Poincaré constant on Λ_L is ∼ L²/(4π²) (gap (2π/L)² → 0). Local
constants do not tensorize through interactions; Term 2 *is* the missing conductance. This is
exactly where the proof dies.

## 4. Steelman and where it turns circular. [CONDITIONAL.]
Add a uniform Dobrushin/conductance contraction: the conditional cell measure depends on its
neighbours with coefficient ϑ < 1 and dyadic influence is summable. Then (Zegarlinski-type)
global Poincaré ≤ C_1·(1 + Σ decays)/(1 − ϑ), uniform in L — a genuine gap, *conditional on*
ϑ < 1. But ϑ < 1 is a strong-coupling/high-temperature statement: it holds at small β and
fails at large β, where the interaction's concave direction makes cell ↔ neighbour dependence
non-contractive. At the Savvidy crossing the cross-cell conductance is the real bottleneck, and
bounding it *is* a lower bound on the coarse transfer-operator gap — the mass gap at the coarse
scale. The mechanism does not dissolve the wall; it relabels it.

## 5. Ledger.
- **Proved:** Gaussian blindness (local functionals insensitive to m for mw ≪ 1, §1); free-
  massless non-tensorization (single-cell Poincaré ⇏ uniform-in-L, §3); the w ∼ ξ vise (§1).
- **Refuted / no-go:** fixed-cell localization + O(1)-V̄ *alone* ⇒ gap (free massless satisfies
  the premises, violates the conclusion: vacuity).
- **Conditional (strong coupling):** localization + Dobrushin contraction ⇒ gap; correct, but
  reproduces the known small-β regime.
- **Open / exact missing step:** a uniform cross-cell conductance bound (ϑ < 1) at large β —
  equivalently a lower bound on the coarse transfer-operator gap across the Savvidy crossing.
  That object is the mass gap restated, so the proposal is circular at the coarse level.

**Answer to the colleague.** w ∼ ξ is circular; single-cell control does *not* extend to a
uniform L→∞ Poincaré constant (free massless is the witness); and a gapless theory realizes
O(1) variation while a heavy gapped one has large variation, so the implication is defeated in
the dangerous direction. The mechanism is a correct strong-coupling inequality wearing an IR
label.
