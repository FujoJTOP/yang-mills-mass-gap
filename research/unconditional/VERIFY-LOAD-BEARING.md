# VERIFY-LOAD-BEARING — auditing (i), (ii), (iii) of the conditional gap

**Verdict.** (i) the identities are CORRECT; one QUALIFICATION is under-stated (block→global localization) plus two idealizations. (ii) CORRECT and non-circular. (iii) an APPROXIMATION (exact iff the marginal is Gaussian), not an identity, with one internal inconsistency in its supporting identity. None is decisively WRONG or makes the theorem vacuous; but the theorem's weight sits openly on H1, which is the gap in curvature language.

## (i) Two-scale law of total variance — CORRECT, qualified
- The identity `Var_μ(F)=E_ν Var_{μ(·|U′)}(F)+Var_ν(E_μ[F|U′])` is EXACT and unconditional. G^E is compact ⇒ standard Borel; P continuous ⇒ ν=P_*μ and a regular conditional probability exist (disintegration). μ is a genuine probability measure. ✓
- `|∇F|²=|∇_∥F|²+|∇_⊥F|²` is exact **iff P is an orthogonal projector** in the bi-invariant metric; a gauge-covariant block transform's derivative generally is not. This P⊥Q tangent split is an IDEALIZATION, not flagged.
- `1/m=1/Λ_UV+1/λ`: CORRECT. The implication (fiber≤1/Λ, coarse≤1/λ ⇒ fine≤1/λ+1/Λ) holds because 1/λ,1/Λ ≤ 1/m. ✓
- **Load-bearing QUALIFICATION.** Step 2 bounds the *global* fiber term by (1/Λ_UV)E_μ[Σ|∇_⊥F|²] "by (H2)", but H2 is stated only for "local f on the block." Promoting block-local Poincaré to the global fiber needs a martingale/Lu–Yau localization. It is standard and available (PROVE-H2's Holley–Stroock is uniform in all boundary data), so it is bridgeable and volume-uniform — but it is never stated, and PROVE-H2 §scope(iii) miscalls the local→global step "open" (conflating the easy Poincaré-martingale with the hard Zegarlinski LSI-contraction). So "the conditional theorem's H2 obligation is discharged" is OVERSTATED: as written the theorem needs a global H2 it does not establish.
- H3's contraction `E_ν[|∇φ|²]≤E_μ[|∇_∥F|²]` is an ASSUMPTION (stated).
- No circularity in the identity. The theorem is conditional on H1, which is openly "the mass gap in curvature language." Not secret, not vacuous: it is a genuine one-RG-step reduction (fine gap ⟸ coarse gap + UV gap + faithfulness).

## (ii) Soft mode lives in coarse, frozen — CORRECT
- The softest physical mode is k_min=2π/L. The factor-2 block keeps |k|≤π/(2a); the coarse lattice (same box L, spacing 2a) carries momenta 2πn/L, n≤N_s/4, so k_min∈range(P)=coarse. Hence a is U′-measurable, ∇_⊥a=0, Var_fiber(a)=0: H2 is silent on the soft-mode AMPLITUDE; its gap is entirely H1. ✓
- So H2 alone ≠ gap, and H2 does not smuggle the soft-mode gap — this is exactly what keeps the theorem non-circular through H2.
- Nuance: the soft mode's *curvature* G″ does contain a fiber force-variance piece Var(∇S|a) (Lemma A bounds it via fiber Poincaré), but the decisive sign is the coarse Savvidy term = H1, not H2. "H2≠gap" holds; the decoupling is about the mode, correctly.

## (iii) G″(0)≈1/Var(a) — APPROXIMATION
- G″(0)=−(log P)″(0) is the peak curvature; 1/Var(a) is the global inverse second moment. For symmetric P equal **iff Gaussian**. Not a coincidence: quadratic Wilson action + quadratic Haar Jacobian (∏sinx/x≈e^{−x²/6}) and quartic≈0 give near-Gaussianity. Rigorous one-sided form: Brascamp–Lieb ⇒ Var(a)≤1/(min V_eff″), equality iff Gaussian.
- Numerics (W7 logs): representative modes 0.274 vs 0.268 (SU2) match ~2%, but ensemble G″ scatters SE≈0.05–0.07 about 0.27 (≈20%) and the mean G″≈0.31 EXCEEDS 1/Var≈0.27 by ~15%. SU(3) 0.485 vs 0.45 (~8%) but the branch-safe rerun gave **0/1500 valid configs** — fragile. So "≈2%" is best-case; honest ~2–20%.
- Use: PROVE-Gpp §5 converts H1 (curvature≥λ) to Var(a)≤1/λ heuristically (sign-level, H1⟺induced mass). NOT used in CONDITIONAL-LATTICE-GAP's proof chain (which uses H1 directly as SLC⇒Poincaré), so the approximation does not weaken the theorem.
- Secondary (W7, not affecting (iii)): W7 writes `V_eff″=E[H]−Var` "exactly … =G″(0)", but its own table has E[H]−Var=+0.02 vs G″=+0.274. PROVE-Gpp has the correct `G″=E[H]−Var+W` (W=Haar Jacobian, +0.25). W7's "exactly" omits W and contradicts its numbers.

## Bottom line
All three survive; none is vacuous or secretly circular. Honest gaps: (i) needs an explicit global fiber localization (bridgeable — standard martingale given uniform-in-boundary block Poincaré — but currently over-claimed as "discharged"); (i)'s orthogonal split is an idealization; H3 is assumed; (iii) is near-Gaussian approximation used only heuristically. The real weight is H1, openly equal to the gap.

*Sources:* CONDITIONAL-LATTICE-GAP.md; PROVE-H2-uvgap.md; PROVE-Gpp-positive.md; W7-Gpp0.md + `w7/` logs; W6-BE-at-2.5.md; VAR-FINAL.md, MARGINAL.md, VAR-TERM.md, MARGINAL-ENSEMBLE.md.
