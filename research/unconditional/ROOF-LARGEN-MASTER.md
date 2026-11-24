# ROOF-LARGEN-MASTER — the large-N master field / Makeenko–Migdal loop equation as the strut

> **VERDICT (no QED).** The large-N route is the **first strut that is renormalon-free by construction**: the `0⁺⁺` gap is a **physical mass (a spectral quantity)**, so it carries **no `u=2` IR renormalon** — it defuses *exactly* the obstruction that killed the condensate floor (`ROOF-CONDENSATE-FLOOR`). Large-N also **supplies a small parameter (`1/N`)**, addressing the "no small parameter at `α_s~O(1)`" wall. **But the leading term — the master-field `0⁺⁺` mass — is not computable:** the **4D Makeenko–Migdal loop equation does not close** (an infinite hierarchy; it closes only in 2D), and the **4D master field is not rigorously constructed** (Singer constructed it in 2D only). So the strut, in large-N form, is a *much better block* (renormalon-free, has a small parameter) than the condensate, but it is **blocked at computability/existence**: "close the 4D loop equation / construct the 4D master field." **Clay 0%.**

---

## §0 · Labels, object, acid

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** lattice evidence.
**Object:** `m_{0⁺⁺}(λ) := inf spec(H) − 0` in the `0⁺⁺` (trivial) sector, `λ = g²N` fixed as `N→∞` ('t Hooft limit). **Acid:** must fail on the massless free field (`m=0`) and be non-abelian-essential.

**Three iron laws carried** (from `SLOT-SYNTHESIS`): non-positivity / nonabelian / no-imported-small-parameter. Large-N's contribution is to the third: it *manufactures* a small parameter `1/N` rather than importing one.

---

## §1 · What large-N gives (the genuine assets)

**`[T]` Factorization + master field (Witten 1979; 't Hooft 1974).** In the 't Hooft limit, expectation values of gauge-invariant operators factorize (`⟨AB⟩ = ⟨A⟩⟨B⟩ + O(1/N²)`), so a single **master field** reproduces all correlators. This makes the `0⁺⁺` mass a **leading-order (`N⁰`) quantity** `m_{0⁺⁺}(λ)`, with a controlled `1/N²` expansion around it — a **small parameter exists**.

**`[T]` The gap is renormalon-free.** A physical mass is a **pole of a propagator / a spectral gap** — scheme-independent — and carries **no IR renormalon**. This is the decisive structural difference from `⟨tr F²⟩`: the condensate *is* the transseries coefficient of `e^{−2s}` (its ambiguity equals its own size, `ROOF-CONDENSATE-FLOOR`), whereas the **mass is a fixed spectral number**. In units of `√σ` (also a physical quantity) the ratio `m_{0⁺⁺}/√σ` is a **ratio of two renormalon-free observables** — no `u=2` ambiguity enters.

**`[V-num]` `m_{0⁺⁺}/√σ` is `N`-independent to leading order** (lattice; Athenodorou–Teper, Lucini–Teper–Wenger): `≈ 3.78 (SU2), 3.55 (SU3), 3.40 (SU4), … → ≈ 3.3 (SU∞)` — a `1/N²` approach to a constant (witness §C). The gap **survives to finite `N`** to within a few percent.

---

## §2 · Why this is a better block than the condensate

The roof round eliminated every known object; the two survivors were *both* blocked by the renormalon (`ROOF-CONDENSATE-FLOOR`) in the condensate coordinate. **In the mass coordinate the renormalon is absent.** So large-N does not just relabel `H1`: it changes the *type* of the unknown from "a renormalon-ambiguous condensate" to "**a renormalon-free physical mass**." That is real progress on the *spec* of the missing strut — it is the first formulation whose target is a clean number.

---

## §3 · `[NG]` The 4D loop equation does not close

The Makeenko–Migdal equation (1979), exact in the 't Hooft limit:
```
∂^μ_x (δ/δσ^{μν}(x)) ⟨W(C)⟩  =  λ ∫_C dy_ν δ^D(x−y) ⟨ W(C_{xy}) W(C_{yx}) ⟩   (+ self-intersection / splitting terms)
```
- **In D=2 it closes**: for a simple loop the RHS vanishes, the area-derivative reduces to `∂_A W = −(λ/2) W`, and the **exact master field** `⟨W(C)⟩ = exp(−λ·Area/2)` solves it (witness §A). This is the only dimension where the equation is closed and the master field is explicit.
- **In D=4 it does not close**: the RHS is a **double-contour** correlator `⟨W(C₁)W(C₂)⟩`, whose own equation brings a **triple**-contour term, and so on. The `0⁺⁺` mass sits at the bottom of an **infinite Makeenko–Migdal hierarchy** — the same infinity that made the DSE route (`SLOT-DSE-FIXPOINT`) die. **The leading gap is not computable from the loop equation alone** (witness §B).

⇒ The loop equation *defines* the strut but does not *evaluate* it. Computing the master-field `m_{0⁺⁺}` needs independent input — which is the wall relabelled.

---

## §4 · `[O]` The 4D master field is not rigorously constructed

Singer (1994) constructed the master field for **2D** Yang–Mills (where it is a free Gaussian, matching §3). **In 4D the existence of the master field is open** — expected (it is the consistent large-N limit implied by factorization), but **not constructed**. So even the *definition* of the leading-order `m_{0⁺⁺}(λ)` is not rigorous in 4D: we have lattice evidence and the factorization *ansatz*, not a master field.

⇒ The strut is blocked at **both** ends: the object's **existence** (4D master field, `[O]`) and its **evaluation** (4D loop-equation non-closure, `[NG]`).

---

## §5 · `[C]/[O]` Finite `N` and the crossover

- **Finite `N`:** `1/N²` corrections; lattice shows `m_{0⁺⁺}/√σ` nearly `N`-independent (witness §C) ⇒ the gap **survives** to `SU(3)`. No phase transition separates large-`N` from finite-`N` (the `1/N²` series is believed convergent/analytic); **rigorous control is absent** `[C]`.
- **Crossover:** at zero `T`, large-`N` YM is a **single confined phase for all `λ`** — the strong-weak **crossover is not a transition** (same fact as **B2**, `ROOF-TRANSFER`). So the large-`N` gap exists for all `λ` but has **no transition to anchor it**; the value at `λ_c` (the finite-`N` crossover coupling) is **not rigorously computable** `[O]`. The `1/N` expansion does not create a phase boundary, so it does not by itself supply the missing strut at `α_s~O(1)`.

---

## §6 · Verdict

| claim | label |
|---|---|
| `0⁺⁺` mass is renormalon-free (a physical spectral quantity) | **[T]** |
| large-N supplies a small parameter `1/N`; `m_{0⁺⁺}/√σ` is `N`-independent | **[T]/[V-num]** |
| 2D Makeenko–Migdal closes; exact master field `e^{−λA/2}` | **[T]** |
| 4D Makeenko–Migdal does **not** close (infinite hierarchy) ⇒ leading gap not computable | **[NG]** |
| 4D master field not rigorously constructed (Singer: 2D only) | **[O]** |
| finite-`N` survival (lattice `N`-independent) | **[C]** |
| crossover: single phase (B2), no anchor; `λ_c` value not rigorous | **[O]** |

**The strut the large-N route names, precisely:** *the master-field `0⁺⁺` mass `m_{0⁺⁺}(λ)` — a renormalon-free physical number — **provided** the 4D loop equation is closed (or the 4D master field constructed) and the `1/N` series is controlled at `N=3` and through the (non-transition) crossover.* This is a **cleaner target than the condensate** (no renormalon) and the **first with a small parameter**, but it is **blocked at computability/existence**, not at ambiguity. **Clay 0%; no QED.**

*Refs:* `ROOF-CONDENSATE-FLOOR.md`, `ROOF-CONFINEMENT-GAP.md`, `ROOF-TRANSFER.md` (B2), `SLOT-DSE-FIXPOINT.md`, `SLOT-SYNTHESIS.md`, `PROFILE.md`. Literature: 't Hooft (1974); Witten (1979); Makeenko–Migdal (1979); Singer, "On the master field in two dimensions" (1994); Athenodorou–Teper (large-N glueball spectrum); Lucini–Teper–Wenger (2004).
