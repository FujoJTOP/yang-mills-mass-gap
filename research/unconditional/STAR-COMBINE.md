# STAR-COMBINE — assembling and auditing the (★) contradiction chain, link by link

> **(★) = "a gapless interacting 4D Wilsonian IR attractor of an asymptotically-free gauge theory is conformal."**
> This file assembles the full contradiction chain and audits every link. **Verdict: it is a valid
> *conditional* proof `(0)+(2)+(3)+(4) ⇒ gapped`; the load-bearing links are (2) and (4b), both `[O]`. The
> load-bearing link (2) is genuinely a *different* statement from `H1` in **form** (RG attractors / conformality,
> not the spectral measure / atom / susceptibility), but of the **same difficulty** — it needs the identical
> nonperturbative strong-coupling IR control. So (★) is a clean re-framing of the wall, not a bypass; its value
> is a different technique pool.** Clay 0%; no QED.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go.

---

## §1 · The chain

- **(0)** Pure YM is asymptotically free: `β(g) = −β₀g³ − β₁g⁵ + …` with `β₀, β₁ > 0`.
- **(1)** AF ⇒ the IR is *not* asymptotically free (the coupling grows toward the IR).
- **(2)** A gapless interacting 4D relativistic vacuum ⇒ **CFT** (scale-invariant) **or a massless dilaton**.
- **(3)** Pure YM has **no massless dilaton**.
- **(4)** A CFT requires an IR fixed point `β(g*)=0`; pure YM has **none** (no Banks–Zaks zero — perturbative `[T]`; and no strongly-coupled fixed point — `[O]`).
- **(5)** ⇒ gapless is impossible ⇒ pure YM is gapped.

---

## §2 · Link-by-link audit

| link | statement | label | exact status |
|---|---|---|---|
| **(0)** | pure YM AF, `β₀,β₁>0` | **[T]** | 2-loop β rigorous at small `g`; coefficients positive for all `N` (witness: `β<0` ∀`g>0`, no zero). |
| **(1)** | AF ⇒ IR not free | **[C]/[O]** | the perturbative `β<0` is `[T]`, but "the IR never returns to `g=0`" needs "no nonperturbative IR fixed point" = **(4b)** — so (1) is partly circular with (4b); it is a *bridge*, not independent content. |
| **(2)** | gapless interacting 4D relativistic ⇒ CFT or dilaton | **[O] ← LOAD-BEARING** | the rigorous direction is **Polchinski `[T]`: scale-invariance + unitarity ⇒ conformal**. The **missing direction is "gapless ⇒ scale-invariance," which is folklore, not a theorem** — a gapless *non-scale-invariant* sector (no Goldstone, scale invariance not a symmetry at all) is not excluded a priori. The "or dilaton" branch is also not exhaustive. **This is the true gap.** |
| **(3)** | pure YM has no massless dilaton | **[C]** | solid: the dilatation current's divergence is the **trace anomaly** `θ=(β/2g)tr F²` (explicit, anomalous breaking) — no exact Goldstone; a spontaneously-broken-scale scenario gives a *massive* dilaton (`~` the glueball), not massless. |
| **(4a)** | no perturbative IR fixed point (no BZ) | **[T]** | `β=−(β₀g³+β₁g⁵)<0` ∀`g>0` for pure `SU(N)` (witness: N=2,3,4, no zero). BZ needs a conformal window with matter; pure YM has none. |
| **(4b)** | no *nonperturbative* IR fixed point | **[O]** | **the gap**: a strongly-coupled conformal IR fixed point at `g*~O(1)` is not excluded by the perturbative β. Proving "no strongly-coupled CFT IR fixed point for pure YM" is itself a frontier problem, equivalent-in-difficulty to the gap. |
| **(5)** | ⇒ gapped | **conditional** | valid only to the strength of (2) and (4b). |

---

## §3 · Is the load-bearing link `H1` in disguise, or genuinely new?

**In *form*: genuinely new.** `H1` is a statement about the **spectral measure ν** of `tr F²` (a susceptibility floor `χ ≤ 1/λ` / the atom `f₁²>0`). Link (2) is a statement about **RG attractors / conformality / scale invariance** — no measure, no atom, no `χ`, no weak-* discontinuity. It is a **different mathematical object**.

**In *difficulty*: the same wall.** Proving (2) rigorously for a non-supersymmetric 4D gauge theory requires exactly the **nonperturbative strong-coupling control of the IR** that `H1` also needs. So (★) does **not** bypass the wall; it is the wall seen from the **RG-attractor side** instead of the **spectral side**.

**Its genuine value:** the RG-attractor coordinate admits a **different technique pool** that the spectral route did not use —
- the **a-theorem** and monotonicity (a would-be gapless CFT IR must have `a_IR < a_UV`, and the *converse* — that an AF flow with no fixed point must gap — is the open content);
- the **conformal bootstrap** applied to the *hypothetical* strongly-coupled fixed point (does a unitary 4D CFT with pure-YM global charges and no relevant coupling exist? — this is the same kind of question the modern bootstrap answers for specific CFTs);
- **dilaton effective field theory** (the "walking vs dilaton" dichotomy from BSM model-building);
- **scale-vs-conformal** theorems (Polchinski and its extensions).

These are not the spectral/measure techniques, and (2) is the cleanest statement the campaign has produced that invites them.

---

## §4 · Honest bottom line

- **The chain is a valid *conditional* proof** `(0)+(2)+(3)+(4) ⇒ gapped`. A gap in (2) or (4b) breaks it.
- **Two load-bearing `[O]` links: (2) "gapless ⇒ CFT/dilaton" and (4b) "no strongly-coupled IR fixed point."** Both are frontier problems; both are the gap.
- **(★) is a clean re-framing, not a bypass:** different object (RG attractors) from `H1` (spectral measure), but the same difficulty (nonperturbative IR control), with a *different and untapped* technique pool.
- **The single most attackable sub-question:** make (2) precise — **"is a gapless interacting 4D relativistic vacuum necessarily scale-invariant?"** This is where the chain lives or dies, and it is a well-posed RG/conformality question, not a measure question.

*Refs:* `IMPOSSIBLE-CONTRADICTION.md` (where (★) emerged), `IMPOSSIBLE-SYNTHESIS.md`, `ONLY-EXIT-SYNTHESIS.md`; witness `w7/star_combine.py`.
