# FALSIFY-OPE — the OPE-coefficients route does NOT falsify the only-exit claim: it is OPE-RHO (finite, blind) or H1 (infinite / renormalon)

> **VERDICT: the only-exit claim STANDS for this route.** Tried to build a *forcing* lower bound on the
> 0⁺⁺ atom from the `tr F²` OPE **coefficients** (operator-specific perturbative data) + exact Ward/anomaly,
> distinct from the previous OPE-as-data attack (`ROOF-OPE`). It fails, structurally, for one reason:
> **the OPE coefficients of the `1/Q^{2k}` terms ARE the linear moments `M_{-k}=∫s^{-k}dν`.** `Π(Q²)=∫dν/(Q²+s)`
> is **linear** in `ν`, so "operator-specific structure" adds no nonlinear constraint, and an **atomless bump
> at the floor reproduces every finite OPE moment** (`Δ_min=0`). The trace anomaly is one of these moments (a
> normalization), not a pole-forcing condition. No forcing exists; the only recoveries are the full expansion
> (`J→∞` = `H1`) and the perturbative coefficients' `u=2` renormalon (`ROOF-OPE` = `H1`).

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.

---

## §1 · The structural fact: OPE coefficients = linear moments

`Π(Q²) = ∫ dν(s)/(Q²+s) = Σ_{k≥0} (−1)^k M_{-(k+1)}/Q^{2k} + …` with `M_{-k}=∫s^{-k}dν`. So the OPE
coefficient of each `1/Q^{2k}` term is `M_{-k}`, a **linear functional of `ν`**. Therefore:
- the "operator-specific, perturbative-UV" structure lives in the **product** `C_n⟨O_n⟩`; to use it you need
  the **condensate** `⟨O_n⟩`, which is nonperturbative IR input — **not** free perturbative data;
- each Wilson-coefficient contribution to `Π(Q²)` is, after the dispersion relation, a function of `ν` that
  is **linear** — so it is in the same class as moments.

## §2 · Why it does not force (witness `w7/falsify_ope.py`, exact)

Reference: atom `f₁²=0.5` at `s₁=1` + continuum on `[2,4]` (mass 0.5). OPE moments `M_{-1..-8}` computed.
**Atomless replacement**: a uniform bump on `[s₁, s₁+2h]` of mass `0.5`, **zero point mass at `s₁`**, same
floor `inf supp = s₁`:

| h | max \|ΔM_{−k}\| over k=1..8 |
|---|---|
| 1e−1 | 2.43e−1 |
| 1e−2 | 3.77e−2 |
| 1e−3 | 3.98e−3 |
| 1e−4 | 4.00e−4 |

As `h→0` the **atomless** bump matches **all finite OPE moments** (to `O(h)`), with the **same floor** and
**zero atom weight**. ⇒ For any **finite** OPE truncation, an atomless measure is feasible ⇒ **OPE cannot
force `f₁²>0`**. The **trace anomaly** `A=M_{-2}` is just one of these moments — a normalization, matched by
the atomless bump too; it adds **no** pole-forcing condition.

## §3 · The three things the OPE could be (none forces)

- **Finite OPE coefficients** (any truncation): linear moments ⇒ **OPE-RHO class**, blind; atom hideable (above).
- **Full OPE** (all `k`, equivalently full `Π(Q²)` at every `Q²`): by Stieltjes inversion this **does** pin the
  atom — but that is **infinite data**, the `J→∞` limit = **`H1`** (the same wall in the OPE coordinate).
- **Perturbative coefficients themselves**: carry the `u=2` IR renormalon (`ROOF-OPE`), ambiguous `~Λ⁴` ⇒ **`H1`**.

**Operator-specific ≠ nonlinear.** The OPE of `tr F²` fixes *which* operators appear, but each contributes
**linearly** to `Π`; there is no sign-definite, condensate-independent structural constraint that forbids the
atomless solution. A positivity/reflection-positivity use of the identity-coefficient sign is an **upper** bound
(OPE-RHO class), not a pole-forcer.

## §4 · Honest bottom line

- **The only-exit claim STANDS for the OPE-coefficients route**: it reduces to **OPE-RHO** (finite, blind) or
  **`H1`** (infinite `J→∞` limit, or the `u=2` renormalon of the coefficients).
- **Why the falsification attempt was the right one to try**: OPE coefficients *are* computable
  perturbative-UV data (a genuinely different object from OPE-as-data / moments); the reason it still fails is
  **linearity of `Π` in `ν`** — the atom is hideable behind a continuum matching any finite linear data, which
  is exactly the `Δ_min=0` theorem.
- No new tool, no forcing. A genuine negative result that **strengthens** the only-exit claim by closing its
  most operator-specific escape. **Clay 0%; no QED.**

*Files:* `research/unconditional/FALSIFY-OPE.md`, `research/unconditional/w7/falsify_ope.py`.
*Refs:* `ROOF-OPE.md`, `TOOL-SYNTHESIS.md` (values→OPE-RHO, zeros/poles→A3), `PNT-OPE-RHO.md`.
