# KAPPA-2D-EXACT — 2D 对刚性退化,而 κ 是**非普适**的微观 LEC:它**不**像 Lüscher 项那样普适,**不能**从 2D / 薄模型转移到 4D ⇒ 4D 的 κ 需要微观芯 = `H1`

> **VERDICT (no QED).** The "单缺失输入 κ = H1" attack was run from the *universality* side: is the
> rigid-string rigidity κ a **universal** EFT coefficient (like the Lüscher term — computable once, then
> transferable), or a **microscopic** low-energy constant? **It is microscopic.** Three exact points:
> **(1)** 2D pure YM is **degenerate** for rigidity (`D−2 = 0` transverse directions ⇒ `κ` multiplies zero
> fields). **(2)** The Lüscher term is universal precisely because it is the zero-point energy of the `D−2`
> free transverse bosons (`−(D−2)π/(6L)` closed); the rigidity `(∂²X)²` is a 4-derivative LEC that is
> **removable by field redefinitions** in the exact Nambu–Goto action (Polchinski–Strominger) — NG generates
> **no** physical θ₂ rigidity. **(3)** A physical thick tube has `κ ~ σ·⟨r²⟩_⊥`, and `⟨r²⟩` is
> **profile-dependent** (at fixed nominal width it spans 12×, witness) ⇒ the O(1) coefficient is the
> microscopic core shape. **Therefore κ does NOT transfer from 2D / thin models; the 4D κ needs the
> microscopic core = `H1`.** Labels `[T]/[NG]/[O]/[V-num]`. **Clay 0%.**

---

## §0 · Labels / the question

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
Question: is κ universal (⇒ compute once in an exact theory ⇒ transfers to 4D) or microscopic (⇒ needs the 4D core = `H1`)?

## §1 · 2D pure YM is DEGENERATE for rigidity `[T]`

2D pure YM is exactly solvable — `W(C)=Σ_R d_R χ_R(U_C)e^{−g²C₂(R)A/2}`, an exact area law with
`σ_R = g²C₂(R)/2`. But the confining "flux tube" in `D=2` is the **entire worldsheet** — there is **no
transverse direction** (`D−2 = 0`). The rigidity `κ(∂²X^i)²` couples only to the `D−2` transverse
embedding fields; with `D−2=0` it is **identically zero / undefined**. The 2D theory has no transverse
profile, hence no κ to compute. (Witness §(2): `c=D−2=0` ⇒ trivial Lüscher term and no κ.) **2D cannot
host the rigidity.**

## §2 · Lüscher universal vs rigidity non-universal `[T]`

The **Lüscher term is universal** because it is the ζ-regularized **zero-point energy of the `D−2` free
transverse bosons** — it depends only on `c=D−2` and the worldsheet topology:
`E₀^closed = −(D−2)π/(6L)`, `E₀^open = −(D−2)π/(24L)` (witness §(1); for `D=4`: `−π/3` and `−π/12`, the 4×
ratio). This universality is the Aharony–Komargodski result (holds through `1/L⁵`).

The **rigidity `(∂²X)²` is the opposite**: it is a **4-derivative low-energy constant**. Polchinski–
Strominger (*Effective string theory*, PRL **67** (1991) 1681) showed that in the exact Nambu–Goto action
the `(∂²X)²` coefficient is **not determined by NG** and can be **removed by field redefinitions**; the
first **unambiguous** physical correction is at **6-derivative (θ₆)** order. So the "universal" (NG-generated)
part of κ is **zero/removable** — a physical κ can only come from the **microscopic thickness** of the tube.

## §3 · A thick tube: `κ ~ σ·⟨r²⟩_⊥`, profile-dependent `[V-num]`

The rigidity of a finite-width tube is set by the second moment of its transverse stress profile:
`κ ≃ c_profile · σ · ⟨r²⟩_⊥`. At a **fixed nominal width** `w`, `⟨r²⟩/w²` is **profile-dependent**
(witness §(3)):

| profile | `⟨r²⟩/w²` |
|---|---|
| Gaussian | 1.000 |
| flat disk | 0.500 |
| exponential | 6.000 |

a **12× spread** at the same `w`. And `⟨r²⟩` itself (the core size) is **not** fixed by the universal
log-slope of `w²(R)` (that slope is rigidity-free — CONSTRUCT-WIDTH-KAPPA-LATTICE) — it is the microscopic
`r₀`. So both the O(1) shape factor **and** the core size are microscopic: κ is **non-universal**.

## §4 · Transfer test: FAILS `[NG]`

Because κ is a non-universal 4-derivative LEC (§2) set by the microscopic core (§3), an "exact κ" computed
in 2D (§1: degenerate, none exists) or in any thin/exact model (NG: rigid part removable) does **not**
transfer to 4D. The Lüscher term transfers (universal, `D−2` only); **κ does not**. This is the sharp
**negative** answer to the "attack κ via universality" route: κ is *not* the kind of coefficient that can be
computed once and reused.

## §5 · Verdict

- **κ is NOT universal** — it is the `(∂²X)²` LEC, removable in NG, profile- and core-dependent for a thick
  tube. It is the **opposite** of the Lüscher term.
- **2D is degenerate** (`D−2=0`); no κ to compute there.
- **4D κ needs the microscopic core = `H1`** (in string coordinates). The universality attack **confirms**,
  rather than breaks, "κ = `H1`": the one input κ is genuinely microscopic and cannot be bootstrapped from
  exact low-dimensional solvability.

*Consistent with:* `CONSTRUCT-KAPPA-DERIVE` (§4: κ is a UV boundary condition), `CONSTRUCT-WORLDSHEET-CFT`
(κ is an irrelevant dim-4 LEC the CFT cannot fix), `CONSTRUCT-WIDTH-KAPPA-LATTICE` (width log-slope
rigidity-free, intercept degenerate with `r₀`). **No QED; Clay 0%.**

*Witness:* `w7/kappa_2d_exact.py` (runs `py -3.12 -I`, numpy only, <1 s): Lüscher coefficients for `D=3..6`,
2D degeneracy, and the profile second moments 1.0/0.5/6.0.
