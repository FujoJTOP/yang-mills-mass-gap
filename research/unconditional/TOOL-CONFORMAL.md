# TOOL-CONFORMAL — conformal geometry is BLIND to the atom weight (a point mass has zero capacity); the only geometric reader is the resolvent level-set radius = A3 (read-not-force)

> **VERDICT (no QED).** **Does not fill the tool blind spot.** Conformal mapping / extremal length /
> harmonic measure applied to the **spectral support** of the KL measure is **provably blind** to the atom
> weight `f₁²`: a point mass has **zero logarithmic capacity, zero harmonic measure, zero extremal length**,
> and two measures with identical support but different atom weight are **geometrically identical**. The
> **only** geometric reader of `f₁²` is the conformal "size" of a level set of the **resolvent** — which is
> exactly the Stieltjes residue `f₁²=−Res Π` = **A3 (read-not-force)**. So conformal geometry gives **no
> forcing, no non-blind tool**; it only supplies a **geometric avatar of A3**. `[NG]` for the class.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid test: must distinguish a true atom (`f₁²>0`) from the gapped-atomless `:φ²:` case (`f₁²=0`); else vacuous.

---

## §1 · Setup

KL measure `ν` on `[m₁²,∞)`, atom `f₁²=ν({m₁²})`, total mass `M=∫dν`. Resolvent (Stieltjes/Cauchy)
`Π(z)=∫ dν(s)/(z−s)`. The atom is a **simple pole** at `m₁²` with residue `f₁²=−Res Π`.

Solvable test measure: `ν = f₁² δ_{m₁} + (1−f₁²) δ_{m₂}`, `m₁<m₂`; `Π(z)=f₁²/(z−m₁)+(1−f₁²)/(z−m₂)`.

---

## §2 · Route A — support geometry is BLIND `[T]`

Every classical conformal invariant is a functional of the **support** only:
- **logarithmic capacity** `cap(E)`: for a single point `cap({p})=0`; `cap({a,b})=|a−b|/4`;
- **harmonic measure** `ω(z,{p})=0` for a point;
- **extremal length** of a point (degenerate set) gives no positive modulus.

Consequence: **the support `{m₁,m₂}` is fixed while `f₁²` varies over `(0,1]`, so every support-geometry
functional is CONSTANT in `f₁²`.** Witness: `cap=0.75` for all `f₁²∈{0.1,…,1.0}` while `−Res Π = f₁²`.
Two measures `0.9δ₁+0.1δ₄` and `0.1δ₁+0.9δ₄` have **identical** support geometry but atom weights differing
9× ⇒ **acid test FAILS** (support geometry gives the same value on `f₁²>0` and `f₁²=0`). This is the
**geometric form of the usc/lsc dichotomy**: the atom is a measure weight on a point; point-geometry assigns
a point **zero size**, so it cannot see a point's weight. (Same structural reason as `PNT-OPE-RHO`.)

---

## §3 · Route B — the resolvent level-set radius reads `f₁²`, but it is A3 `[T]`

The resolvent does encode `f₁²`: near `m₁`, `Π(z)≈f₁²/(z−m₁)`, so the **preimage of the level circle
`|Π|=R`** is a small contour around `m₁` of radius `≈ f₁²/R`. Witness: the preimage radius equals `f₁²/R` to
`0.05%` for `f₁²∈{0.1,…,1.0}` (`ratio 1.0000–1.0005`). So the **conformal "size" of a resolvent level set is
a geometric reader of `f₁²`** — but this is **exactly the residue** `f₁²=−Res Π` = **A3**. It is the same
"pole vs branch cut" information (`WI-MODIFIED-FUNCTIONAL`, `WI-ESCAPE`), in a conformal coordinate. It
**changes the representation, not the content**.

---

## §4 · Verdict and the exact step

- **Support geometry: blind** (point ⇒ zero capacity; acid-fails). **No forcing, no non-blind tool.**
- **Resolvent geometry: = A3** (the level-set radius is the residue; read-not-force).
- **The atom weight is not a geometric property of the support** — it is a measure-theoretic weight on a
  point, and every point-geometry invariant is zero on a point. **Conformal/extremal-length/harmonic-measure
  tools therefore cannot lower-bound `f₁²`**, for the same reason `∫ρφ`-integrals cannot (`PNT-OPE-RHO`).
- **Exact missing step `[O]`:** a forcing functional would require a geometric invariant **positive on a
  point mass and vanishing on a continuum** — which **does not exist**: a point has zero capacity, and any
  invariant that is positive there is reading the **residue** (= A3), i.e. the resolvent, not the geometry.
- **Genuine (small) yield:** the **conformal radius of a resolvent level set is a clean geometric avatar of
  A3** — a positive structural identification (A3 = "the preimage size of a level circle"), not a new lever.

---

## §5 · Witness

`w7/tool_conformal.py` (`py -3.12 -I`, numpy, runs clean): `(A)` cap constant `0.75` vs residue `f₁²`;
`(B)` preimage radius `=f₁²/R` (ratio `1.0000–1.0005`); acid test two measures same support / 9× weight
difference. Verdict line printed.

*Refs:* `PNT-OPE-RHO.md` (integral no-go, the same dichotomy), `WI-MODIFIED-FUNCTIONAL.md` /
`WI-ESCAPE.md` (A3 pole reading), `PROFILE.md` (usc atom / lsc data). Labels `[T]/[NG]/[O]/[V-num]`;
Clay 0%; no QED.
