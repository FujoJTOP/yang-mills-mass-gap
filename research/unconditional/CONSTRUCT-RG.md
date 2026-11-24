# CONSTRUCT-RG — the best rigorous 4D gauge RG assembled today reaches β<1/48 (IR gap) and β→∞ (UV control) but never the crossover; the exact missing estimate is the Savvidy-cone sign (★)=H1

> **VERDICT (no QED).** The most complete *rigorous* 4D pure-gauge RG constructible from the literature is a
> **two-ended flow**: exact Polchinski/Wetterich flow organised as the **Bałaban multiscale** medium-field
> expansion, with the **marginal log resummed by asymptotic freedom** (BBS's mechanism), the **IR gap**
> supplied by **Shen–Zhu–Zhu** at `β<1/48`, and the **large-β confinement** by **Cao–Nissim–Sheffield**
> (area law). It is rigorous on **both ends** — `β<1/48` (IR gap) and `β→∞` (UV stability) — and the plain
> **area-law** reaches `β~N`. But **no piece controls the sign of the effective-action Hessian at the
> crossover `β≈2.5`**, where the Savvidy cone is high-momentum and the `k→0` goal mode sits *inside* it.
> The best rigorous **gap** reach today is **`β<1/48≈0.0208`** (peer-reviewed; `β<0.0272` unrefereed). The
> exact missing estimate is one inequality — `Var_μ(a_{k→0}) ≤ 1/λ` — i.e. `(★)=H1`. **Clay 0%.**

---

## §0 · Labels

**[T]** 真定理 (unconditional in its stated domain) · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for a named class. Free-massless acid test: every step is checked to fail on the massless free field.

---

## §1 · The assembled RG map

The constructible object is the composition of four rigorous layers. Write the flow with `t=log(Λ_0/Λ)`.

**(1) Exact Polchinski flow** (the RG *equation*, rigorous, gauge-covariant regulator, modified Ward/Slavnov–Taylor):
```
∂_t V_t[A] = ½ ∫ ( δV_t/δA · C_t · δV_t/δA  −  δ²V_t/δA² · C_t ),
```
with `V_t` the effective action, `C_t = ∂_t(regulator covariance)`. On the BRST-restricted gauge-invariant operator space (`ghost-# 0`, dim 4) only `{Tr F², Tr F F̃}` survive; `Tr F F̃` is a total derivative and drops, so there is **one** coupling `g_j²` and `Z_j` — **no infinite counterterm tower** (contrast `φ⁴₄`, where the whole `λ_j`-tower must be tracked).

**(2) Bałaban multiscale medium-field expansion** (the RG *organisation*, `[T]`): split `A = A^{(s)}+A^{(l)}`; the **small field** is a Gaussian expansion about a background `A_j` (effective action with coupling renormalization and bounded effective potential); the **large field** is an R-operation (localisation + exponentiation + exponential suppression). This proves cutoff-uniform **UV stability** and convergence of the RG transformations (`β→∞`), in finite volume, background/axial gauge, `SU(2)` (Bałaban, *CMP* 95:17, 96:223, 98:17, 109:249, 116:1, 119:243, 122:175,355).

**(3) Marginal-log resummation by AF** (the BBS *mechanism*, ported): on the running trajectory `g_j² = 1/(2β₀ j)+O(j⁻²)`, so the marginal sum is absolutely convergent:
```
Σ_{j≥1} g_j⁴ = ζ(2)/(4β₀²) < ∞ ,   δg_J = O(J^{−3/2} log J) → 0 ,
```
(BBS, *LNM* 2242 (2019), marginal-log resummation; AF-analogue verified in `w7/construct_rg_witness.py` §A: `Σ_{j≤2e5} g_j⁴ = 190.69`, `ζ(2)/(4β₀²)=190.69`, rel. diff `3×10⁻⁶`). The BBS **method** — do **not** absolute-sum the BPHZ trees (false by the renormalon), prove `T_uni` as a **flow trapping region** — also ports:
```
‖R_{j+1}‖_X ≤ e^{−c} ‖R_j‖_X + C g_j⁴ ,   c,C independent of j,a,L .
```
This is the genuine BBS transplant: the renormalon is dodged by **tracking the flow, not summing the series** (C-BBS-MARGINALLOG §2).

**(4) SZZ IR gap at strong coupling** (the *gap* on the IR end): at `β<1/48` the lattice measure satisfies log-Sobolev + Poincaré uniformly in `L`, hence exponential clustering ⇒ volume-uniform **mass gap** (Shen–Zhu–Zhu, *CMP* 400:805 (2023), arXiv:2204.12737; `[T]`, all `d,N`, `λ=g²N>96`).

**The two-ended structure.** (1)–(3) control the action's *value/norm* uniformly up the UV tower `g_j→0`; (4) gives the gap at small `β`. They are both rigorous **where rigorous**, and they do not fight: this is the most complete honest object.

---

## §2 · Rigorous validity range — what holds where

| layer | rigorous statement | range | does it control the **sign** at β=2.5? |
|---|---|---|---|
| Polchinski flow | exact RG equation, gauge-covariant regulator, modified Ward | all scales (as an equation) | no — it *is* the equation; the sign is the unknown |
| Bałaban multiscale | cutoff-uniform UV stability, β→∞ | `β→∞` (UV), `SU(2)`, finite `V`, axial/background gauge | no — controls value/norm, **not** `δ²V`'s sign |
| AF log resummation | `Σ g_j⁴<∞`, `δg_J→0` | `g_j→0` (UV tower) | no — UV only |
| BBS trap | `‖R_{j+1}‖≤e^{−c}‖R_j‖+Cg_j⁴` | where the large-field Hessian is **convex** | **no** — YM's large-field Hessian is **indefinite** (`−D²+2g ad F_B`, `λ_min=−g|B|`), so `e^{−c}` becomes amplification on the Savvidy cone |
| SZZ IR gap | Poincaré/LSI ⇒ exponential clustering ⇒ gap | `β<1/48`, `L`-uniform, `a`-uniform | yes — **but only at `β<1/48`** |
| CNS area law | Wilson area law | `β<c_d` absolute, and into `β~N` (Hooft) | no — **area law, not the transfer-matrix gap**; 1-form channel |

The key structural fact (C-BBS-MARGINALLOG, ABC-CONSTRUCTIVE): BBS's one-step contraction survives into the large field **because** `φ⁴`'s large field is *more* convex (`−Δ+m²+12λφ² ⪰ m²`). YM's is *less* convex: `−D²+2g ad F_B ⪰ −2g|F_B|`. So `e^{−c}` (contraction) becomes **amplification** on `|k|<√(2g|F_B|)`. The trap cannot close there.

---

## §3 · The best β this program rigorously reaches today

- **Gap (peer-reviewed): `β<1/48 = 0.0208`** ⟺ `λ=g²N>96` (SZZ). This is the **only** proven 4D gauge uniform gap, and it is **strong coupling**.
- **Gap (unrefereed, 2026): `β<195/7168 ≈ 0.0272`** (gauge-orbit-curvature preprint — unrefereed).
- **Confinement/area-law (not the gap): `β<c_d` (absolute) and into `β~N`** (Cao–Nissim–Sheffield, arXiv:2505.16585, arXiv:2509.04688).
- **UV control (not the IR gap): `β→∞`** (Bałaban).
- **The crossover `β≈2.5`** (Hooft `λ≈7.3`, `β_W≈2.47`): **reached by none** of the above for the *gap*.

**Distance to the crossover: `λ 96 → 7.3 ≈ 13× in λ`** (witness §C). The barrier is an **analyticity radius** of the absolute methods, not a thermodynamic transition (`BASE-PHASE-RIGOROUS` §1). Constant-optimisation buys `<0.5` orders; it cannot close `~2` orders.

---

## §4 · The exact missing estimate to reach β=2.5

The RG program fails at exactly one place: the sign of the effective-action Hessian on the **typical-field, high-momentum cell**. There `Q≈Vκ_typ` so the coercivity multiplier `e^{Q/4g²}=O(1)`; the WSTC exponent `→1`; AF has `g≈1.26` (not small); and `N/12≈0.17` is short of `2g|F|` at the high-momentum cone that carries 90–97 % of the negative directions (W6-BE). The Hessian on the `k→0` mode is indefinite (`PᵀHP` min `≈−0.3`).

The exact identity (IBP, periodic) is
```
G″(0) = E_μ[V″|0] − Var_μ(V′|0),   V = S_W + Φ,
```
and the required inequality is
```
(★)   Var_μ(⟨v,∇S⟩|0) ≤ E_μ[H_vv|0] + W − δ
   ⟺ Var_μ(a_{k→0}) ≤ 1/λ
   ⟺ m_eff²(k→0) ≥ δ            =  H1 .
```
No exponential is available (`μ~g` ⇒ weight `O(1)`), and `W_geo≥0` (Prékopa) is a **drift, not a Hessian sign** and is anti-correlated with `|F|`. The sharpest universal bound is `G″ ≳ −0.08` (one variance bound short, `ABC2 §4`), so the cell's sign is the conditioned joint statement `Δ_dyn + X ≥ λ − W_geo` — a cross-covariance rearrangement inequality, not a positivity, and `SU(N)` has no lattice FKG.

**This is exactly the `H1` invariant, in the RG coordinate.** It is not a missing estimate "of size ε"; it is a missing **lower bound on a conditioned variance** — the irreducible core the whole campaign reduces to.

---

## §5 · Literature audit (numbers verified against the corpus, citations checked)

- **Bałaban** multiscale RG: *CMP* 95:17; 96:223; 98:17; 109:249; 116:1; 119:243; 122:175,355 — UV stability `β→∞`, `SU(2)`, finite `V`; **not** the IR gap. `[T]` (UV), `[NG]` (gap).
- **Bauerschmidt–Brydges–Slade**, *LNM* 2242 (2019): `φ⁴₄` trapping region + marginal-log resummation. **Ports**: the log (via AF) and the flow-not-sum method; **does not port**: the large-field convexity (indefinite YM Hessian). `[T]` (mechanism), `[O]` (gauge sign).
- **Shen–Zhu–Zhu**, *CMP* 400:805 (2023), arXiv:2204.12737: `β<1/48`, `L`-uniform gap. `[T]`.
- **Cao–Nissim–Sheffield**, arXiv:2505.16585 / 2509.04688: area law, absolute β and `β~N`; **gap not implied** (1-form vs 0-form). `[O]` at the bridge.
- **Chandra–Chevyrev–Hairer–Shen**: 2D *Publ. Math. IHÉS* 136:1 (arXiv:2006.04987); 3D *Invent. Math.* 237:541 (arXiv:2201.03487); Hairer framework *Invent. Math.* 198:269. 4D is **critical** (`[g]=0`), no gauge-covariant BPHZ ⇒ `[NG]` at 4D (see CONSTRUCTION-METHODS §(d)).
- **Magnen–Rivasseau–Sénéor**, *CMP* 155:325 (1993): UV-renormalized YM₄ with an **IR cutoff**, axial gauge, `SU(2)`, trivial sector. `[C]` (IR cutoff not removed), `[NG]` (gap).
- **Chatterjee**, *CMP* 366:203 (2019), arXiv:1502.07719; Chatterjee–Jafarov arXiv:1604.04777: strong-coupling large-N lattice; the weak-coupling master field is unbuilt, the constructed scaling limit is Gaussian-Proca (arXiv:2401.10507). `[NG]` (continuum gap).

*Correction logged:* `C-BBS-MARGINALLOG §1` writes `Σ_j g_j⁴ = ζ(4)/(16β₀⁴)`; for `g_j²=1/(2β₀ j)` the closed form is `ζ(2)/(4β₀²)` (verified numerically to `3×10⁻⁶`). The claim `Σ g_j⁴<∞` is unaffected; only the constant was off.

---

## §6 · Honest verdict

- **The best rigorous 4D gauge RG today is two-ended:** UV-stable to `β→∞` (Bałaban), gapped at `β<1/48` (SZZ), with the marginal log resummed (BBS/AF) and area law at large β (CNS).
- **It never reaches the crossover `β≈2.5` for the gap.** The best rigorous **gap** reach is `β<1/48 ≈ 0.0208` (peer-reviewed), an analyticity wall `~13×` in `λ` short of `β_W≈2.47`.
- **The exact missing estimate is one inequality, `(★) = Var_μ(a_{k→0}) ≤ 1/λ`** — the sign of the Savvidy-cone Hessian on the typical-field high-momentum cell — which is exactly `H1`.
- **Building the flow itself is done; closing its sign is the Clay step.** No QED. `Clay 0%`.

*Files:* `research/unconditional/CONSTRUCT-RG.md`; witness `research/unconditional/w7/construct_rg_witness.py` (runs clean, `py -3.12 -I`).
*Refs:* `ABC-CONSTRUCTIVE`, `ABC2-CONSTRUCTIVE`, `C-BBS-MARGINALLOG`, `CONSTRUCTION-METHODS`, `BASE-PHASE-RIGOROUS`, `BASE-INFINITE-VOL`, `ROOF-CONFINEMENT-GAP`.
