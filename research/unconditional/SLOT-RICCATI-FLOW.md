# SLOT-RICCATI-FLOW — the `k→0` sign-definite self-consistency in the flow/Riccati coordinate: the exact Wetterich flow projects to a **non-closing Riccati** in `κ(k)`; `κ≡0` is its free/Gaussian solution and the `κ²` term is **marginal**; a positive IR fixed point requires a **sign-definite source at `κ=0`** — absent at the coarse mode at β≈2.5 (Savvidy negative mode is transverse-Landau-gapped & high-momentum). Exact failing step = *no sign-definite source at `κ=0`* = `H1`. **Not a new species; the flow-coordinate dialect of `H1`.** Clay 0%.

> 2026-11-24 · ROUTE: Riccati / running-mass self-consistency (task slot: PROFILE §3 P4). Attempt to write the `k→0` sign-definite flow for the `0⁺⁺` mass `κ(k)` that **excludes the trivial (massless) trajectory**, or pin its exact failure.
> **判词 (verdict up front).** The exact Wetterich/Polchinski flow of the coarse marginal (E1-FRG, `[T]`) projects at the symmetric vacuum `a=0` to an **exact Riccati-in-scale for `κ(k):=V_{eff,k}''(0)`** — eq. **(R)** below — that I derive and **verify numerically** (instantaneous identity, rel. err `5.5·10⁻⁴`, `[V-here]`). Its structure is decisive and negative:
> - **(i) `κ≡0` is a flow solution** — exactly the free/Gaussian trajectory (source `C(k)≡0`), acid-clean; **[T]**.
> - **(ii) the `κ²` term is marginal.** Alone it lets `κ` decay only logarithmically toward `k→0` (`[V-here]`, witness: `κ:1→0.505` over `k∈[0.02,1]`) and **can never hold a nonzero IR fixed point**. A **positive** fixed point requires the **linear coefficient and the source to have a definite sign at `κ=0`** — i.e. an IR-*relevant* direction. That is exactly MS-GAP-EQUATION's `F′(0)>1` / the `k→0` tachyon.
> - **(iii) no sign-definite source at `κ=0` exists at the coarse mode at β≈2.5.** The source `C=−β₀φ₃−σ₀λ+β₃−σ₄` and the linear coefficient `B=3σ₂−2β₁` are **sign-indefinite**, controlled by the coarse self-energy `Π(0)`; `Π(0)≥λ ⟺ Var(a;κ)≤1/λ ⟺ H1`. The Savvidy negative mode that would supply `C>0` is **transverse-Landau-gapped and 90–97 % high-momentum** at 2.5 (W6-BE), so the coarse block is positive and the source keeps the **gapless** sign. **[V-doc/T]**.
> - **(iv) the Riccati does not close** — `φ₃=V'''(0)`, `λ=V''''(0)`, `∂²σ`, `∂²β` carry their own flows (E1-FRG §4); closing it is a truncation = the DSE vertex Ansatz (PROOF-FRG-DSE §2). **[V-here]**.
> - **Second face (residue flow):** the multiplicative transport `f₁²=Z_O²f₁²` is a **linear** flow `∂_k f₁²=2γ_O f₁²` — sign-preserving (`Z_O²>0`), but `0` is an **absorbing** fixed point; it cannot manufacture the atom. **[T] identity / [O] at crossover**.
>
> **The route reduces to Schema C (gap equation `κ=F(κ)`) + Schema A (Banach–Riccati) — the flow-coordinate dialect of `H1`.** It is neither a new species nor a closure. **Exact failing step: `κ=0` is the (marginal) flow solution and no sign-definite source at `κ=0` exists without `Π(0)≥λ` = `H1`.** Scorecard `2/5`; **Clay 0 %; no QED.**

*Labels*（宪章 §2）: **[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** 具名类 no-go · **[V-here]** 本轮新算 · **[V-doc]** 本会话核实原文 · **[V-num]** 数值见证（仅结构，非物理）· **[R]** 外部据报。**不在未证上敲 QED。**

---

## §1 · The exact flow and its `a=0` projection (the Riccati)

### 1.1 The flow (`[T]`, E1-FRG §1)
Let `μ_k` be the scale-`k` exact-RG measure, `a(A)=⟨v,\log U⟩` the coarse mode, `p_k(a)` its marginal, `V_{eff,k}=−\log p_k`. The **exact** drift–diffusion (Fokker–Planck) flow `∂_k p_k=−∂_a(β_k p_k)+∂_a²(σ_k p_k)` gives the viscous Hamilton–Jacobi equation `[T]` (E1-FRG eq. §1):
```
  ∂_k V_k(a) = −σ_k V_k''  +  σ_k (V_k')²  +  (2σ_k' − β_k) V_k'  +  (β_k' − σ_k'')      (*)
```
with `β_k(a)=E_k[−ȧ|a]` (drift), `σ_k(a)=½E_k[⟨ȧ,ȧ⟩_c|a]≥0` (conditional half-variance = the one-loop-per-scale kernel, built from `(Γ_k^{(2)}+R_k)^{-1}(k∂_kR_k)(Γ_k^{(2)}+R_k)^{-1}`), primes `=∂_a`. `(* )` is `∂_k`, integrating from UV `k=Λ` down to IR `k→0`.

### 1.2 The Riccati (`[V-here]`, derived + numerically verified)
Define the `0⁺⁺`/coarse-curvature **running mass**
```
   κ(k) := V_{eff,k}''(0) = m_eff²(k; k→0)     ( = G″  of the corpus;  κ=m₁² at k=0 ).
```
At a `Z₂`-symmetric vacuum `a=0` (`V_k'(0)=0`), differentiate `(*)` twice in `a` and evaluate at `0`. The algebra (checked symbolically, then against the PDE — below) yields the **Riccati-in-scale**:
```
┌───────────────────────────────────────────────────────────────────────────────┐
│  ∂_k κ  =  2σ₀ κ²  +  (3σ₂ − 2β₁) κ  −  β₀ φ₃  −  σ₀ λ  +  (β₃ − σ₄)        (R) │
│  A = 2σ₀ ≥ 0 ,   B = 3σ₂ − 2β₁ ,   C(k) = −β₀φ₃ − σ₀λ + β₃ − σ₄              │
└───────────────────────────────────────────────────────────────────────────────┘
```
where, all at `a=0`: `σ₀=σ_k(0)`, `σ₂=σ_k''(0)`, `σ₄=σ_k''''(0)`, `β₀=β_k(0)`, `β₁=β_k'(0)`, `β₃=β_k'''(0)`, and `φ₃:=V_k'''(0)`, `λ:=V_k''''(0)` (quartic). `(R)` is a **first-order nonlinear (Riccati) ODE in the scale** — the "running-mass flow" the task asks for.

**Numerical witness `[V-num]`** (`research/unconditional/w7/riccati_flow_witness.py`). (A) *Instantaneous identity check*: take the exact `(*)` RHS at a test potential `V₀` and compare `∂_kκ = ∂_a²( rhs )[0]` against the RHS of `(R)` computed from `V₀`'s `κ,φ₃,λ` and the model `σ,β`. Two configurations (with `σ₂,σ₄,β₀,β₃,φ₃` all nonzero to exercise every term):
```
  [1] d_k kappa = +0.600074  vs  (R) = +0.600000   rel.err 1.2e-4
  [2] d_k kappa = +0.210115  vs  (R) = +0.210000   rel.err 5.5e-4
```
i.e. `(R)` is the correct `a=0` projection of the exact E1 flow to discretisation order.

### 1.3 The Gaussian closure — the exactly closed anchor `[V-here]`
Set `σ=σ₀` const, `β=b₁ a` (so `σ₂=σ₄=0`, `β₀=β₃=0`), and `V_k` quadratic (`φ₃=λ=0`). Then `(R)` closes **exactly**:
```
   ∂_k κ = 2σ₀ κ² − 2b₁ κ ,      fixed points:  κ=0   and   κ* = b₁/σ₀ .        (G)
```
This is the unique self-contained Riccati submanifold; everything else couples to the unclosed moments. `(G)` is where the whole question localises.

---

## §2 · Gapless-trajectory test: is `κ≡0` a solution? is there a sign-definite source?

### 2.1 `κ≡0`: the free/Gaussian flow solution (`[T]`, acid-clean)
In `(R)`, `κ≡0` is a flow solution **iff `C(k)≡0`** along the trajectory. For the free/Gaussian marginal: `φ₃=λ=0` and, in the closed Gaussian case, `β₀=β₃=σ₄=0 ⇒ C=0`. Hence
> **`κ≡0` is exactly the free-field / massless trajectory of the flow.** Free field: `Γ_k^{(2)}=k²`, `V_{eff,k}''=0` for all `k`, `κ≡0`, `(*)`→`0`; the acid test **passes** (no manufactured gap). **[T]**, matches E1-FRG §2, MS-GAP-EQUATION §4, BANACH §5.

For the **interacting** marginal `λ≠0`, `C`'s leading term is `−σ₀λ ≠ 0`, so `κ≡0` is *not* a literal solution — **but then `(R)` does not close** (§1.3, §4); `λ` carries its own flow. So the honest statement is: **the gapless trajectory is a solution precisely in the Gaussian closure, and the interacting positive branch is a separate (unclosed) object.**

### 2.2 The crux: the `κ²` term is **marginal** — only the linear/source term decides (`[V-here]` + `[V-num]`)
Drop the linear/source terms: `∂_kκ = 2σ₀κ²`. Integrating, `κ(t)=1/(C−t)`, `t=\ln k`; as `k→0` (`t→−∞`), `κ→0` **but only logarithmically** (`κ∼1/(2σ₀|\ln k|)`). Witness: from `κ(1)=1`, `σ₀=0.5`, over `k∈[0.02,1]`, the exact solution gives `κ(0.02)=1/(2−0.02)=0.505` — *still O(½)*. **The `κ²` term can never hold a nonzero IR fixed point; it is marginal.** Therefore:
> **A positive IR fixed point is produced *only* by a definite sign of the linear coefficient `B=3σ₂−2β₁` and/or the source `C` at `κ=0`** — a genuinely **IR-relevant** direction. Witness `[V-num]` on `(G)`:
```
  b1 = 0.0  : only fixed point κ*=0 ; κ(0.02)=0.505   (marginal ⇒ GAPLESS)
  b1 =+0.4  : fixed point κ*=+0.8 ; κ(0.02)=0.880      (definite linear term ⇒ positive FP)
  b1 =−0.4  : fixed point κ*=−0.8 (unphysical); κ(0.02)=0.272 → 0   (⇒ GAPLESS)
```
This is the flow-coordinate image of MS-GAP-EQUATION §2.3: the positive branch exists **iff `F′(0)>1`** (a pitchfork / linear IR instability), equivalently **a sign-definite source at `κ=0`**.

### 2.3 Is the source at `κ=0` sign-definite? — **No, at the crossover** (`[V-doc]` + `[R]`)
- **β=0 / β<1/48: YES.** `β=0`: `V=Φ` convex (Haar), Prékopa `[T]` ⇒ `V_{eff,k}″≥N/12` uniformly in `k,L` (E1-FRG §3) — a **flow-invariant** positive fixed point `κ*≥N/12`. `β<1/48`: P3 residue theorem `[T]` (`PROOF-P3-UNCONDITIONAL`) gives `f₁²≥r₀>0`. The source is sign-definite positive here.
- **β≈2.5 at the coarse mode: NO.** The source sign is set by the coarse self-energy `Π(0)`: `κ>0 ⟺ Var(a;κ)≤1/λ ⟺ Π(0)≥λ` (**`H1`**, MS-GAP-EQUATION §2, TRANSPORT-CORE Lemma 2.2 lossless). The only microscopic generator of `C>0` is the **Savvidy tachyon** `λ_min=−g|B|`. But `[T]` (PROOF-SAVVIDY-EXACT-SIGN) `λ_{n,s}=k_∥²+(2n+1)g|B|+2s\,gB` is **transverse-Landau-gapped** at `√(g|B|)`, and at β=2.5 `[V-doc]` (W6-BE) the **coarse plane-wave block is positive** (`+0.40β`, 0 % negative) and **90–97 % of the negatives are high-momentum**, orthogonal to the `k→0` coarse mode. A fast-mode negative cannot gap the slow mode (PROOF-IR-H1 §4). Hence `F′(0)≤1`, `C` keeps the **gapless** sign, and `κ=0` stays the flow solution.
> **No sign-definite source at `κ=0` at the coarse mode at β≈2.5.** `[V-doc/T]`. This is the exact crux.

**Sign-definite invariant — does the flow preserve the sign?** The answer is **no**: since `C(k)` and `B(k)` are not sign-definite (they are functionals of `σ,β` = the coarse self-energy), the region `{κ≥0}` is **not invariant** under `(*)` at 2.5. A sign-preserving invariant `{κ≥κ₀>0}` would require `C(k)≥δ>0` (and `B` compatible) for all `k` on the coarse mode — that is precisely `H1`. **[V-here/O]**

---

## §3 · Uniform attainment across β≈2.5

The fixed point of `(G)` is `κ*(k)=b₁(k)/σ₀(k)` — **k-dependent** through `σ₀,β₁`; a *uniform* `κ*≥λ` for all `k→0` and `(L,a)` requires `Var(a;κ)≤1/λ` **uniformly**, i.e. a uniform coarse Poincaré constant on the `k→0` mode `[O]`.
- "Uniform across β≈2.5": β is a *lattice coupling*, not a flow parameter (category distinction, PROOF-FRG-DSE §2). The only rigorous finite-β object is the unclosed lattice measure; the only rigorous continuum object needs `a→0` uniform control = the **transport wall** (PROOF-IR-H1, TRANSPORT-MACHINE Stage 3).
- The Riccati value at the crossover is **not bounded uniformly**: the `σ,β` that set `C`, `B` are conditional supertraces over the *full* tangent space, whose Savvidy-negative part (`≈−1.1β`, W6-BE) is anti-localized to typical field × high momentum — the complement of the `k→0` gap cell (TRANSPORT-CORE §4, "anti-localization of positivity"). **Uniform attainment FAILS.** **[O]**

---

## §4 · Acid test + the three iron laws

**Free-field acid test (`[T]`, passes non-vacuously).** Free massless: `σ,β` of the Gaussian field, `λ=φ₃=0`, `C=0`, `κ≡0` is the unique flow solution, `(R)`→`0`. The tool **fails on the free field in the correct way** (returns `0`, not a false `N/12`): unlike a lattice `Γ₂` species whose positivity is the free-field Haar floor (TRANSPORT-CORE §3.2), here `σ_k` is the free *variance* and the IR fixed point is genuinely `κ=0`. **PASS.**

**Three iron laws（宪章 §2 / BRIEF §6）:**
1. **Label / no-QED.** Every line above is labelled `[T]/[V-here]/[V-doc]/[V-num]/[O]/[R]`. No step is called proven that is not. The route is **`[O]` at the crux**, **`[T]` for the `β=0`/`β<1/48`/free-field/generic-flow parts**. **No QED.**
2. **Free-field acid test.** Done above — the constructed flow/Riccati **fails on the free field** (`κ≡0`), so it is not vacuous; its `β=0` strengthening (`N/12`, flow-invariant) is a real theorem. **PASS.**
3. **No new name hiding the same wall.** The Riccati is the **flow-coordinate dialect of `H1`**: the positive fixed point `⟺ F′(0)>1 ⟺ Π(0)≥λ ⟺ Var(a;κ)≤1/λ ⟺ H1` (§2.3, MS-GAP-EQUATION §2, TRANSPORT-CORE Lemma 2.2 *lossless*). It is **not** a sixth species; it is Schema C + Schema A re-expressed. Stated as such, not disguised. **PASS (honest negative).**

**Backdrop (BRIEF §4, "三条硬结论").** The failure is the same three facts in flow dress: *anti-localization* (the positive `σ,β` compensation lives at high `k`/`‖F‖`, the sign at typical field × `k→0`); *two-coordinate wall* (positivity is Schur-concave/upper-bound, the atom needs the Schur-convex thimble side); *transport ≡ confinement* (uniform `κ*` across β≈2.5 is nonperturbative RG transport = confinement).

---

## §5 · The second face — the multiplicative residue / transport flow

The task's alternative "multiplicative transport resource flow" is `f₁²(μ)=Z_O(μ,μ₀)² f₁²(μ₀)` (registry A1; INVENT-TRANSPORT-INEQ). Its **flow form** is *linear*:
```
   ∂_k f₁²  =  2γ_O(k) f₁² ,    γ_O = −d\ln(β/2g)/d\ln μ = 2b₀g²+…     (M)
```
- **Sign-definite invariant: YES, trivially.** `Z_O²>0` (it is a *square*), so the flow **preserves the sign** of `f₁²`; the free field (`a≡0` by `Z₂`) is sent to `0`: acid-clean `[T]`.
- **But `0` is an absorbing fixed point.** If `f₁²(k₀)=0` (gapless at the crossover), `(M)` ⇒ `f₁²≡0` for all `k`. A *linear* flow (`A=0` Riccati — the degenerate/marginal case) `Z_O²` **multiplies a sign it cannot manufacture** (OTHER-F1SQ-TRANSPORT item 2). It is a **scale-covariance within one phase**: `f₁²(2.5)=f₁²_{ν_b}(β_eff)` at every block `b` (PROOF-TRANSPORT-RG §2), so it **presupposes** the single phase `(CP)` it would transport. **[T] identity, [O] the crossover bite.**
> So the residue flow is the **linear limit of `(R)`** (`A=0`): it reproduces exactly the marginal/gapless case of §2.2 — sign preserved, but no source, no positive fixed point without `(CP)=H1`.

---

## §6 · Verdict, exact failing step, scorecard

### 6.1 What is constructed (`[T]/[V-here]`)
- The **exact flow `(*)`** (E1-FRG) and its **`a=0` projection `(R)`** — a genuine Riccati-in-scale for `κ(k)`, **derived and numerically verified** (rel. err `5.5·10⁻⁴`). `[V-here]`
- The **`β=0` fixed point** `κ*≥N/12` as a **flow invariant** (`[T]`, E1-FRG §3), and the **β<1/48** residue fixed point (`[T]`, P3).
- The **α-residue flow `(M)`**: sign-preserving, free-field→0, linear feedback. `[T]`

### 6.2 The exact failing step (`[O]`)
> **No sign-definite source at `κ=0` at the coarse mode at β≈2.5.**
> Concretely: `κ≡0` is the (Gaussian) flow solution; the `κ²` term is **marginal** (§2.2); a positive IR fixed point needs `B` and `C` to take the **definite (IR-relevant) sign** at `κ=0` — a **sign-definite source**, equivalently `F′(0)>1` ⟺ a `k→0` Savvidy tachyon ⟺ `Π(0)≥λ` ⟺ `Var(a;κ)≤1/λ` ⟺ `H1`. Which sign the coarse mode carries is *itself* the open input. The tachyon is **transverse-Landau-gapped & 90–97 % high-momentum** at β=2.5 (`[V-doc]` W6-BE); the coarse block is positive, so the source keeps the **gapless** sign (`C(κ=0)` in the free/Gaussian direction). Additionally `(R)` **does not close** — it couples to `φ₃,λ,∂²σ,∂²β` (E1-FRG §4), the DSE vertex hierarchy (PROOF-FRG-DSE §2). **This is the precise missing step; it equals `H1`.**

### 6.3 Five-point scorecard (`2/5`, consistent with E1-FRG)
| # | condition | grade | basis |
|---|---|---|---|
| 1 | nonperturbative | **PASS** | exact `(*)` + exact `(R)`; no loop expansion; all β |
| 2 | sign lower bound `κ(k→0)≥λ` at 2.5 | **FAIL** | provable fixed point only gapless; source not sign-definite at coarse mode |
| 3 | non-circular | **PASS as identity, FAIL at closing** | positive FP ⟺ `F′(0)>1` ⟺ `H1` (lossless) |
| 4 | non-SUSY + fails free field | **PASS (clean)** | free field `κ≡0`; `β=0` `N/12` real but non-dynamical at 2.5 |
| 5 | uniform `(L,a)` `≥λ` | **FAIL** | `κ*(k)=b₁/σ₀` k-dependent; uniform ⟺ transport wall |

### 6.4 Relation to the corpus (no new species — iron law 3)
This route is **Schema C (`κ=F(κ)`, MS-GAP-EQUATION) seen as a scale flow**, and **Schema A (Banach–Riccati, BANACH-FIXED-POINT) seen through the `k`-variable**: the "neutral `k→0` mode" (BANACH §2) *is* the marginal `κ²` term; the "Savvidy makes `T` expansive" (BANACH §3) *is* the sign-indefinite source `C`; the "gapless contraction attractor" (BANACH §5, MS-GAP-EQUATION §4) *is* `κ≡0`. The residue flow is Stage 2 of TRANSPORT-MACHINE (`[T]`, sign-preserving, single-phase covariant). All three terminate at the **same** `H1` = `I2` = `H*` (PROFILE-FACET-ANCHOR §2.8; PROFILE §3). **Representational freedom is illusory; the flow coordinate is the same wall.**

### 6.5 Honest bottom line
- **Constructed OR exact failure?** Both, in the honest sense: the flow/Riccati **is constructed** (`(R)`, verified) and yields **real theorems** at `β=0`, `β<1/48`, and on the free field — but the **crux is an exact failure**: no sign-definite source at `κ=0` at the crossover, `(R)` non-closing, `κ≡0` the marginal gapless flow solution. **`[O]`, not `[T]`.**
- **Does it reduce to `H1`?** Yes — the fixed-point condition is `Π(0)≥λ`, lossless (TRANSPORT-CORE Lemma 2.2). **Not a sixth species.**
- **`∀class ≠ ∀possible`**: this is a no-go/`[O]` for the *flow/Riccati class* (Wetterich LPA-type + multiplicative transport), not a global impossibility theorem. A different, non-flow `k→0` sign-definite object is not excluded.

**Clay = 0 %. No QED. This is the flow-coordinate dialect of `H1`, with its single unfilled slot — a sign-definite source at `κ=0` — named exactly.**

---

## §7 · References
**External `[R]` (standard; not re-verified this session):**
- C. Wetterich, *Phys. Lett. B* **301** (1993) 90 (exact RG / `k∂_kΓ_k=½STr[(Γ^{(2)}+R)^{-1}k∂_kR]`).
- J. Polchinski, *Nucl. Phys. B* **231** (1984) 269 (exact flow, mass-term form).
- J. Berges, N. Tetradis, C. Wetterich, *Phys. Rept.* **363** (2002) 223 (hep-ph/0005122).
- N. Dupuis, L. Canet, A. Eichhorn, W. Metzner, J. M. Pawlowski, M. Reichert, D. F. Litim, *Phys. Rept.* **910** (2021) 1 (arXiv:2006.04853) — convexity `Γ_k^{(2)}+R_k⪰0`.
- D. F. Litim, *Phys. Rev. D* **64** (2001) 105007 (optimised regulator / threshold functions).
- G. K. Savvidy, *Phys. Lett. B* **71** (1977) 133; N. K. Nielsen, P. Olesen, *Nucl. Phys. B* **144** (1978) 376 (Savvidy spectrum `λ_{n,s}`).
- J. M. Cornwall, *Phys. Rev. D* **26** (1982) 1453; A. C. Aguilar, D. Binosi, J. Papavassiliou, *Phys. Rev. D* **78** (2008) 025010; Cyrol et al., *Phys. Rev. D* **94** (2016) 054005 (DSE/FRG gluon gap equation).
- W. T. Reid, *Riccati Differential Equations* (Academic, 1972) (Riccati comparison).

**Internal（absolute paths; labels carried）:**
- `D:\Dev\math\research\unconditional\E1-FRG.md` — the exact viscous-HJ flow `(*)` (`[T]`); convexity `≥0`; Prékopa `N/12` flow-invariant; §4 non-closure.
- `D:\Dev\math\research\unconditional\MS-GAP-EQUATION.md` — Schema C: `κ=F(κ)`, `κ=0` always, `F′(0)>1` pitchfork = `k→0` tachyon, absent at coarse mode at 2.5.
- `D:\Dev\math\research\unconditional\BANACH-FIXED-POINT.md` — Schema A: constrained-curvature Riccati `r′+r²=W−E`; neutral `k→0` mode; Savvidy expansive; gapless attractor.
- `D:\Dev\math\research\unconditional\PROOF-FRG-DSE-GAP-EQUATION.md` — DSE/FRG gluon gap eq; truncation/vertex; transfer fails.
- `D:\Dev\math\research\unconditional\W6-BE-at-2.5.md` — coarse plane-wave `+0.40β`, 0 % neg; negatives 90–97 % high-momentum; `G″(0)≈−0.165` proxy.
- `D:\Dev\math\research\unconditional\PROOF-SAVVIDY-EXACT-SIGN.md` — `λ_min=−g|B|`, transverse-Landau-gapped.
- `D:\Dev\math\research\unconditional\PROOF-IR-H1.md` §4 — fast mode cannot gap slow mode; `k→0` coarse Poincaré = H1.
- `D:\Dev\math\research\unconditional\TRANSPORT-CORE.md` — integrated `Γ₂` ⟺ `G″⪰ρ` ⟺ H1 (*lossless*, Lemma 2.2); anti-localization of positivity.
- `D:\Dev\math\research\unconditional\TRANSPORT-MACHINE.md` — machine factorization; Stage-2 multiplicative transport `[T]`; Stage-3 = H1/I2.
- `D:\Dev\math\research\unconditional\PROOF-TRANSPORT-RG.md` §2 — residue identical at every block (single-phase covariance).
- `D:\Dev\math\research\unconditional\OTHER-F1SQ-TRANSPORT.md` — `Z_O²` multiplies a sign it cannot manufacture; fails at I2.
- `D:\Dev\math\research\unconditional\PROOF-P3-UNCONDITIONAL.md` — `f₁²>0` at β<1/48 (I1).
- `D:\Dev\math\research\unconditional\PROFILE-FACET-ANCHOR.md` §2.3, §2.8 — Schema C/E gapless obstruction; residual DOF = `k→0` sign structure = H1.
- `D:\Dev\math\research\unconditional\PROFILE.md` §3 — the single open slot = `k→0` sign-definite self-consistency = H1.
- `D:\Dev\math\docs\new-tools-registry.md` §A1 — multiplicative residue transport `f₁²=Z_O²f₁²`, `γ_O=2b₀g²+…`.
- `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` §4, §6 — three hard conclusions; three iron laws.
- `D:\Dev\math\research\unconditional\w7\riccati_flow_witness.py` — numeric witness (`[V-num]`): `(R)` identity `5.5·10⁻⁴`; marginal `κ²`; fixed-point/source sign; run `py -3.12 -I`.

*Status: **FLOW/RICCATI slot attack.** Constructed the exact `a=0` Riccati `(R)` (verified) and the multiplicative residue flow `(M)`; both **fail at the same crux** — `κ≡0` is the marginal (Gaussian) flow solution, no sign-definite source at `κ=0` at the coarse mode at β≈2.5, `(R)` non-closing. `β=0`/`β<1/48`/free-field parts are `[T]`; the crossover is `[O]`. **The route reduces to Schema C + Schema A = the flow dialect of `H1`. Not a new species. Clay 0%; no QED.***
