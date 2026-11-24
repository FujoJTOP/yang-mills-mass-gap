# P3-CYLINDER-RESIDUE — the cylinder restatement of the P3 residue theorem, and the transfer-matrix/Källén–Lehmann residue bridge

> 2026-11-24 · Closes the **two audit gaps** of `PROOF-P3-UNCONDITIONAL.md` (see `P3-AUDIT-CORRECTED.md` items 4 and 5) to their honest depth.
> **Gap 1 (BLOCKING):** the file works on the four-torus `T⁴_L` (time periodic of length `L`) yet takes `f₁²_L = lim_{t→∞} e^{m₁ t} C_L(t)`; with periodic time the limit does not exist. → **§1: rigorous cylinder `T³_L×ℝ` restatement, exact residue definition, limit order.**
> **Gap 2:** the file asserts `f₁² = |⟨0|O|G⟩|² = a₀(1−O(u))` from `w_tube>0` "by modulus-square". → **§2: a transfer-matrix/Källén–Lehmann residue lemma; the standard parts proved, the polymer identification reduced to named inputs, and the isolation split off.**
> **Labels:** `[T]` theorem (proved here or in cited standard work) · `[C]` conditional theorem (precise input named) · `[O]` observation/factual · `[NG]` no-go · `[R]` required-but-missing step / open input · `[J]` judgment. Charter §2: **no `QED` on anything unproven; no fake closure.**

---

## VERDICT (up front)

**Gap 1 — the finite-volume time-periodic limit: CLEANED.** `[T]` The theorem is restated on the **cylinder** `T³_L × ℝ` (spatial torus, infinite time). The reflection-positive transfer matrix `T_L = e^{−aH_L}` on the gauge-invariant Hilbert space `H_phys` is, at finite spatial volume, **compact** (continuous kernel on a compact group) ⇒ its spectrum is **purely discrete**, `0 ≤ E_{0,L} < E_{1,L} < E_{2,L} ≤ … ↑ ∞`. The connected `tr F²` correlator has the exact Källén–Lehmann form `C_L(t)=Σ_{n≥1} f_{n,L}² e^{−E_{n,L}t}` with `f_{n,L}²=|⟨n,L|O|Ω_L⟩|²≥0`. The residue is **defined spectrally**, `f₁²_L := ν_L({E_{1,L}})=f_{1,L}²`, and — *because the finite-`L` spectrum is discrete* — equals `lim_{t→∞} e^{E_{1,L}t} C_L(t)`. **Order of limits:** `t→∞` at **fixed** spatial `L` (time first); the spatial size `L` enters **only as the uniformity axis** `inf_L f₁²_L`, never as an `L→∞` limit needed for the residue. On the four-torus this is impossible; on the cylinder it is exact.

**Gap 2 — the spectral-residue bridge: PROVED for the standard parts, REDUCED to two named inputs.** `[T]/[C]/[R]` The bridge splits cleanly:
- **Transfer-matrix / KL content** — spectral representation, positivity `f_{n,L}²≥0`, finite-`L` atom extraction, "large-`t` rate = lightest atom" — is **standard** `[T]` (spectral theorem for a compact self-adjoint `T`, made positive by Osterwalder–Schrader reflection positivity; Osterwalder–Seiler '78, Seiler LNP 159).
- **Leading-tube → leading-atom identification** — that the P3 tube weight `a₀ u^{4t/a}` is the *same object* as the leading KL coefficient — is **not** "modulus-square" alone; it needs two inputs: **(R1)** the *dressed-mass / amplitude resummation* (the "genuine remaining gap" of `P3-AUDIT-CORRECTED` §1.3) `[R]`, and **(R2)** the *full-representation KP estimate* `[T]` (standard, Seiler, closable). Given (R1)+(R2), the identification holds **at each finite `L`**, unconditionally in structure.
- **Isolation `I2`** — that the leading atom is an **isolated spectral point** (gap `Λ−m₁²≥δ>0`), *not a continuum threshold*, and that the atom survives the `L→∞` limit with a **uniform** gap — is **exactly the open input `I2`** `[R]`. **It is NOT needed for the finite-`L` residue floor** (where discreteness already isolates every atom); it **IS needed** for (i) the statement "atom, not a continuum threshold" *in infinite volume*, and (ii) the infinite-volume extraction `f₁²_∞=lim_{t→∞}e^{m₁t}C_∞(t)`.

**Net.** `[J]` The residue theorem can be honestly stated on the cylinder as a **finite-`L`, uniform-in-`L`** theorem: for each `T³_L`, `f₁²_L ≥ r₀(β)>0`. This is **I1 only**, and it needs **not** `I2`. The "isolated atom / not a threshold / survives to `L=∞`" upgrades are `I2`, still open. **No falsetto closure.**

| Object | Standard (transfer-matrix/KL) | Needs P3 polymer | Needs dressed-mass (R1) | Needs `I2` |
|---|---|---|---|---|
| `C(t)=Σ f_n²e^{−E_nt}`, `f_n²≥0` | **[T]** | — | — | — |
| finite-`L` atom extraction `lim_t e^{E₁t}C` | **[T]** (discreteness) | — | — | — |
| rate `=E_{1,L}=inf supp ν_L` | **[T]** | — | — | — |
| `f₁²_L=a₀(1−O(u))` (value, finite `L`) | **[T]** structure | **[T]** geometry | **[R]** | — |
| atom **not threshold**, `L=∞` | — | — | — | **[R]** |
| `f₁²_∞=lim_t e^{m₁t}C_∞`, uniform `L→∞` | — | — | — | **[R]** |

---

## §1 The cylinder: setup, transfer matrix, Källén–Lehmann measure, exact residue, limit order

### 1.1 Why `T⁴_L` fails — the periodic-time obstruction `[T]/[O]`

On `T⁴_L` every direction, including time, has period `L`. The transfer matrix in the time direction gives the **periodic** (Gibbs-on-a-circle) two-point function: with `T|n⟩=e^{−E_n a}|n⟩`, `Ω` the vacuum (`E_0=0`),
```
C_L^{per}(t) = Σ_{n≥1} f_{n,L}² · [ e^{−E_{n,L} t} + e^{−E_{n,L}(L−t)} ] / ( 1 − e^{−E_{n,L} L} ),   0 ≤ t ≤ L.
```
This is **even about `t=L/2`**: the direct propagator `e^{−E_n t}` is balanced by its image `e^{−E_n(L−t)}` (the heat kernel on a circle of circumference `L`). It is bounded and has **no decay direction**; `t` ranges over `[0,L]` and never reaches `∞`. Hence `lim_{t→∞} e^{m₁ t} C_L(t)` is **not defined** on `T⁴_L` — the expression only makes sense as `L→∞` inside the time circle, which is a different (joint) limit that fixes nothing. This is exactly audit item 4. `[T]`

### 1.2 The cylinder and the transfer matrix `[T]`

**Geometry.** Spatial lattice `Λ_s = (aℤ/Lℤ)³` (three-torus, side `L`, `N_s=(L/a)³` sites, spatial links `E_s`, `|E_s|=3N_s`), compact configuration space `X_L = G^{E_s}` (`G=SU(N)`), Haar measure. Time is **infinite** (`ℝ`, discrete spacing `a`): the cylinder `𝒞_L := T³_L × ℝ`. Wilson action `S=β Σ_p (1 − N⁻¹ Re tr U_p)`, `β=2N/g₀²`.

**Transfer matrix.** By Osterwalder–Schrader reflection positivity of `(S,Haar)`, there is a **self-adjoint contraction** `T_L = e^{−aH_L}` on the **physical (gauge-invariant) Hilbert space**
```
H_phys,L := L²(X_L)^{G}   (wavefunctions invariant under time-independent spatial gauge transformations; = ker of the Gauss-law constraint),
```
with `0 ≤ T_L ≤ 1`, `H_L ≥ 0`, and (Perron–Frobenius, positive kernel) a **unique** top eigenvector `Ω_L` (the vacuum, `T_L Ω_L = Ω_L`, `E_0=0`). [Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Seiler, LNP **159** (1982), Ch. IV; Osterwalder–Schrader, *CMP* **31** (1973) 83.] `[T]`

**Compactness ⇒ discrete spectrum.** The kernel `T_L(U,U')` is a continuous function on the compact space `X_L×X_L` (finite sum/integral of characters of `U,U'`), hence **Hilbert–Schmidt**, hence **compact**; the restriction to the closed invariant subspace `H_phys,L` is compact. Therefore
```
spec(T_L) = { 1 = λ_{0,L} > λ_{1,L} ≥ λ_{2,L} ≥ … } ∪ {0},   λ_{n,L} ↓ 0,
E_{n,L} := −(1/a) log λ_{n,L} ∈ (0,∞),   0 = E_{0,L} < E_{1,L} ≤ E_{2,L} ≤ …,   E_{n,L} ↑ ∞.
```
This is the **decisive** structural fact: **at finite spatial volume the spectrum is purely discrete** — there is *no continuum*, every spectral point is an **atom**. (Compactness of `T_L` holds **only** for `L<∞`; as `L→∞` the spectrum may acquire a continuous part — see §1.4, §3.) `[T]`

### 1.3 Källén–Lehmann measure and the exact residue `[T]`

Let `O := (2N)^{−1} Σ_{p∈A1g}(χ_f+χ_f̄)(U_p)` be the `A1g`/zero-momentum gauge-invariant single-plaquette density on one time slice (`O∈H_phys,L`, `O=O†`). The **connected** cylinder correlator (infinite time; vacuum expectation, vacuum piece subtracted) is
```
C_L(t) := ⟨Ω_L| O T_L^{t/a} O |Ω_L⟩ − ⟨Ω_L|O|Ω_L⟩² ,   t ∈ aℤ_{≥0}.
```
Insert the spectral resolution of the compact `T_L` and use `f_{n,L}²:=|⟨n,L|O|Ω_L⟩|²`:
```
C_L(t) = ∫_0^∞ e^{−s t} dν_L(s) = Σ_{n≥1} f_{n,L}² e^{−E_{n,L} t},   ν_L := Σ_{n≥1} f_{n,L}² δ_{E_{n,L}} ,  ν_L ≥ 0.
```
This is the **Källén–Lehmann (spectral) representation** on the lattice. `f_{n,L}²≥0` is the spectral theorem (`f_{n,L}²` is a modulus-square of a matrix element); `ν_L≥0` is Osterwalder–Schrader reflection positivity. [Källén, *Helv. Phys. Acta* **25** (1952) 417; Lehmann, *Nuovo Cim.* **11** (1954) 342; OS'78; Seiler LNP 159.] `[T]`

**Exact residue definition.**
```
m_{1,L} := E_{1,L} = inf supp ν_L          (lightest O-coupled energy above the vacuum),
f₁²_L   := ν_L({m_{1,L}}) = f_{1,L}² = Σ_{E_{n,L}=E_{1,L}} |⟨n,L|O|Ω_L⟩|²  .
```
Two harmless precisions. (i) `m_{1,L}` is the lightest energy with **nonzero `O`-overlap** (the lightest state *coupled to `tr F²`*); physically this is the lightest `0⁺⁺` glueball, and the P3 tube shows it is O-coupled (§2), so `m_{1,L}` is not a "dark" state. (ii) If `E_{1,L}` is degenerate, `f₁²_L` is the **total** atom weight at that level — a sum of non-negative terms. This is the robust quantity: `ν_L({m_{1,L}}) ≥ |⟨G_L|O|Ω_L⟩|²` for *any single* eigenstate `G_L` at the level, so a lower bound on `ν_L({m_{1,L}})` is the natural statement; the single-state Clay form `|⟨0|O|G⟩|²` coincides with it when the lightest level is non-degenerate (expected — PF pins only the vacuum, not the next level, but the ground glueball is generically simple). `[T]/[O]`

**Two equivalent descriptions of `f₁²_L` (finite `L`).** Since `E_{2,L}>E_{1,L}` at finite `L` (discrete, distinct),
```
f₁²_L = lim_{t→∞} e^{E_{1,L} t} C_L(t)         (the extraction formula — valid at each finite L),          (*)
f₁²_L = ν_L({E_{1,L}})                          (the spectral definition — needs no limit).
```
Both are legitimate **because the finite-`L` spectrum is discrete**. The extraction `(*)` is the general Källén–Lehmann residue of the atom at `E_{1,L}`; the spectral form is what the compactness of `T_L` guarantees directly. `[T]`

### 1.4 Order of limits — exactly which one the residue needs `[T]`

There are two independent parameters: time `t` (cylinder direction) and spatial size `L` (uniformity parameter). **The residue needs time first.**
```
[fixed L]  ──t→∞ (or: exact finite-L spectral projection)──▶  f₁²_L  ──inf over L──▶  uniform bound inf_L f₁²_L ≥ r₀ > 0.
```
- **Time first, at fixed `L`.** The extraction `(*)` requires `t→∞` with `L` held **finite**. On the cylinder this is available; on `T⁴_L` it is not (§1.1). If one insists on the finite slab `T³_L×[−T,T]`, the `T→∞` limit is the ground-state projection `T_L^{T/a}→|Ω_L⟩⟨Ω_L|` and **must precede** any `L→∞` — it is exactly the spectral projection used above, so no boundary-layer subtlety survives.
- **`L` is a uniformity axis, not a limit.** Theorem 2.5's content is `inf_L f₁²_L ≥ r₀(β)` over **finite** `L`; it never takes `L→∞`. The **infinite-volume residue** `f₁²_∞ = lim_{t→∞} e^{m₁t} C_∞(t)` is a *separate, stronger* statement that requires (a) `f₁²_L` to converge/`liminf`, and (b) the atom at `m₁` to remain **isolated** in the infinite-volume spectrum — i.e. **`I2`** (§3). It is **not** part of the residue theorem and is **not** claimed. `[T]/[R]`

Thus the correct order for the *residue theorem* is **`L` fixed → `t→∞` → `inf_L`**, i.e. time-first with `L` only as uniformity; the reverse (`L→∞` before time extraction) is what the four-torus forced and is exactly the ill-posed object.

---

## §2 The transfer-matrix / Källén–Lehmann residue lemma

### 2.1 Statement

**Lemma (KL residue bridge).** Work on the cylinder `T³_L×ℝ` of §1, `T_L=e^{−aH_L}` compact on `H_phys,L`, `O` gauge-invariant. Then:

(a) **[T] Spectral representation.** `C_L(t)=Σ_{n≥1}f_{n,L}²e^{−E_{n,L}t}`, `f_{n,L}²≥0`, `ν_L≥0`, with `0<E_{1,L}<E_{2,L}≤…↑∞` discrete (§1.2–1.3).

(b) **[T] Finite-`L` atom extraction.** Because `E_{2,L}>E_{1,L}`, `f₁²_L:=ν_L({E_{1,L}})=lim_{t→∞}e^{E_{1,L}t}C_L(t)`. In particular the large-`t` exponential of the connected correlator is the lightest **atom** `E_{1,L}`, and its coefficient is the atom weight `f₁²_L`.

(c) **[T] Rate = support infimum.** `−lim_{t→∞}(1/t)log C_L(t)=E_{1,L}=inf supp ν_L` (given `f₁²_L>0`).

(d) **[C] Polymer identification (P3).** In the KP regime `β<β₁` the connected-`tr F²` expansion computes, **exactly** and **uniformly in `L` and `t`**,
```
C_L(t) = a₀(β) u(β)^{4t/a} (1 + O(u)) + O(e^{−δ' t}),   a₀(β) > 0  (L-independent),
```
the leading term being the **four-walled tube** (`m₁^{(0)}=4 ln(1/u)/a`, factor-4 corrected geometry; Gupta §14.2.2, Münster). Matching the two *exact* representations of the same `C_L(t)` — (a–c) and the polymer sum — gives, **at each finite `L`**,
```
E_{1,L} = 4 ln(1/u)/a + O(u/a)   (dressed),      f₁²_L = a₀(β)(1 − O(u)) ≥ ½ a₀(β) =: r₀(β) > 0,
```
uniformly in `L`. `[C]` — conditional on **(R1)** the dressed-mass / amplitude resummation (`P3-AUDIT-CORRECTED` §1.3) and **(R2)** the full-rep KP estimate (`[T]`, Seiler, closable).

### 2.2 Proof of the standard parts (a)–(c) `[T]`

*Self-adjointness and positivity.* OS reflection positivity of `(S,Haar)` gives `T_L` self-adjoint and `0≤T_L≤1`, and makes the OS Hermitian form positive semidefinite on `H_phys,L`; hence the inner product under which `T_L` is self-adjoint is genuine, `H_L≥0`, and the spectral theorem applies. This is the lattice Osterwalder–Seiler construction. `[T]`

*Compactness.* §1.2: continuous kernel on compact `X_L` ⇒ Hilbert–Schmidt ⇒ compact; restriction to `H_phys,L` compact. Compact self-adjoint ⇒ orthonormal eigenbasis `{|n,L⟩}` of `H_phys,L`, real eigenvalues `λ_{n,L}↓0`, `E_{n,L}=−(1/a)log λ_{n,L}∈(0,∞)`, `E_{n,L}↑∞`. `[T]`

*(a).* Expand `C_L(t)=⟨Ω_L|O T_L^{t/a}O|Ω_L⟩−⟨Ω_L|O|Ω_L⟩²` in `{|n,L⟩}`; `T_L^{t/a}|n,L⟩=e^{−E_{n,L}t}|n,L⟩`; the subtracted vacuum term removes `n=0`. Coefficients `f_{n,L}²=|⟨n,L|O|Ω_L⟩|²≥0`. `[T]`

*(b).* `e^{E_{1,L}t}C_L(t)=f₁²_L + Σ_{n≥2}f_{n,L}²e^{−(E_{n,L}−E_{1,L})t}`. Since `E_{n,L}−E_{1,L}≥E_{2,L}−E_{1,L}=:δ_L>0` (finiteness ⇒ `δ_L>0`, though possibly `L`-small), the tail `→0` as `t→∞` **at fixed `L`**. Hence the limit is `f₁²_L`. (Note: `δ_L>0` at each finite `L` is all (b) uses; *uniformity* `δ_L≥δ>0` is **not** used and is **not** available — that is `I2`.) `[T]`

*(c).* From (b), `(1/t)log C_L(t) = −E_{1,L} + (1/t)log[f₁²_L + o(1)] → −E_{1,L}` when `f₁²_L>0`; and `ν_L` is supported on `{E_{n,L}}`, so `inf supp ν_L=E_{1,L}`. `[T]`

### 2.3 The polymer identification (d) — why "modulus-square" was insufficient, and what supplies it `[C]/[T]`

The audit objection (item 5): "*from `w_tube>0`, `a₀>0` one cannot pass to `f₁²=a₀(1−O(u))` by 'modulus-square' alone; one needs a transfer-matrix/KL proof that the leading exponent is the lowest relevant spectral state, not a continuum threshold; that the amplitude is the square of a spectral projection matrix element; and that it is volume-uniform.*" The bridge answers each clause, and identifies precisely where standard theory ends and inputs begin:

1. **"Amplitude = modulus-square of a spectral projection"** — *this is (a) itself*: `f_{n,L}²=|⟨n,L|O|Ω_L⟩|²` is the spectral theorem. So the modulus-square is **not** an extra assumption; it is *what `f₁²` is*. `[T]` The defect in the original file was to say "`f₁²` is a modulus-square, so positive" **without** connecting the polymer number `w_tube=a₀u^{4t/a}` to this modulus-square. That connection is the content of (d): a₀ is the *coefficient* of the leading exponential of the same `C_L(t)`, hence equals `ν_L({E_{1,L}})=f₁²_L`. `[J]`

2. **"Leading exponent = lowest relevant spectral state, not a continuum threshold"** — At **finite `L`** this is **(b)/(c)** and is **automatic**: compactness makes the spectrum discrete, so *every* spectral point is an atom; the leading exponential of the exact `C_L(t)` is necessarily the lightest atom `E_{1,L}`, never a threshold density. `[T]` The phrase "not a continuum threshold" is therefore **vacuously true at finite `L`** and becomes a *real* statement only for `L=∞` — which is exactly `I2` (§3). This is the sharp place where the audit's "非连续阈" demand lives: **satisfied for the residue theorem (finite `L`), open for the infinite-volume atom (`I2`).**

3. **"The leading polymer tube *is* that atom"** — the polymer expansion is, in the KP regime, an **exact, absolutely convergent** representation of `C_L(t)` (not merely asymptotic): `C(t)=Σ_{p,q}Σ^conn_{P∋p,q}w(P)` with `Σ_{P∋x}|w(P)|e^{a|P|}≤a` (`u·K<1`, `K≤CN`). Since it computes the *same* function as the exact spectral sum (a), its leading exponential must be `E_{1,L}` and its coefficient must be `f₁²_L`. The minimal-cost object is the **four-wall tube** (factor-4 geometry: four plaquettes per time slice, `m₁^{(0)}=4ln(1/u)/a`; the "one-per-slice" object has weight **zero** — uncontract'd source links — not merely a loose bound [Gupta §14.2.2; Münster]). Hence `E_{1,L}=4ln(1/u)/a+O(u/a)` and `f₁²_L=a₀(1−O(u))`. `[T]` for the *geometry/uniqueness* (Lemma 2.1 of `P3-AUDIT-CORRECTED`); `[C]` because the **dressed-mass / amplitude resummation (R1)** and the **full-rep KP (R2)** must be in place for the expansion to be uniformly convergent and for the `O(u)` remainder to be `t`-uniform (this is `P3-AUDIT-CORRECTED` §1.3 `[R]`).

4. **"Volume-uniform"** — `a₀(β)` is built from `L`-independent positive fundamentals `c_f>0` and the `f⊗f̄→1` singlet recoupling `>0`; the KP criterion is **independent of `L`**; hence `f₁²_L≥r₀(β)` with `r₀` `L`-independent. `[C]`

**Sharp statement of what the bridge proves vs. imports.**
- Proved here / standard `[T]`: the residue **is** the leading atom weight (a); at finite `L` it **is** an atom, not a threshold, and its **value** is the leading coefficient (b,c); the leading polymer tube is the unique leading concentration with positive amplitude (the geometry `[T]`, the positivity `[T]`).
- Imported `[R]`: the dressed-mass/amplitude resummation (R1) — *existence of a `t`-uniform leading exponential with a positive amplitude floor `1−O(u)`* — the one non-printed step of `P3-AUDIT-CORRECTED` §1.3. Full-rep KP (R2) is `[T]` (factorial domination of `c_r`, Seiler).
- **Not** required for the residue identification: `I2` (§3). The cost grading is a statement about polymer **amplitudes**, not transfer-matrix **eigenvalue spacing** (`PROOF-P1P2` §3; `PROOF-EDGE-E2`); it cannot and need not produce a uniform gap `δ` for `f₁²_L≥r₀`.

### 2.4 Consistency check with the audit's own warning `[J]`

`PROOF-P1P2` §4 already flags: the `O(u)` remainder bound `e^{m₁t}R_L=O(u)` is *"a `t`-uniform multiplicative bound, satisfied by any `λ₂/λ₁∈[1−Cu,1)` — hence it certifies the residue I1, not the isolation I2."* This is **consistent** with the lemma above and sharpens it: the `O(u)` bound fixes the **coefficient** `f₁²_L` (the value, clause (d)) but is **compatible with `δ_L=E_{2,L}−E_{1,L}→0`** as `L→∞`. That is precisely why (b) uses only `δ_L>0` *at finite `L`* (not a uniform `δ`): the coefficient is a value, the gap is a separate object. **The residue identification does not secretly smuggle in a gap.** `[T]`

---

## §3 What remains open — `I2`, and the honest boundary

### 3.1 `I2` (isolation) is the exact open input for the infinite-volume statements `[R]`

`I2` (from `PROOF-ISOLATION`, `GAP-G2-edge`, `PROOF-EDGE-E2`, `mass-spectrum` Theorem A) is
```
I2:   Λ_L − m₁²_L ≥ δ > 0  uniformly in L,   where Λ_L := inf supp(ν_L − f₁²_L δ_{m₁,L}),
```
equivalently a uniform two-exponential / two-mass asymptotic `C_L(t)=A_L e^{−m₁t}(1+η_L(t))`, `|η_L(t)|≤B e^{−δt}`, `A_L≥a₀`. `[R]`

**What `I2` is needed for** (and **not** needed for):
- **Needed** for the claim that `E_{1,L}` is an **isolated atom in infinite volume**, i.e. that the glueball persists as a stable particle with a gap to the next spectral point (the excited `0⁺⁺*` at `≈1.5m₁` or the multi-glueball threshold `≥2m₁`), rather than dissolving into a continuum that touches `m₁`. `[R]`
- **Needed** for the infinite-volume extraction `f₁²_∞=lim_{t→∞}e^{m₁t}C_∞(t)` to exist and equal the atom weight: if the infinite-volume continuum touches `m₁` (a threshold), the atom weight is `0` and the limit fails or diverges. `[R]`
- **NOT needed** for the **finite-`L` residue floor** `f₁²_L≥r₀(β)`: at each finite `L`, compactness makes `δ_L>0` automatically, and the *value* `f₁²_L` is bounded by the polymer amplitude (clause (d)), independent of how small `δ_L` is. **The residue theorem (I1) stands without `I2`.** `[T]/[C]`

This matches `PROOF-P3-UNCONDITIONAL §5`'s own correction: "the cost grading gives I1 only, not I2." The cylinder/spectral bridge makes the reason transparent: *I1 is a statement about the atom's **weight** (a coefficient), I2 about the **spacing** to the next atom; the polymer expansion controls the former, not the latter.*

### 3.2 The remaining named inputs, precisely `[R]`

To promote Theorem 2.5 from `[C]` to `[T]` one must supply, in addition to §1–§2:
1. **`I2` (isolation)** — the uniform-in-`L` second-level gap / two-mass bound. Open (`PROOF-ISOLATION`, `PROOF-EDGE-E2`). **Not** needed for I1.
2. **Dressed-mass / amplitude resummation (R1)** — existence of a `t`-uniform effective exponential and a positive amplitude floor `1−O(u)` (Glimm–Jaffe–Spencer two-point analysis adapted to strong-coupling gauge theory). Standard-in-principle, unwritten (genuine gap of `P3-AUDIT-CORRECTED` §1.3).
3. **Full-rep KP (R2)** — factorial domination `|c_r|≤K'β^{|r|}/|r|!` × poly(`dim r`) ⇒ finite `C_g(β,N)`. `[T]`, Seiler LNP 159, closable.

**Not touched by this note (out of scope, unchanged):** continuum `a→0`, transport to `β≈2.5` (the analyticity wall), and any conversion of the strong-coupling residue into a gap theorem. Clay remains open.

### 3.3 Bottom line `[J]`

- **Gap 1 (cylinder): closed `[T]`.** The theorem lives on `T³_L×ℝ`; `f₁²_L` is defined spectrally as `ν_L({E_{1,L}})` and — because the finite-`L` transfer matrix is compact — equals `lim_{t→∞}e^{E_{1,L}t}C_L(t)`. **Limit order: `t→∞` at fixed `L`, then `inf_L`; `L→∞` is never taken for the residue.** The four-torus expression was ill-posed; the cylinder one is exact.
- **Gap 2 (bridge): reduced to named inputs.** The transfer-matrix/KL parts (spectral representation, positivity, finite-`L` atom extraction, rate = atom) are standard `[T]`. The leading-tube ↔ leading-atom **identification** is exact at finite `L` *given* (R1) dressed-mass resummation `[R]` and (R2) full-rep KP `[T]`; "modulus-square" alone was insufficient because it did not connect `w_tube` to the spectral coefficient — the bridge (matching the two exact representations of the same `C_L(t)`) supplies the connection. The **isolation `I2`** ("atom, not threshold", "survives to `L=∞`") is **separate and open**, and is **not** required for the finite-`L` residue floor.

*Refs.* Källén, *Helv. Phys. Acta* **25** (1952) 417; Lehmann, *Nuovo Cim.* **11** (1954) 342; Osterwalder–Schrader, *CMP* **31** (1973) 83, *CMP* **42** (1975) 281; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440–471 (transfer matrix, OS positivity, uniform strong-coupling gap); Seiler, *Gauge Theories as a Problem of Constructive QFT and Statistical Mechanics*, LNP **159** (1982) Ch. IV (transfer matrix, compactness, character expansion); Kotecký–Preiss, *CMP* **103** (1986) 491–498; Glimm–Jaffe–Spencer, *Ann. Math.* **100** (1974) 585 (two-point / particle structure); R. Gupta, *Introduction to Lattice QCD*, arXiv:hep-lat/9807028 §14.2.2, Eqs. (14.20)–(14.21) (four-wall minimal tube); Münster, *Phys. Lett. B* **95** (1980) 59 / *Nucl. Phys. B* **190** (1981) 439, erratum B200 (1982) 536 (`m=−4 ln u+…`); Shen–Zhu–Zhu, *CMP* **400** (2023) 805 / arXiv:2204.12737 (strong-coupling gap, gap only). Internal: `PROOF-P3-UNCONDITIONAL`, `P3-AUDIT-CORRECTED` (items 4, 5; §1.3 (R1); Lemma 2.1 factor-4), `PROOF-P1P2` (§3, §4), `PROOF-P3`, `PROOF-ISOLATION` (`I2`), `GAP-G2-edge` (residue–level-gap lemma), `PROOF-EDGE-E2` (bottom-functional blindness), `PROOF-F1SQ` (gapped free field `f₁²=0`), `BS-RP-STRUCTURE`/`BS-THRESHOLD` (KL measure, `H_phys`, thresholds), `mass-spectrum`/STATUS (Theorem A, `I1∧I2⟺` sharp two-exponential), `ATOM-strong-coupling-extend` (analyticity wall).
