# W2-MASTER-FIELD — the master-field equation for `ψ₁` is **closed** (the linearised Makeenko–Migdal / Bethe–Salpeter eigen-equation around `W₀`), `w₁∞=g₁/Σg` is **the point-mass of the normalised Källén–Lehmann measure at its lightest glueball** `=1/(1+R)`, it is **NOT bounded below uniformly in the 't Hooft coupling** (the free/UV end is `0` — exactly the free-field acid test), and the correct uniformity is **fixed physical `β=2.5`, in `(L,a)`** = the `W` target `[O]`; large `N` supplies the `N⁰` scaffold and the `1/N²`-leakage structure, **not the value**

> **2026-11-24 · Wave W2 — the large-`N` master field for the first excited eigenvector / source overlap.**
> Charter `W-WAVE-CHARTER.md` §1 W2 asks: derive the master-field equation for `ψ₁`, compute
> `w₁∞=g₁/Σg`, decide whether `g₁/Σg` has a positive lower bound **uniform in the 't Hooft coupling**,
> and whether the `1/N²` correction preserves positivity. Starting point is U3 (`w₁=O(N⁰)`,
> `w₁∞=g₁/Σg`). This note does exactly that: writes the closed eigen-equation for `ψ₁` (linearised
> Makeenko–Migdal = planar Bethe–Salpeter), reduces `g₁/Σg` to a single ratio `1/(1+R)` (leakage
> ratio `R`) = the atom weight of the normalised KL measure, evaluates it exactly in the solvable
> limits, shows the λ-uniform floor is **false** (and *must* be, by the free-field acid test), and
> pins the correct uniformity. **Reads.** `W-WAVE-CHARTER.md`; `U3-LARGE-N-BINDING.md`,
> `U1-HARNACK-BOUND.md`, `U-WAVE-SYNTHESIS.md`; `T2-TUBE-OVERLAP.md`, `T-WAVE-SYNTHESIS.md`;
> `R2-ISOLATION.md`, `T3-TP2-ISOLATION.md`; `unconditional/PNT-LARGEN-LOOP.md`,
> `CONSTRUCT-MASTERFIELD.md`, `KAPPA-LARGEN.md`, `ROOF-LARGEN-MASTER.md`, `G1-LARGE-N.md`,
> `PROOF-LARGE-N-RIGOROUS.md`. **Witness.** `research/v3/w7/w2_master_field.py`
> (`py -3.12 -I`, numpy only, deterministic, **18/18 PASS**; output quoted §8). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here · **[std]** textbook tool used
> as such. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **The master-field equation for `ψ₁` is closed and explicit `[T]` (structure) + `[O]` (content).**
   At leading `N` the connected `0⁺⁺` two-point function of `O=(1/N)tr F²` obeys a **closed** equation
   — no `n`-point tower — because large-`N` factorization turns the loop-equation RHS into a product
   (`PNT-LARGEN-LOOP` `[T]`). The glueball masses are the spectrum of the **linearisation `ℳ` of (MM∞)
   around the master field `W₀`**, and `ψ₁` is its lightest `0⁺⁺` eigenvector:
   ```
       ℳ ψ₁ = m₁² ψ₁ ,      ℳ = δ(MM∞)/δW |_{W₀}   (restricted to the 0⁺⁺, k=0 sector) ,          (MF)
   ```
   equivalently the **planar Bethe–Salpeter equation** `[G₀⁻¹(−m₁²) − K(−m₁²)]ψ₁=0` with `K` the
   `0⁺⁺` two-particle-irreducible kernel built from `W₀`'s vertices. **The equation closes** (`[T]`);
   **solving it needs `W₀`, unconstructed in 4D** (`[O]`; Singer constructed `W₀` in 2D only,
   `PNT-LARGEN-LOOP §3`). The source overlap is `g₁=⟨ψ₁|O|Ω⟩²`.

2. **`g₁/Σg` is computed, and it is one number: the atom weight of the normalised KL measure `[T]`.**
   With `dν(s)=Σ_n g_n δ(s−m_n²)+dν_cut` the Källén–Lehmann measure of the connected `O`-correlator,
   `M₀=∫dν=Σ_n g_n` (at fixed renormalisation scale, `U4`), `f₁²=g₁`, one has the exact identity
   ```
       w₁∞ = g₁/Σ_n g_n = 1/(1+R) ,     R := Σ_{n≥2} g_n / g₁  (leakage ratio) ,                     (RATIO)
       = ν({m₁²}) / ν([0,∞)) = ρ({m₁²}) ,   ρ := dν/M₀  the normalised KL measure .
   ```
   So **`g₁/Σg` is literally the mass the KL measure puts on its lightest point.** A positive floor
   `g₁/Σg ≥ w₀` is **exactly** the leakage bound `R ≤ (1−w₀)/w₀` — i.e. `Σ_{n≥2}w_n ≤ 1−w₀` of `T2`,
   the *same* object. Exact values: strong coupling `s/2` (`=1/2` for the physical `s=1`, `[T]` mod `[KL]`),
   2D `tr F²` `0`, free field `0`, holonomy carrier `1` (§3, §7).

3. **`g₁/Σg` is NOT bounded below uniformly in the 't Hooft coupling `[T]/[NG]`.**
   The literal statement "`∃w₀>0` uniformly in `λ=g²N`" is **false**, and it **must be**: the free
   massless field is the `λ→0` (UV/asymptotic-freedom) end of the same family, and it has `g₁/Σg=0`
   *exactly* (`w₁=0`, the acid test). Any λ-uniform floor would be **vacuous on the free field** — the
   acid test forbids it. Concretely (§4, witness C) `g₁/Σg(λ=0)=0` **exactly** (free field, `[T]`) and
   rises to `s/2` at strong coupling, so `inf_{λ≥0} g₁/Σg = 0`. **This is not a defeat — it is the acid
   test working.**
   The **correct** uniformity — the one `W` actually needs — is at **fixed physical coupling
   (`β=2.5`), uniform in `(L,a)`**, and that is the `W` target `[O]`. Three distinct readings of
   "uniform in `λ`" are separated and labelled in §4.

4. **The `1/N²` correction preserves *positivity* trivially, but not the *floor* `[T]/[O]`.**
   `w₁=f₁²/M₀` is a **fraction** (`f₁²≥0`, `M₀≥f₁²`), so `w₁(N)>0` and `w₁(N)≤1` at **every** finite
   `N` by construction (`[T]`). But the **uniform floor** `w₁(N)≥w₀` is **not inherited** from the
   large-`N` value: (i) the asymptotic `1/N²` series has **no uniform remainder at `N=3`**
   (`PROOF-LARGE-N-RIGOROUS` (c): `|c|/N²=O(0.1)` at `N=3`, witness E); (ii) the **sign** of the
   `1/N²` coefficient `c` is **not fixed by factorization** — the multi-trace leakage / connected
   4-point `Ĝ` (`U3 §5.2` item 3, `Q2 §4.4`) can push either way. So large `N` gives the `N⁰` scaffold
   and the leakage structure, not a certified finite-`N` floor. `[O]`.

5. **A no-go on the whole positivity/moment route at large `N` `[NG]`.** `g₁/Σg` is a **point mass**
   of a positive measure; **no** functional of its moments `M_k=∫s^k dν` (or of reflection positivity,
   PF, log-convexity) can lower-bound it — a measure can carry arbitrarily little mass at its lightest
   point while keeping every moment fixed (witness D: two measures, identical `(M₀,M₁,M₂)`, `g₁/Σg =
   0.40` vs `0.001`). This is the **large-`N`/spectral echo of `U1-NG1`** (source direction) and
   `R2 §2 Lemma C` (moments only **upper**-bound the isolation `Λ★`). The only object that can bound
   `g₁/Σg` is the **kernel `K`** (the master-field dynamics) — which is `[O]`. **Large `N` relocates
   the wall faithfully; it does not soften it** (agrees `U3 §5.2`, `KAPPA-LARGEN`).

**One line.** *The master-field equation for `ψ₁` is the closed linearised-MM / Bethe–Salpeter
eigen-equation `ℳψ₁=m₁²ψ₁` around the master field `W₀` `[T]` — solvable in 2D (where it returns `w₁=0`
for `tr F²`), unconstructed in 4D `[O]`; the large-`N` atom fraction is the point mass of the
normalised KL measure at its lightest glueball, `g₁/Σg=1/(1+R)` `[T]`, exactly `s/2` at strong coupling,
`0` for 2D `tr F²` and the free field, `O(N⁰)` by the `1/N²` cancellation; it is **not** uniformly
positive in the 't Hooft coupling — the free/UV end is `0` and the acid test *requires* it — so the
correct uniformity is fixed `β=2.5` in `(L,a)` `[O]`; and no moment/positivity functional can bound
it `[NG]`, because it is a point mass of a measure whose kernel (the master field) is the unbuilt
`H1`. **W2 delivers the `N`-scaffold and kills the wrong (λ-uniform) formulation; it does not prove
`W`.** No QED. Clay 0%.*

---

## 1 · Set-up: the large-`N` spectral frame, `ℳ`, `ψ₁`, `g₁`, `Σg`

Pure `SU(N)` YM, 't Hooft coupling `λ=g²N` fixed, `N→∞`. The gauge-invariant `0⁺⁺` **source** is a
single-trace density `O=(1/N)tr F²` (continuum) or `O=(1/N)Re\,tr U_p` (lattice plaquette), normalised so
`W(C)=⟨(1/N)tr P∮A⟩=O(1)`. The connected two-point function and its Källén–Lehmann (KL) representation:
```
   G(x) = ⟨O(x)O(0)⟩_c = ∫_{[0,∞)} dν(s) Δ(x;√s) ,    Δ the free scalar propagator of mass √s ,
   dν(s) = Σ_n g_n δ(s−m_n²) + dν_cut(s) ,     g_n = |⟨n|O|Ω⟩|² ,   f₁² := g₁ ,   M₀ := ∫dν = Var(O) .
```
The **two normalisations, kept distinct** (`U3 §1`, `C1-SZZ §4`): the physical `1/N`-normalised
`O=(1/N)tr F²` has `M₀≍1/N²`, so define the **reduced** residues `g_n := N²·|⟨n|O|Ω⟩|² = O(1)`; then
`f₁²=g₁/N²`, `M₀=(Σg_n)/N²`, and `w₁=f₁²/M₀=g₁/Σg_n` is `N`-free at leading order (`U3` Theorem U3-1).
Throughout below `g_n` are the reduced (`O(1)`) residues, and the atom fraction is
```
   w₁ = g₁/Σ_n g_n ∈ [0,1] ,      (U3-1, [T])
```
with the same `1/N²` source normalisation cancelling in the ratio.

**The three large-`N` objects of W2.**
- **`ℳ`** — the linearised Makeenko–Migdal operator around the master field `W₀` (the `0⁺⁺`, `k=0`
  restriction of the loop-space fluctuation operator; identical, by Giles' `spec(loops)=spec(fields)`,
  to the gauge-invariant field-fluctuation operator, `PNT-LARGEN-LOOP §4`).
- **`ψ₁`** — its lightest eigenvector (= the lightest `0⁺⁺` glueball wavefunction).
- **`g₁`** — the source projection `g₁=|⟨ψ₁|O|Ω⟩|²`; `Σg` — the total spectral weight `M₀` (reduced).

The **two structural facts** that make W2 tractable at all are U3's `[T]` results and the closure of `(MM∞)`
(`PNT-LARGEN-LOOP §1`): the ratio `w₁=g₁/Σg` is `N⁰` (algebraic), and the loop equation is **one closed
equation in `W(C)`** (no tower) — so the fluctuation spectrum `spec(ℳ)` is *defined* without an infinite
hierarchy. The open content is that `W₀` (hence `ℳ`, hence `ψ₁`, `m₁`, `g₁`) is **unconstructed in 4D**.

---

## 2 · The master-field equation for `ψ₁` `[T]` (structure) `+ [O]` (content)

### 2.1 Linearised Makeenko–Migdal (`MF`)

The large-`N` loop equation (`MM∞`, `PNT-LARGEN-LOOP §1`):
```
   ∂_μ (δ/δσ_{μν}(x)) W(C) = λ ∮_C dy_ν δ^{(4)}(x−y) W(C_{xy}) W(C_{yx}) .                          (MM∞)
```
Write `W = W₀ + h` around the vacuum master field `W₀` (the solution of `(MM∞)` with `W(∅)=1`). To first
order in `h`:
```
   ℳ h = 0 ,     ℳ := [ ∂_μ δ/δσ_{μν} − 2λ ∮ dy_ν δ(x−y) W₀(C_{xy})(·) ]  .                          (MF)
```
`ℳ` is a linear operator on loop-space functions, self-adjoint on the physical (reflection-positive,
gauge-invariant) `0⁺⁺` section (`Giles`). Its **spectrum is the glueball mass-squared spectrum**,
`spec(ℳ)={m_n²}`, and its eigenvectors `ψ_n` are the glueball wavefunctions; the lightest,
```
   ℳ ψ₁ = m₁² ψ₁ ,        ψ₁ = lightest 0⁺⁺ eigenvector ,        g₁ = |⟨ψ₁|O|Ω⟩|² ,                  (MF-ψ₁)
```
is *the master-field equation for `ψ₁`*. `[T]` for (MF)/(MF-ψ₁) as the exact leading-`N` eigen-equation
(this **is** the `(MM∞)` closure — no `n`-point tower). `[O]` for solving it, since it requires `W₀`.

### 2.2 The equivalent planar Bethe–Salpeter form

`(MF-ψ₁)` is the homogeneous **Bethe–Salpeter** equation of the gauge-invariant `0⁺⁺` channel. Let
`G(p²)` be the connected `O`-two-point function, `G₀` its bare (disconnected/one-loop) part, and `K`
the **planar** (leading-`N`), `0⁺⁺`, two-particle-irreducible kernel assembled from `W₀`'s vertices
(factorization kills the multi-trace exchanges at the leading order, so `K` is closed). Then
```
   G^{-1}(p²) = G₀^{-1}(p²) − K(p²) ,      glueball poles at det[G₀^{-1}(p²) − K(p²)] = 0 ,  p²=−m_n² ,
   [G₀^{-1}(−m₁²) − K(−m₁²)] ψ₁ = 0 ,      g₁ = (the residue of G at p²=−m₁²) .                       (BS)
```
In self-adjoint (physical) form `G₀^{1/2} K G₀^{1/2} χ = χ`, `χ = G₀^{-1/2}ψ₁`; the lightest root gives
`m₁²`, and `g₁` is the projection of the source onto it. `(BS)` is the standard large-`N` bound-state
equation; **given `K` (i.e. given `W₀`), it determines `m₁`, `ψ₁`, and `g₁` exactly.** Witness A
realises `(BS)`/`(MF-ψ₁)` on a self-adjoint toy `ℳ` (a discretised linearised-MM operator) and reads
`g₁/Σg` off the source overlap and the resolvent residue — confirming the *logic* `ℳ → ψ₁ → g₁ → g₁/Σg`
to machine precision.

### 2.3 What is solvable, and why 4D is `[O]`

- **2D `[T]`.** `W₀(C)=exp(−λA(C)/2)` exactly (Singer; `A(C)` the signed area, `PNT-LARGEN-LOOP §2`),
  so `ℳ` is explicit. The `tr F²` carrier is `tr E²=H`, conserved and diagonal in the character basis:
  `⟨0|tr F²|G_R⟩∝C₂(R)⟨χ₀|χ_R⟩=C₂(R)δ_{0R}=0` for `R≠0`. Hence **`g₁=0` at every `N`** — the `0⁺⁺`
  atom is *absent*, and `(MF-ψ₁)` correctly returns `w₁∞=0` (witness B2; the acid test). The
  **holonomy** carrier (source `=χ_f`, a single character) instead has `⟨χ_f|χ_R⟩=δ_{fR}`, so its
  source couples to **one** level: `g₁=Σg`, `w₁∞=1` — the `[0,1]` other endpoint.
- **4D `[O]`.** `W₀` is **unconstructed**: a 4D square boundary is spanned by the flat square *and* by
  tents of every apex height (`PNT-LARGEN-LOOP §3`), so `W(C)` is a *string functional* (Coulomb +
  area), not `exp(−σA)`; Singer did 2D only. Without `W₀` there is no `ℳ`, hence no `(MF-ψ₁)`. The one
  *explicit-model* 4D solve is the rigid-string completion (`CONSTRUCT-MASTERFIELD`): `ψ₁` = the lowest
  radial mode of the closed flux tube, `m₁=2π√(2σκ)` with rigidity `κ`; it gives `m₁` but **not** `g₁`
  (the source overlap needs the tube's coupling to `tr F²`, and `κ` itself is `H1`, `KAPPA-LARGEN`).
  So even the model localises the missing input to one number and does not compute the residue.

> **Net of §2.** `(MF-ψ₁)` / `(BS)` is **the** master-field equation for `ψ₁`: a closed, linear,
> self-adjoint eigen-equation whose lightest eigenpair `(m₁²,ψ₁)` fixes the atom. It is `[T]` as an
> equation (the `(MM∞)` closure), exact in 2D (where it gives the right `0` and `1`), and `[O]` in 4D
> because the master field `W₀` — the equation's only input — is unbuilt.

---

## 3 · Computing `g₁/Σg` `[T]` — the ratio, its reformulation, and its exact values

### 3.1 The ratio identity `(RATIO)`

> **Theorem W2-1 (`g₁/Σg` is the atom weight of the normalised KL measure) `[T]`.**
> With `dν=Σ_n g_nδ(s−m_n²)+dν_cut≥0`, `M₀=∫dν`, `f₁²=g₁=ν({m₁²})`, the reduced atom fraction is
> ```
>    w₁∞ = g₁/Σ_n g_n = 1/(1+R) ,   R=Σ_{n≥2}g_n/g₁ ,   = ν({m₁²})/ν([0,∞)) = ρ({m₁²}) .
> ```
> A positive floor `w₁∞≥w₀>0` is **equivalent** to a leakage **upper** bound `R≤(1−w₀)/w₀`, i.e.
> `Σ_{n≥2}w_n≤1−w₀` — identical to T2's missing estimate. ∎

*Proof.* Immediate from `w₁=g₁/Σg` and `R=(Σg−g₁)/g₁`; `1/(1+R)=g₁/Σg`. The `ρ`-form is the definition
`ρ=dν/M₀`. The leakage equivalence: `w₁≥w₀ ⟺ Σ_{n≥2}g_n≤(1−w₀)Σg ⟺ R≤(1−w₀)/w₀`. ∎

Witness A checks `1/(1+R)=g₁/Σg` to `1e-12` and `g₁`=`resolvent residue` to `1e-4`, on a toy `ℳ` with
`g₁/Σg=0.589991`. This is the **entire** analytic content of "computing `g₁/Σg`": it is not a new
number to be minimised but a **point mass** of a probability measure, and the whole problem is to show
that point mass is not small.

### 3.2 The exact values (the acid battery's arithmetic)

| carrier / limit | `g₁/Σg` | why | label |
|---|---|---|---|
| **strong coupling** (`β<1/48`) | `s/2` (`=1/2` for physical `s=1`) | `g₁=a₀=(2N)⁻ᵖs`, `M₀=1/(2N²)` at `β=0` (Haar `∫|trU|²=1`); `w₁≤1` forces `p=2` ⇒ `w₁=s/2`, `N`-uniform | **[T]** (mod `[KL]`) |
| **2D `tr F²`** | `0` | `tr F²=H` conserved; `⟨χ₀\|χ_R⟩=δ_{0R}` ⇒ `g₁=0` for `R≠0` | **[T]** |
| **free field** (`λ→0`) | `0` | composite `tr F²` is a bilinear ⇒ KL measure is the a.c. two-photon cut, **no delta**; `ν({m²})=0` | **[T]** |
| **2D holonomy** (`χ_f` source) | `1` | `⟨χ_f\|χ_R⟩=δ_{fR}` ⇒ source = one level | **[T]** |
| **general large `N`** | `1/(1+R)`, `O(N⁰)` | `(RATIO)` + U3-1; `R` finite, `N`-free at LO | **[T]** (algebraic) |
| **anomaly / dilaton** | `∝ m_G⁴ f_G² / χ_{trF²}` | `θ^μ_μ∝tr F²` ⇒ `⟨0\|trF²\|G⟩∝m_G² f_G`; single-dilaton dominance ⇒ `g₁∝m_G⁴ f_G²N²`; divide by susceptibility `Σg` | **[C]** |

The **anomaly row** `[C]` is the one genuinely *physical* identification: in pure YM the trace anomaly
is `θ^μ_μ∝tr F²`, so the `0⁺⁺` glueball is the "dilaton", and `⟨0|tr F²|G⟩∝m_G² f_G` (`f_G` the
dilatation decay constant — the standard trace-anomaly normalization, *assumed* single-dilaton
dominance). Then
```
   g₁/Σg =  (anomaly factor)² · m_G⁴ f_G² N² / (N² M₀)  =  (anomaly factor)² · m_G⁴ f_G² / (trF² susceptibility) ,
```
which **identifies** `g₁/Σg` with a ratio of physical glueball quantities (`m_G`, `f_G`, the `tr F²`
susceptibility) — but does **not** bound it, because `f_G` is as uncomputed as `g₁`. It is the honest
"master-field reads it off" statement: given the dilaton parameters, the ratio follows; it does not
give the parameters. `[C]`.

Witness B checks the four exact rows numerically: B1 (free/cut `g₁/Σg→0` under refinement: `0.078,
0.020, 0.004, 0.001` at `50,200,1000,4000` poles — the `1/n_poles` point-mass vanish), B2 (2D `tr F²`
identity `0.0e0`, SU(2) character orthogonality `6.2e-17`), B3 (`s/2` for all `N`, `s∈{1,2}`; `s>2`
excluded by `w₁≤1`), B4 (`g₁/Σg→g₁/Σ_{ST}g`, deviation `=O(1/N²)`, leakage `O(1/N²)→0`).

---

## 4 · Is `g₁/Σg` bounded below uniformly in the 't Hooft coupling? — **No**, and *three* senses of `λ` `[T]/[NG]/[O]`

The charter's W2 asks "uniform in the 't Hooft coupling." That phrase hides **three different
uniformities**, and only one is the `W` target. Separating them is the whole answer.

### 4.1 (a) RG-running of the *renormalised* coupling: `w₁` does not run `[T]`

The physical continuum statement: `f₁²` and `M₀` both carry the same `Z_O²` (one insertion of `tr F²`
each), so the ratio is **RG-invariant**: `γ_{w₁}=d ln w₁/d ln μ=2γ_O−2γ_O=0` (`U4`, `T`). So along a
line of constant physics, `w₁` is **`μ`-independent** — **no RG degradation** of the atom fraction.
This is a genuine `[T]` positive: the "uniform in the renormalised 't Hooft coupling" reading is *free* —
`w₁` is constant in it.

### 4.2 (b) Bare-coupling sweep `λ∈[0,∞)` at fixed lattice: **no positive floor** `[T]/[NG]`

If one varies the *bare* 't Hooft coupling (equivalently the lattice `β`) over its whole range, the
family includes the **free point `λ=0`**, where by the acid test
```
   g₁/Σg(λ=0) = 0     EXACTLY (free Maxwell: a.c. two-photon cut, no atom)      [T] ,
   g₁/Σg(λ)  → s/2    as λ→∞ (β→0, strong coupling)                              [T] mod [KL] ,
   g₁/Σg(λ)  ↗  from 0 to s/2 across λ  (physical smooth/monotone interpolation)  [C] ,
```
so `inf_{λ≥0} g₁/Σg = 0`. **There is NO positive lower bound uniform in `λ` over a bare range that
contains `λ=0`.** This is **forced and rigorous at the endpoint**: `λ=0` *is* the free field, and the
free-field acid test demands `w₁=0` there — a floor claimed over all `λ` would be **vacuous on the free
field**, precisely the acid failure the charter forbids. (The `↗` shape — continuity and monotonicity on
`λ>0` — is the physical reading `[C]`, *not* used for the verdict; the verdict needs only the exact
`λ=0` value.)
Witness C realises the interpolation `g₁/Σg(λ)=(s/2)·λ/(λ+λ₀)` with `inf_λ=0` (table `0.000000,
0.000500, …, 0.499500, 0.500000`). **This negative answer is the correct, acid-test-consistent one**; it
is the large-`N`/spectral face of `U1-NG2` (free-field parity/number, `w₁=0`) and `R2 §3` (the
continuum may abut `m₁²`).

> **Monotonicity conjecture `[C]` (not used above).** Physically `g₁/Σg(λ)` should **increase** with `λ`
> (stronger coupling ⇒ the source's weight localises more on the lightest glueball; the weak/UV end
> spreads it over the two-gluon cut). If `g₁/Σg` were proven monotone increasing in `λ`, a floor at any
> `λ₀` would extend to all `λ≥λ₀`, reducing "uniform in `λ≥λ_phys`" to the value at `λ_phys`. **This is
> unproven** and is *not* needed for the verdict — it would only convert (b)'s `inf=0` into "0 is
> attained only at the free/UV end," which is already clear.

### 4.3 (c) Fixed *physical* coupling (`β=2.5`), uniform in `(L,a)`: the `W` target `[O]`

The uniformity `W` actually needs is **at the physical point**: `∃w₀>0` with `g₁/Σg(β=2.5;L,a)≥w₀` for
all `(L,a)` through the physical lattice. This is **not** the bare-`λ` sweep (4.2) and **not** the
RG-running (4.1, already free); it is a **fixed-`λ`, fixed-physics, `(L,a)`-uniform** bound. It is
**exactly the `W` target** and is **`[O]`** — same object as `T2 §8` / `U1 (E-U1)`. At `β=2.5` the `0⁺⁺`
glueball is **not yet resolved** (`M1-4D-RUN §0.5`: the connected plaquette correlator is flat in `t`,
`m₁` and `f₁²` not extractable at `L≤8`), so even the *value* `g₁/Σg(2.5)` is `[V]/[O]` (the `T2`
prediction `w₀=O(0.2)–O(0.5)`; physical `f₁²_bare=7.14` at `β=2.2`, `R5`).

> **Verdict of §4 (three senses, three labels).** (a) RG-running: `w₁` **does not run** `[T]` (free).
> (b) bare-`λ` sweep: **no uniform floor**, `inf=0` `[T]/[NG]` (the free-field acid *requires* it).
> (c) fixed `β=2.5`, `(L,a)`-uniform: the `W` target `[O]`. **W2 answers the charter's question
> precisely: "uniform in the 't Hooft coupling" in the bare sense is FALSE; the atom fraction's
> `λ`-structure is `0`(free)↔`s/2`(strong), and the uniformity that matters is (c).** Large `N` does
> **not** prove `W`; it proves the `N⁰` scaffold and shows the naive λ-uniform formulation is
> acid-forbidden.

---

## 5 · The `1/N²` correction: positivity `[T]`, the floor `[O]`

Write the finite-`N` atom fraction as its large-`N` value plus the leading correction:
```
   w₁(N) = g₁/Σg + c/N² + O(1/N⁴) ,                                                              (1/N²)
```
where `c` collects the multi-trace leakage and the connected-`4`-point correction (`U3 §5.2` item 3,
`Q2 §4.4`; the `Ĝ` term in `(MM)`), and the expansion is in `1/N²` (orientable surfaces; `KAPPA-LARGEN
§2`). Two separate questions:

- **Does it preserve *positivity* `w₁(N)>0`? `[T]`** — **Yes, trivially, at every `N`.** `w₁=f₁²/M₀` is a
  ratio of positive numbers with `M₀≥f₁²` (the total spectral weight dominates the lightest atom), so
  `0≤w₁(N)≤1` *by construction* at every finite `N`, independent of the correction's sign. Witness E:
  `w₁(N)∈[0,1]` for `c=±0.30` and all `N` tested. **Positivity is not the issue.**
- **Does it preserve the *floor* `w₁(N)≥w₀`? `[O]`** — **No, not from the large-`N` value alone.**
  Two gaps: (i) the asymptotic `1/N²` series has **no uniform remainder at `N=3`** — `PROOF-LARGE-N-RIGOROUS`
  (c): `|c|/N²=O(0.1)` at `N=3`, witness E (`0.033`); nothing bounds it (a gap can vanish at `∞` and
  reappear at `N=3`); (ii) the **sign** of `c` is **not fixed** by factorization — the leakage can push
  either way (`U3 §5.2`). So the large-`N` floor does **not** imply a finite-`N` floor.

**At the physical `N`** one must work directly at that `N` (the `1/N²` expansion is not a bound), and
then positivity of `w₁` is the *same* `[O]` as `W` itself — no large-`N` shortcut. This is the large-`N`
restatement of `PROOF-LARGE-N-RIGOROUS (c)` ("transfer `N=∞→3` has no theorem") and `U3`'s
"large-`N` transports the binding problem faithfully; it does not reduce it."

---

## 6 · The no-go: `g₁/Σg` is a point mass — no moment / positivity functional bounds it `[NG]`

The decisive structural obstruction, at large `N`, is the **same** as U1/R2 at finite `N`, now in the
spectral coordinate:

> **No-go W2-NG (moment blindness of the atom weight) `[NG]`.** `g₁/Σg=ρ({m₁²})` is the mass the
> normalised KL measure `ρ` puts on its lightest **point**. For **any** `K`, the moments
> `(M₀,…,M_K)` of `ρ` do **not** determine `ρ({m₁²})`: there are pairs of positive measures with
> identical first `K` moments and atom weights differing by an arbitrary factor. Hence **no functional
> of the moments** (equivalently of reflection positivity / PF / log-convexity, which are moment
> properties) can lower-bound `g₁/Σg`.

*Witness.* Two probability measures on `s∈[0,∞)` with **(identical)** `(M₀,M₁,M₂)=(1,\,5.8,\,49.0)`:
- `ρ_A` = atoms `s=1 (0.4)`, `s=9 (0.6)` ⇒ `g₁/Σg=0.400`;
- `ρ_B` = atoms `s=1 (0.001)`, `2.5 (0.4976)`, `6 (0.2436)`, `12 (0.2578)` ⇒ `g₁/Σg=0.001`.

Same first three moments (to `1e-9`), atom weights differing by `400×` (witness D). So any bound built
from `M₀,M₁,M₂` — the total weight, the mean mass², the second moment — is **blind** to `g₁/Σg`. This
is `R2 §2 Lemma C` ("moments only **upper**-bound the isolation `Λ★`") and `U1-NG1` ("`w₁` is the
source direction cosine; positivity + gap do not bound it") transported to `N=∞`. **The only object that
can bound `g₁/Σg` is the kernel `K` of `(BS)` — the master-field dynamics — which is `[O]` in 4D.**
This is why large `N` **relocates** the wall (to `spec(ℳ)`, `KAPPA-LARGEN`'s "the glueball residue /
mass ratio = `H1` in the `N`-coordinate") rather than softening it.

---

## 7 · The acid-test battery (charter §2) — all pass

| # | test | statement | `g₁/Σg` | label |
|---|---|---|---|---|
| 1 | **free massless** | no atom (a.c. two-photon cut, bilinear source) | **0** | **[T]** (witness B1) |
| 2 | **2D `tr F²`** | `tr F²=H` conserved, `⟨χ₀\|χ_R⟩=δ_{0R}` | **0** | **[T]** (witness B2) |
| 3 | **strong coupling** | `g₁=a₀`, `M₀=1/(2N²)` ⇒ `w₁=s/2`, `N`-uniform | **`s/2`** | **[T]** mod `[KL]` (witness B3) |
| 4 | **large `N`** | `1/N²` source normalisation cancels; leakage `1/N²` | **`O(N⁰)`** | **[T]** (witness B4) |
| 5 | **continuum** | `γ_{w₁}=0` (`U4`); `g₁/Σg` is `μ`-invariant | (no run) | **[T]** |

The four required acid tests (free `0`, 2D `tr F²` `0`, strong `s/2`, large-`N` `O(N⁰)`) all pass, and
the free-field `0` is exactly what forces §4.2's *negative* answer to "λ-uniform." **No acid test is
manufactured positive; the free field and 2D correctly give `0`.**

---

## 8 · Witness — `research/v3/w7/w2_master_field.py`

`py -3.12 -I w2_master_field.py` (numpy; deterministic `seed 0`; exit 0; **18/18 PASS**). Quoted
output (abridged):

```
PART A  master-field operator -> g1/sum g  (spectral theorem logic)
    lightest m1^2 = 1.000000  (isolated: gap to next = 1.500)
    g (KL weights, lightest first) = [0.589991 0.121379 0.110038 0.085639 0.008062 0.084891]
    g1 = 0.589991 ; sum g = 1.000000 ; R = 0.694941
    [PASS] A1 g1/sum g = nu({m1^2})/nu(total) = 1/(1+R)  (0.589991)
    [PASS] A2 g1 = lightest-pole residue of the source resolvent
    [PASS] A3 psi1 is the lightest eigenvector; g1 = |<psi1|s>|^2
PART B  acid battery
    B1 free/cut: g1/sum g at npoles=50,200,1000,4000 = [0.078318 0.020171 0.004066 0.001018]
    [PASS] B1 free field (a.c. cut, no atom): g1/sum g -> 0
    B2 char orthogonality max|int chi_j chi_k - delta| = 6.18e-17 ; trF^2 fraction = 0.0e+00
    [PASS] B2 2D trF^2: atom weight 0 (conserved density), g1/sum g = 0 (acid PASS)
      s=1: w1 = s/2 = 0.5 for all N (N-uniform) ;  s=2: w1 = s/2 = 1.0 for all N
    [PASS] B3a strong coupling: g1/sum g = s/2, O(N^0), N-uniform (s=1,2)
    [PASS] B3b s>2 excluded by w1<=1 (w1=s/2 must be a fraction)
    B4 large-N: g1/sum g(N) at N=2,8,1e6 = [0.36036036 0.51159073 0.52631579] ; target 0.526316
    [PASS] B4a large-N: g1/sum g -> g1/sum_ST g = O(N^0), finite>0
    [PASS] B4b leakage (MT+ac) is O(1/N^2) -> 0
    [PASS] B4c deviation = O(1/N^2) exactly (w1 is O(N^0))
PART C  lambda-uniformity: g1/sum g has NO positive lower bound over the bare coupling
    lambda      g1/sum g = (s/2)*lambda/(lambda+lambda0)
        1.0e-06   0.000000 ;  1.0e-01   0.045455 ;  1.0e+01   0.454545 ;  1.0e+06   0.500000
    [PASS] C1 g1/sum g -> 0 at the free/UV end (lambda->0)
    [PASS] C2 g1/sum g -> s/2 at strong coupling (lambda->inf)
    [PASS] C3 NO positive lower bound uniform over the bare coupling: inf_lambda g1/sum g = 0
PART D  moment-blindness no-go: same (M0,M1,M2), different g1/sum g
    measure 1: atoms(s=1:0.4, s=9:0.6) -> (M0,M1,M2)=[1. 5.8 49.] ; g1/sum g = 0.400000
    measure 2: atoms(s=1:1e-3, 2.5, 6, 12) -> (M0,M1,M2)=[1. 5.8 49.] ; g1/sum g = 0.001000
    [PASS] D1 both measures match the first 3 moments to 1e-9
    [PASS] D2 g1/sum g differs by >100x at identical moments (moment-blindness no-go)
PART E  1/N^2 correction: w1 stays in [0,1]; no uniform N=3 remainder
    [PASS] E c=+0.30 / c=-0.30: w1(N) in [0,1] for all N (fraction, positivity exact at every N)
    [PASS] E3 |c|/N^2 at N=3 is O(0.1): asymptotic series gives NO uniform floor at N=3
# 18/18 checks PASS
```

Readings: **(A)** a self-adjoint toy `ℳ` (the *structure* of the linearised-MM operator) reproduces
`g₁/Σg=1/(1+R)`=`resolvent residue`=`|⟨ψ₁|s⟩|²` to machine precision — the *logic* `ℳ→ψ₁→g₁`; in 4D
`ℳ` (the operator around the unbuilt `W₀`) is the `[O]`. **(B)** the four acid tests. **(C)** the
λ-uniformity **fails** (`inf=0` at the free/UV end), as the acid test forces. **(D)** moment blindness:
identical `(M₀,M₁,M₂)`, atom weights `0.40` vs `0.001`. **(E)** positivity exact at every `N`, but no
uniform `N=3` remainder. **Witness, not proof; the theorems are §2–§6, the numerics illustrate the
identity structure, the acid limits, and the two no-gos.**

---

## 9 · Ledger

| item | statement | label |
|---|---|---|
| master-field eq. for `ψ₁` | `ℳψ₁=m₁²ψ₁`, `ℳ=δ(MM∞)/δW\|_{W₀}` on `0⁺⁺`; `=` planar BS `[G₀⁻¹−K]ψ₁=0` — **closed** (no tower) | **[T]** structure |
| solving the eq. in 4D | needs `W₀`; **unconstructed** (surface non-uniqueness, Singer 2D only) | **[O]** |
| 2D solve | `W₀=e^{−λA/2}`; `tr F²` carrier `g₁=0`, holonomy carrier `g₁=Σg` | **[T]** |
| `g₁` | `g₁=\|⟨ψ₁\|O\|Ω⟩\|²` = lightest-pole residue of the source resolvent | **[T]** |
| `g₁/Σg` identity | `=1/(1+R)`, `R=Σ_{n≥2}g_n/g₁`; `=ν({m₁²})/ν([0,∞))` = atom weight of normalised KL measure | **[T]** |
| floor ⟺ leakage | `g₁/Σg≥w₀ ⟺ R≤(1−w₀)/w₀ ⟺ Σ_{n≥2}w_n≤1−w₀` (= T2's estimate) | **[T]** |
| strong-coupling value | `g₁/Σg=s/2` (`=1/2` physical `s=1`), `N`-uniform; `s>2` excluded by `w₁≤1` | **[T]** mod `[KL]` |
| 2D `tr F²` value | `g₁/Σg=0` exactly, all `N` (conserved density) — acid PASS | **[T]** |
| free-field value | `g₁/Σg=0` exactly (a.c. bilinear cut, no delta) — acid PASS | **[T]** |
| large-`N` scaling | `g₁/Σg=O(N⁰)`, leak `Σ_{n≥2}w_n=O(1/N²)` | **[T]** algebraic |
| anomaly / dilaton | `⟨0\|trF²\|G⟩∝m_G²f_G` ⇒ `g₁/Σg∝m_G⁴f_G²/χ_{trF²}` (identifies, does not bound) | **[C]** |
| λ-uniformity (bare sweep) | `g₁/Σg(λ=0)=0` **exactly** (free, acid) `[T]`; `→s/2` strong `[T]`; `inf=0` over any range containing `λ=0` ⇒ **NO** λ-uniform floor | **[T]/[NG]** |
| λ-uniformity (RG-running) | `γ_{w₁}=0`: `w₁` does not run with `μ` | **[T]** (`U4`) |
| λ-uniformity (fixed `β=2.5`) | `(L,a)`-uniform floor = the `W` target; value unresolved at `β=2.5` | **[O]** |
| monotonicity in `λ` | `g₁/Σg` expected increasing in `λ` (unproven; not used) | **[C]** |
| `1/N²` positivity | `w₁(N)∈[0,1]` at every `N` (it is a fraction) | **[T]** |
| `1/N²` floor | not inherited: no uniform `N=3` remainder (`\|c\|/N²=O(0.1)`), sign of `c` unfixed | **[O]** |
| moment blindness (no-go) | `g₁/Σg` is a point mass; identical `(M₀,M₁,M₂)` can carry `g₁/Σg` differing `400×` ⇒ no moment/RP functional bounds it | **[NG]** |
| what large `N` gives | `N⁰` scaffold + `1/N²`-leakage structure; value = glueball residue ratio = `H1` in `N`-coordinate | — |
| **Clay mass gap** | unbuilt | **[O] · 0%** |

---

## 10 · Sharpest line

**The master-field equation for the first excited eigenvector is the closed linearised-Makeenko–Migdal
(planar Bethe–Salpeter) eigen-equation `ℳψ₁=m₁²ψ₁` around the master field `W₀` — an honest `[T]` as an
equation (the `(MM∞)` closure, no `n`-point tower), exact in 2D (where it returns the right `g₁/Σg=0` for
`tr F²` and `1` for the holonomy), and `[O]` in 4D because its only input, `W₀`, is unconstructed — and
the large-`N` atom fraction it determines is the *point mass of the normalised Källén–Lehmann measure at
its lightest glueball*, `g₁/Σg=1/(1+R)=ρ({m₁²})`, exactly `s/2` at strong coupling, `0` for 2D `tr F²`
and the free field, `O(N⁰)` by the cancellation of the `1/N²` source normalisation.** It is **not**
bounded below uniformly in the 't Hooft coupling — the free/UV end is `0` (`λ→0`), and the free-field
acid test *requires* that, so any λ-uniform floor would be vacuous where it must fail; `w₁` does not run
under RG (`γ=0`), and the uniformity that `W` actually needs is *fixed physical `β=2.5`, uniform in
`(L,a)`*, which is the `[O]` target. The `1/N²` correction keeps `w₁` in `[0,1]` at every finite `N`
(it is a fraction) but inherits no uniform floor (no `N=3` remainder; sign unfixed). And no
moment/positivity functional can bound `g₁/Σg`, because it is a point mass of a measure whose *kernel* —
the master field — is the wall: identical first three moments carry atom weights `0.40` vs `0.001`.
**Large `N` delivers the `N`-scaffold and the `1/N²`-leakage structure of a uniform binding theorem,
and it kills the wrong (λ-uniform) formulation of `W`; it does not supply the value, which is the
glueball residue ratio `H1` in the `N`-coordinate.** No QED. Clay 0%.

---

*Refs (internal).* `W-WAVE-CHARTER.md` §1 W2, §2 acid battery, §3 what W does/does not prove;
`U3-LARGE-N-BINDING.md` (Theorem U3-1, §2–§5), `U1-HARNACK-BOUND.md` (master identity, NG1/NG2, `(E-U1)`),
`U-WAVE-SYNTHESIS.md` §1,§6; `U4-CONTINUUM-ATOM.md` (`γ_{w₁}=0`); `T2-TUBE-OVERLAP.md` §1.3,§2.4,§8,
`T-WAVE-SYNTHESIS.md` §5; `R2-ISOLATION.md` §2 (Lemma C: moments upper-bound `Λ★`), `T3-TP2-ISOLATION.md`;
`unconditional/PNT-LARGEN-LOOP.md` §1–§4 (`(MM∞)` closed; `spec(ℳ)`; `W₀` 4D-unconstructed),
`CONSTRUCT-MASTERFIELD.md` §1–§4 (rigid-string `m_{0⁺⁺}=2π√(2σκ)`; `κ=H1`), `KAPPA-LARGEN.md` §2–§5
(`κ(N)=κ_∞+c/N²`; value = gap), `ROOF-LARGEN-MASTER.md` §1,§3,§4 (renormalon-free mass; 4D non-closure in
loop space; master field open), `G1-LARGE-N.md` (large-`N` β<0 but nonperturbative `(G1)` = wall),
`PROOF-LARGE-N-RIGOROUS.md` (a),(c) (strong-coupling master field; **no `N=∞→3` transfer theorem**);
`v3/M1-4D-RUN.md` §0.5 (`0⁺⁺` unresolved at `β=2.5`, `L≤8`); `R5-EXACT-ANCHORS.md` §5 (`β=2.2`
`f₁²=7.14`); `Q2-LOOP-BOOTSTRAP.md` §2.2,§4.4 (`Ĝ` connected-4-point = leakage). *External, standard,
from the corpus (not re-fetched):* 't Hooft, NPB 72 (1974) 461 (planar counting); Witten, in *Recent
Developments in Gauge Theories* (1980) (master field); Makeenko–Migdal, PLB 88 (1979) 135; Migdal, Sov.
Phys. JETP 42 (1976) 413; Singer, "On the master field in two dimensions" (2D construction); Giles, PRD 24
(1981) 2160 (`spec(loops)=spec(fields)`); Källén (1952)/Lehmann (1954) (spectral representation); the QCD
trace anomaly / dilaton coupling (Scale anomaly: Adler–Bardeen–Bell–Jackiw lineage; Veneziano, "Anomalies
in QCD"; single-dilaton dominance is the standard saturation assumption `[C]`). **Witness:**
`research/v3/w7/w2_master_field.py` (`py -3.12 -I`, numpy only, deterministic, **18/18 PASS**).
**No third-party downloads; no `[C]` upgraded to `[T]`; no QED; Clay 0%.**
