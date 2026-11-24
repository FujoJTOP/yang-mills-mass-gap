# P3-POLYAKOV-CENTER — the constant temporal mode is the Polyakov holonomy; the Z(N) center sector does **not** change the physical periodic floor (it is the Haar-carried P2 value ≈0.222, bracketed by P1's exact 0.2051); the ~0.48 twist floor is an explicitly centre-broken 't Hooft-flux boundary, not the periodic value

> **2026-11-24 · P3 · the Polyakov / centre-sector analysis of the periodic constant-mode curvature.**
> Reads: `D-WAVE-SYNTHESIS.md`, `D5-UNIVERSALITY.md`, `D3-TRANSFER-CONSTANT.md`,
> `SM4-CONSTANT-MODE.md`, `PROOF-GAUSS-INVARIANT-PF.md`, `docs/paper-v2/ch04-gap.tex`,
> `docs/paper-v3/main.tex §6.8`, `research/probe-center/REPORT.md`.
> **Anchors (now landed):**
> **P1** — on the periodic torus the **global constant-holonomy slice has `S_W ≡ 0`** (every
> plaquette trivial: a constant commuting field is flat), so the measure there is **pure Haar**
> and `G″ = 0.2051` exact (single-link Haar marginal; the L-fold Haar convolution → `1/Var(ξ₃)=0.2688`).
> **P2** — the periodic `k=0` mode **saturates**: `c = 0.22196`, `q = 0.849`, `Var(a) → 4.45`
> (⇒ `1/Var = 0.2247`), a positive L-independent floor just above the Haar value.
> Witness: `research/v3/w7/p3_polyakov_center.py` (`py -3.12 -I [quick|full]`, numpy + the exact
> pure-Cartan transfer matrix of `d3_transfer.py`; deterministic, no Monte Carlo, **no third-party
> sources**; output quoted below).
>
> **Labels (charter §2).** **[T]** true theorem / exact identity · **[V]** numerical witness ·
> **[C]** conditional · **[O]** exact missing step · **[NG]** no-go for a named class.
> **No QED. Clay 0%.**

---

## 0 · Verdict (no QED; Clay 0%)

**The Z(N) centre sector does *not* change the periodic constant-mode floor.** Three exact/verified
reasons, in order of strength:

1. **[T] The constant-holonomy slice is flat in *every* centre sector.** On the global
   constant-holonomy slice (all links equal up to a centre factor, a constant commuting gauge field)
   **every plaquette is the identity**, so `S_W ≡ 0` — and multiplying one direction by the centre
   element `z∈Z(N)` keeps the field commuting and the plaquettes trivial. The measure on the slice is
   therefore **pure Haar in every centre sector**, and the constant-mode curvature is the Haar value
   `0.2051` in every sector (`[T]`; numerically `all |½tr P_plaq − 1| = 0`, `S_W=0`, for holonomy
   `I`, `−I`, `(−I,−I)`). This is P1, and P1 is **centre-blind by construction**.

2. **[T]+[V] The periodic measure is exactly centre-independent (`θ`-independent).** Haar is
   right-invariant: with `B=e^{iθσ₃/2}` and `U` Haar, `UB` is Haar for every `B`, so the L=1
   minimal-torus partition function and the plaquette distribution are **exactly** `θ`-independent
   (numerically `∫ρ_long da = 26.0980322` for `θ∈{0,π/2,π,3π/2,2π}`, spread `2.7×10⁻⁸`). The
   flat-chart `G″(a)` `θ`-spread (0.730 at `B=I` vs 0.013 at `B=−I`, L=1) is a **chart artifact**:
   the exponential chart is centred at `U=I` and `J(2π)=0` at the `−I` edge. It is not a change of the
   physical floor.

3. **[T] Finite-L centre symmetry is unbroken.** The 1D strip transfer matrix is gapped
   (`Δ = −log(λ₁/λ₀)>0` at `θ=0,π,2π`: `0.644, 3.176, 1.096` at `β=2.5`), and a gapped finite-range
   1D kernel has **no spontaneous symmetry breaking**; at `θ∈{0,2π}` the marginal is symmetric under
   `a→−a`, so `⟨tr P⟩=0`. The physical ensemble is the **centre-neutral** one — no sector is selected.

**Consequence for the question asked.** The physical periodic floor is **the Haar-carried P2 value
≈0.222** (bracketed below by P1's exact `0.2051`, above by the L→∞ Haar CLT `0.2688`), **not** a
centre-sector floor and **not** the `≈0.48` twist value. The fixed-boundary twist of D5/D3 (`θ=π`) is
an **actionful 't Hooft-flux sector** (`B=iσ₃`, `B²=−I`, twist plaquette `=−I`, `S_W/twist-plaquette
= 2β = 5.0 > 0`), i.e. an **explicitly centre-broken boundary**, not a flat Polyakov sector; its
higher floor is a boundary artifact. Decompactification in 4D is unchanged `[O] = H1`; the centre
sector neither rescues nor dooms it, because the floor it controls is centre-neutral.

**One line.** *The constant temporal mode is the Polyakov holonomy. The centre symmetry is a symmetry
of the periodic measure; the constant-holonomy slice is flat (`S_W≡0`) in every sector, and the
periodic measure is θ-independent and centre-unbroken at these L. So the centre sector does not change
the floor: the physical periodic floor is the Haar-carried P2 value ≈0.222 (P1: 0.2051 exact), not the
explicitly-broken twist floor ≈0.48.* **No QED. Clay 0%.**

---

## 1 · Identification: the constant temporal mode **is** the Polyakov / holonomy

On a periodic torus the temporal zero mode of `A₀` **is** the holonomy — the Polyakov loop. In the
corpus's flat exponential chart `U_l=e^{iξ_l·σ/2}` the constant mode is the uniform Cartan zero mode
```
   v_l = e₃/√L ,   a = ⟨v,ξ⟩ = (1/√L) Σ_l ξ_{l,3} ,
```
which is the Cartan/temporal component of the log of the holonomy `∏_l U_l` (SM4 §1 calls this the
"Cartan / holonomy class"). The boundary holonomy `B = e^{iθσ₃/2}` is the **spatial** holonomy (a fixed
spatial twist / 't Hooft flux); `a` and `B` are the two holonomy variables of the strip.

### 1.1 The Z(N) centre structure

- **Centre group.** `Z(N) = {e^{2πik/N}·I} ⊂ SU(N)` (for SU(2): `{I, −I}`). It acts by multiplying
  **all temporal links on one time slice** by `z∈Z(N)`; each temporal plaquette contains one temporal
  link `U` and one `U†`, so it picks up `z·z̄ = 1` — the Wilson action is **Z(N)-invariant** at every
  `β`, every volume, every boundary.
- **Charge.** The Polyakov loop `P(x)=tr∏_t U_t(x) → z·P(x)` carries **N-ality 1** (centre-charged).
  The constant mode `a` (its holonomy log) is therefore defined only **up to the centre orbit**. All
  gauge-invariant local operators — `tr F²`, glueballs, the vacuum, `0⁺⁺` flux — have **N-ality 0**
  (centre-neutral): they are **exactly blind** to the centre sector (`research/probe-center/REPORT.md`
  §4: centre projection is identically trivial in the adjoint, `max|W_A−1|=0`).
- **Centre-neutral sector.** The physical content is the centre-orbit sum
  `ρ_neutral(a) ∝ Σ_{z∈Z(N)} ρ(z·a)`; the constant-mode curvature that any gauge-invariant observable
  can see is the curvature of `ρ_neutral`. The question "does the centre sector change the floor?" is
  precisely "does `G″[ρ_neutral]` differ from `G″[ρ]` for a definite sector?"

---

## 2 · Does the centre sector change the floor? **No.**

### 2.1 [T] P1's flat slice is flat in every centre sector

The global constant-holonomy slice is a constant commuting gauge field. Every plaquette is a
commutator `U_μU_νU_μ†U_ν† = I`; multiplying one direction by the centre element `z` keeps all links
commuting (the centre is central) and keeps every plaquette `= I`. So `S_W ≡ 0` on the slice **for
every centre sector**, the measure is **pure Haar**, and `G″ = 0.2051` (single-link Haar marginal;
the L-fold Haar convolution → `1/Var(ξ₃)=0.2688` per D3/D4). The centre sector cannot change a
quantity that does not depend on it. Witness `[A]`:
```
   sector                          |½tr P_plaq − 1|     S_W/L²
   holonomy = I    (trivial)            0.00e+00        0.00e+00
   holonomy = −I   (centre)             0.00e+00        0.00e+00
   holonomy = −I × −I                   0.00e+00        0.00e+00     -> S_W = 0, every sector
```

### 2.2 [T]+[V] The periodic measure is exactly θ-independent (right-Haar)

With `B` fixed and `U` Haar, `P=UB` is Haar (right-invariance), so the L=1 partition function and the
plaquette distribution are **exactly** `θ`-independent. Witness `[B]`:
```
   θ = 0, π/2, π, 3π/2, 2π :   ∫ρ_long da = 26.0980322   (spread 2.7e-08)
```
The flat-chart `G″(a)` *does* vary with `θ` (0.730 at `B=I` vs 0.013 at `B=−I` at L=1), but this is a
**chart artifact**: `a=ξ₃` is a chart coordinate (centred at `U=I`, with `J(2π)=0` at the `−I` edge),
not a centre-covariant variable. A chart reparametrisation cannot change a physical floor; the
chart-invariant statements are `R=G″·Var→1` and `χ<∞` (D5), which are centre-neutral.

### 2.3 [V] A *fixed non-periodic* boundary does move the value — but that is explicit breaking

At a fixed boundary holonomy the sector value does move (exact pure-Cartan transfer matrix, `β=2.5`;
`d3_transfer.py` calibrated against full-3D at L=1 to `O(0.15)` on the action part). Witness `[C]`:
```
    L |  θ=0  G(a0*) |  θ=π  G(a0*) | θ=2π  G(a0*)
    4 |  0.7440      |  0.7029      |  0.4045
    5 |  0.5553      |  0.3686      |  0.4807
    6 |  0.4577      |  0.4665      |  0.4074
    7 |  0.4012      |  0.4579      |  0.0824      (θ=2π a0* wanders: flat/double-well chart)
    8 |  0.3653      |  0.4383      |  0.3887
```
The value moves by `~0.2–0.4` across sectors but stays **positive** everywhere; `θ=0` decreases toward
the Haar floor (`→~0.22`, matching D3's `c=+0.22±0.03` and P2's `0.22196`), while `θ=π` sits higher
at `≈0.44–0.47` (D5's `≈0.48`). **All these are explicitly broken boundaries.** In particular `θ=π`
is an **actionful 't Hooft-flux** sector: `B=iσ₃` with `B²=−I`, twist plaquette `=−I`,
`S_W/twist-plaquette = β(1−½tr(−I)) = 2β = 5.0 > 0` (witness `[A]`). It is *not* a flat Polyakov
sector; its higher floor is a boundary-flux artifact.

### 2.4 [T] Finite-L centre symmetry is unbroken — no selection, no SSB

The 1D strip transfer matrix is gapped (`[T]`, D3): `Δ=0.644 (θ=0), 3.176 (θ=π), 1.096 (θ=2π)` at
`β=2.5`, `ξ<∞`. A gapped finite-range 1D Gibbs kernel has **no spontaneous symmetry breaking**, and at
`θ∈{0,2π}` the marginal is symmetric under `a→−a`, so `⟨tr P⟩=0`: the Polyakov order parameter
vanishes. The physical ensemble is **centre-neutral** at these `L` — every sector enters with equal
weight. (In 4D the centre *does* break at high temperature / large volume; but on the flat slice
`S_W≡0` in the chosen sector too, so the floor there is still Haar `0.2051` — even a broken centre
does not move the constant-holonomy floor.) Centre-sector SSB would affect decompactification only if
it re-pinned the constant-mode measure; it does not, because that measure is Haar-carried.

---

## 3 · Which floor is the physical periodic one? Haar / P1 / P2 — **not** a centre-sector floor

| candidate floor | value | status |
|---|---|---|
| **P1** constant-holonomy slice, pure Haar (L=1 marginal) | **0.20510** | **[T]** exact; centre-blind (§2.1) |
| Haar CLT, L→∞ `1/Var(ξ₃)=9/(4π²−6)` | 0.26883 | **[T]** exact |
| **P2** periodic `k=0` saturation `c`; `1/Var(a)`, `Var→4.45` | **0.22196** (0.2247) | **[V]** |
| D5/D3 fixed-boundary twist `θ=π` | ≈0.44–0.48 | **[V]** explicitly-broken 't Hooft flux |

**Answer.** The physical periodic floor is the **Haar-carried P2 value ≈0.222**, bracketed by P1's
exact `0.2051` (single-link Haar) and the L→∞ Haar CLT `0.2688`. They agree because they are the same
object — a measure/entropic curvature of the (centre-neutral) holonomy marginal — at different L/action
order: P1 is the flat-slice Haar value, P2 adds the nonzero `S_W`-off-slice fluctuation variance
(`Var 3.72→4.45`, softening `0.2688→0.2247`). **No centre-sector floor exists:** the fixed-boundary
values (0.74→0.365 at `θ=0`; ≈0.48 at `θ=π`) are explicit-breaking boundary artifacts, and the `θ=π`
one is an actionful 't Hooft flux, not a flat Polyakov sector. In 4D the statement `χ<∞` — that this
centre-neutral floor stays positive as `L→∞` — is still **`H1` `[O]`**; P3 does not change it.

---

## 4 · Acid test (free field) — passes

- **Free massless Gaussian** (`Φ≡0`, `S=κ|ξ|²/2`): `ρ(a)` Gaussian ⇒ `G″=κ→0`, `Var=1/κ`. The
  **centre-orbit average of a Gaussian is the same Gaussian** (a Gaussian is centre-symmetric), so
  `G″` stays `κ→0`: the centre projection manufactures **no** floor. `[T]`
- **β=0 pure Haar**: `G″=0.2051` at L=1 in **every** centre sector (pure Haar is centre-blind), rising
  to `0.2688` at large L. The discriminator is the Haar Jacobian `Φ=−log J` (absent for the free
  Gaussian). `[T]+[V]`
- **Explicitly-broken boundary** (`θ=π`) gives a different, higher value — non-vacuously, and it is
  correctly attributed to boundary flux, not to a physical centre floor.

So the tool returns `0` on the free/Gaussian field, the Haar value `0.2051` (sector-independent) on the
flat `β=0` Haar slice, and a positive centre-neutral `≈0.222` at the periodic interacting head —
discriminating, and never manufacturing a centre-sector floor.

---

## 5 · Comparison with P1 / P2 and the twist; scope ledger

**P1 ↔ P3.** P1's `S_W≡0` / `G″=0.2051` is exactly the flat-slice statement; P3 shows it holds in
**every** centre sector (§2.1) and that the periodic measure is `θ`-independent (§2.2), so **P1 is
already centre-neutral** — the centre sector cannot move it.

**P2 ↔ P3.** P2's saturation `c=0.22196`, `Var→4.45` is the **physical periodic** (centre-neutral,
unbroken) floor; it brackets P1's `0.2051` from above and the Haar CLT `0.2688` from below. It is not
a centre-sector floor.

**Twist ↔ P3.** The D5/D3 twist `θ=π` floor `≈0.48` is an **explicitly centre-broken 't Hooft-flux
boundary** (`S_W>0`), *higher* than the physical periodic floor; it is the floor of a *fixed* spatial
holonomy, not of the periodic torus.

| item | statement | label |
|---|---|---|
| constant temporal mode = Polyakov holonomy log (`a=Σξ_{l,3}/√L`) | §1; SM4 §1 | **[T]** |
| `S_W≡0` on the constant-holonomy slice in **every** centre sector | §2.1, witness `[A]` | **[T]**+**[V]** |
| periodic L=1 measure & partition fn **exactly** `θ`-independent (right-Haar) | §2.2, `26.09803` all `θ` | **[T]**+**[V]** |
| flat-chart `G″(a)` `θ`-spread (0.730↔0.013) is a **chart artifact** | §2.2 | **[T]** |
| fixed-boundary sector values move (0.74→0.365 at `θ=0`; ≈0.48 at `θ=π`) but stay >0 | §2.3, witness `[C]` | **[V]** |
| `θ=π` = actionful 't Hooft flux (`B²=−I`, `S_W=2β>0`), not a flat centre sector | §2.3 | **[T]**+**[V]** |
| finite-L centre **unbroken** (`Δ>0`, `⟨tr P⟩=0`), physical ensemble centre-neutral | §2.4 | **[T]**+**[V]** |
| physical periodic floor = Haar-carried P2 `≈0.222` (P1 `0.2051`, Haar CLT `0.2688`) | §3 | **[T]**+**[V]** |
| **no centre-sector floor**; the centre does not change the periodic floor | §0,§2,§3 | **[T]**+**[V]** |
| free Gaussian `G″=κ→0`, orbit avg keeps `κ→0`; β=0 Haar sector-independent | §4 | **[T]** |
| 4D `χ<∞` (uniform-in-L centre-neutral floor) | unchanged | **[O] = H1** |
| **Clay mass gap** | unbuilt | **0%** |

---

## 6 · Witness — `research/v3/w7/p3_polyakov_center.py`

`py -3.12 -I p3_polyakov_center.py [quick|full]` (exit 0, deterministic). `[A]` constant-slice
plaquette check (all `S_W=0`; `θ=π` 't Hooft flux `S_W=2β`); `[B]` full-3D L=1 `θ`-independence via
`sm4_constant_mode.rho_long`; `[C]` exact pure-Cartan transfer-matrix `G″(a₀*)` at `θ=0,π,2π`,
`L=4..8` (`d3_transfer.build_Z/gap_from_Z`); `[D]` transfer-matrix gaps (`bulk_spectrum`); `[E]` acid
test. Quoted output is in §2. No third-party sources; no QED.

---

*Refs.* `docs/paper-v3/main.tex §6.8` (decompactification, the Haar floor `0.2688`, twist-specificity);
`D-WAVE-SYNTHESIS.md`; `D5-UNIVERSALITY.md` (§4 boundary dependence); `D3-TRANSFER-CONSTANT.md`
(gapped transfer matrix `[T]`, `θ=0` saturation at `1/Var`); `SM4-CONSTANT-MODE.md` §1,§3 (Cartan
holonomy class, `θ=π` boundary); `PROOF-GAUSS-INVARIANT-PF.md` (gauge-invariant transfer matrix,
`H_phys`); `research/probe-center/REPORT.md` §4 (exact centre blindness of the neutral/adjoint
sector); `research/lattice-cpb/` (periodic SU(3) MC, `⟨P⟩≈0.17` at `β=2.6` — the unbroken-centre
regime). Standard inputs (right-Haar invariance, 1D no-SSB) cited from the corpus; **no third-party
sources downloaded; no QED.**
