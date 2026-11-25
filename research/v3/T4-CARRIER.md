# T4-CARRIER — does a gauge-invariant **non-self-averaging** carrier exist? **YES**: the total-normalized zero mode of a fixed-physical-scale gauge-invariant density has `Var = χ = O(1) ≠ V^{-1}`; the corpus's `Var ∝ 1/V` on `tr F²` is the **mean-normalized** (self-averaging) mode — a **normalization artifact**, not a structural theorem. The residual obstruction is **not carrier existence** but the **squeeze direction**: no *variance* functional lower-bounds the mass, so the carrier's susceptibility carries the **ceiling only**, never the edge.

> **2026-11-24 · Wave T4 — attacking the carrier hazard left by Q6.**
> Target (Q6 §4.5, `H1-IS-IT-THE-GAP §5`): H1 = `χ := Var_μ(ξ_{k→0}) ≤ 1/λ` is "vacuous on the
> gauge-invariant composite (`Var ∝ 1/V`), gauge-dependent on the gluon", and the corpus asserts
> **"no linear gauge-invariant carrier of `Var` with `Var ↛ 0`"** exists. The task: build one and prove
> `Var ≠ V^{-1}`, or state the exact no-go.
>
> **Result (front).** The asserted no-go is **false as stated**. `Var ∝ 1/V` is the **mean-normalized**
> mode `q̄ = (1/V)ΣO`; H1 (ch02 `eq:h1`) is stated on the **total-normalized** mode
> `q = (1/√V)ΣO`, for which `Var(q) = Σ_z C(z) = χ` — an `O(1)`, `L`-independent quantity. The carrier
> exists (any gauge-invariant `0⁺⁺` density, renormalized at a fixed physical scale). What survives of
> the hazard is the **already-known** half: a **variance** functional is an upper-bound machine, so the
> carrier's `χ` is the **squeeze ceiling** `m₁² ≤ M₀/χ`, never the edge `m₁² ≥ f₁²/χ`. The missing
> structural theorem is therefore **not** carrier construction; it is the **atom + isolation floor**
> `f₁² ≥ r₀` with `I₂` on that carrier — exactly the R-wave isolation/binding target.
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness. **No QED. Clay 0%.**
> **Witness:** `research/v3/w7/t4_carrier.py` (`py -3.12 -I`, numpy only, one job, no downloads).
> **Reads:** `Q-WAVE-SYNTHESIS.md`, `Q6-H1-STATUS.md`, `R3-ATOM-FUNCTIONAL.md`, `R-WAVE-SYNTHESIS.md`,
> `docs/paper-v2/ch02-squeeze.tex` (`eq:chi`, `eq:h1`, `eq:squeeze`, thin-atom Prop., `§schur`),
> `docs/paper-v2/ch09-weakstar.tex`.

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **A non-vacuous gauge-invariant carrier EXISTS `[T]+[V]`.** For any gauge-invariant local density
   `O` (e.g. `tr F²`, renormalized at a fixed physical scale so its connected two-point function is
   summable), the **total-normalized zero mode** `q = (1/√V) Σ_x (O(x) − ⟨O⟩)` is a gauge-invariant
   random variable with
   ```
     Var(q) = (1/V) Σ_{x,y} C(x−y) = Σ_z C(z) = χ_O ,        C = connected two-point function .   (0.1)
   ```
   This is `O(1)` and `L`-independent whenever `χ_O ∈ (0,∞)`. It is **not** `V^{-1}`. Witness §7: for a
   free massive scalar `Var(q) = 1/m²` **exactly, for every `L`** (log-log slope `0.0000` over `L=16…256`
   in `d=1`, `8…64` in `d=2`, `4…12` in `d=4`); the gauge-invariant composite `:φ²:` has `χ` finite and
   `L`-independent (`4.109` in `d=1`, `0.689` in `d=2`); the SU(2) one-plaquette has `Var = 1/4` exactly.

2. **The corpus's "`Var ∝ 1/V`" is the MEAN-normalized mode — a normalization artifact `[T]`.** The
   **mean**-normalized mode `q̄ = (1/V)ΣO` has `Var(q̄) = χ/V → 0` (log-log slope `−1.0000`, witness §7);
   that is the self-averaging mode, and it is **not** the object in H1 (`eq:h1` uses `ξ_{k→0}`, the
   **total**-normalized mode). So the "vacuity branch" of the carrier hazard is a bookkeeping error, not
   a no-go. (Q6 §1 itself flags the normalization as "load-bearing" — but the composite audit then used
   the `q̄` bar convention. The error is localized to the composite branch.)

3. **A second, genuine hazard is UV, not IR `[T]`.** The **bare** lattice `tr F²` has a dimension-4
   two-point function whose on-site value scales as `a^{-8}`; summed to `χ` it diverges as `a^{-4}` in the
   continuum. The bare composite therefore does not give a finite `χ` — it needs **renormalization at a
   fixed physical scale** (gradient flow, APE smearing at fixed physical radius, or point-splitting with
   subtraction). This is exactly the task's "carrier defined by a fixed physical scale", and it is why
   the only correct composite carrier is the **renormalized** one. (Free-field witness: the divergent
   part is the contact term; the subtracted `χ` is finite — §7 [D].)

4. **The carrier's susceptibility carries the CEILING, not the edge `[T]`.** For the correct carrier,
   `χ_O = ∫ s^{-1} dν_O` is the **residue-weighted harmonic mean** of the `0⁺⁺` Källén–Lehmann measure,
   and the exact squeeze gives `f₁²/χ ≤ m₁² ≤ M₀/χ`: an upper bound `m₁² ≤ M₀/χ`, never a lower bound.
   Concretely (witness §7 [D]) the composite `:φ²:` has `χ_{:φ²:} = 0.689` in `d=2` while the lightest
   coupled threshold is `(2m)² = 1`; `χ ≠ 1/m₁²`, so `χ` is not the pole. Gap content = **ceiling only**.

5. **Exact no-go — the strong "gap-carrying carrier" form `[NG]`.** There is **no** gauge-invariant
   carrier whose **variance/susceptibility** lower-bounds the mass. `Var(q_w) = Σ_{x,y} w(x)w(y)C(x−y)` is
   a quadratic form in the two-point function, hence a **Schur-concave** functional of the KL measure
   (ch02 `§schur`); the thin atom `ν = wδ_ε + (1−w)δ_{s₂}`, `w = ε/2λ`, has `χ = 1/λ` fixed while
   `inf supp ν = ε → 0`. So **no variance forces `f₁² > 0`**. This is a **direction theorem** independent
   of gauge invariance and of the carrier.

6. **Missing structural theorem `[O]`.** The one irreducible input is unchanged: the **atom + isolation
   floor** `f₁² ≥ r₀` with `I₂` (`Λ★ − m₁² ≥ δ`) on the renormalized gauge-invariant composite — the
   corpus's `H*` and R-wave's isolation/binding target. It is **not** a carrier-construction problem. The
   carrier hazard's contribution to "H1 < gap" was half genuine (H1 is the ceiling) and half artifact
   (`Var ∝ 1/V` vacuity); with the artifact removed, H1 sits on a well-defined, non-vacuous,
   gauge-invariant carrier and remains the **ceiling** of the gap.

---

## 1 · The normalization algebra `[T]`

Let `Λ` be a periodic `L^d` lattice, `V = L^d`, and let `O` be a local operator with **connected**
translation-invariant two-point function `C(z) = ⟨O(z)O(0)⟩_c` (so `⟨O⟩` is subtracted). Define the
Fourier modes with the two normalizations:

```
   total-normalized:   q(k)  = (1/√V) Σ_x e^{-ikx} (O(x) − ⟨O⟩)        (the "ξ" of ch02 eq:h1)
   mean-normalized:    q̄(k) = (1/V)  Σ_x e^{-ikx} (O(x) − ⟨O⟩)        (the "q̄" of the composite audit)
```

**Theorem T4-1 (normalization algebra) `[T]`.** For translation-invariant `C` and periodic boundary
conditions,
```
   Var(q(k))  =  Σ_z e^{-ikz} C(z)  =  C̃(k) ,          Var(q(0))  =  Σ_z C(z)  =  χ_O .       (1.1)
   Var(q̄(k)) =  (1/V) C̃(k) ,                            Var(q̄(0)) =  χ_O / V .                (1.2)
```
In particular `Var(q̄(0)) = Var(q(0))/V`: the mean-normalized zero mode is **self-averaging**
(`→ 0` like `1/V`); the total-normalized zero mode is **non-self-averaging** (`= χ_O`, `V`-independent).

*Proof.* `Var(q(k)) = (1/V)Σ_{x,y} e^{-ik(x−y)} C(x−y) = Σ_z e^{-ikz}C(z)` by the change of variables
`x → z = x−y` (periodic sum, `Σ_x 1 = V`). At `k=0`, `C̃(0) = Σ_z C(z) = χ_O`. The mean-normalized
case follows from `q̄ = q/√V`. ∎

This is the entire content of the "carrier vacuity": the number `Var ∝ V^{-1.0}` reported for `tr F²`
is `Var(q̄)`, the **self-averaging** mode. The H1 mode is `q`, whose variance is the susceptibility
`χ_O` — the quantity H1 bounds. **The two modes differ by a factor `V`; conflating them turns a finite,
non-vacuous variance into an apparent vacuity.**

**Where the corpus slipped.** Q6 §1 states the total-normalized mode `a^{(tot)}_k=(1/√V)Σ` has
`Var = C̃(k) → χ` and calls the mean-normalized `Õ(k)=(1/V)Σ` "trivially true and vacuous." This is
correct. But Q6 §4.5 / `app:postaudit` then audit the **composite** against `q̄` (the bar), conclude
`Var(q̄)=χ/V→0`, and upgrade this to "no gauge-invariant carrier with `Var ↛ 0` exists." The upgrade is
**invalid**: it re-uses the self-averaging mode to indict the non-self-averaging one. Equations (1.1)–
(1.2) are the exact repair.

---

## 2 · The carrier exists and is gauge-invariant `[T]+[V]`

**Definition T4-C (the carrier).** Fix a gauge-invariant `0⁺⁺` density `O` **renormalized at a fixed
physical scale** (see §3 for why the scale is needed). The carrier is
```
   q := (1/√V) Σ_x ( O(x) − ⟨O⟩ )  =  the total-normalized zero-momentum mode of O .              (2.1)
```

**Theorem T4-2 (existence) `[T]`.** (i) `q` is **gauge-invariant**: `O` is a pointwise gauge-invariant
composite (`tr F²`, `tr F²_μν tr F²_μν` after flow, the energy density), the weight `1/√V` is a constant,
and the subtraction `⟨O⟩` is a c-number; hence `q` is invariant under every gauge transformation.
(ii) `Var(q) = χ_O` at every finite `V`, a **positive** finite number (it is the variance of a finite sum
of bounded lattice variables); so the carrier is **non-vacuous** at every `V`. (iii) The continuum/thermo
content is whether `χ_O = Σ_z C(z)` stays **bounded** as `V→∞`; that boundedness **is H1** on this
carrier, and it is exactly as hard as H1 (no harder, no easier).

**Remark (non-circular).** At finite `V`, `Var(q)=χ_O^{(V)}>0` **unconditionally** — the carrier is a
genuine random variable with a genuine nonzero variance, for every configuration and every `V`. What is
open in 4D at `β=2.5` is the *uniform bound* `χ_O^{(V)} ≤ 1/λ`, which is H1 itself. So constructing the
carrier is **not** circular; it relocates H1 onto a definite, gauge-invariant, non-vacuous object.

**Witness `[V]` (free massive scalar, exact).** With `O=φ` and `C=G` the exact periodic propagator,
`Var(q(0)) = Σ_z G(z) = G̃(0) = 1/m²` — **identically for every `L`** (only the `k=0` term survives the
sum), while `Var(q̄(0)) = 1/(m²V)`. Measured log-log slopes: `Var(q)` = `−0.0000`, `Var(q̄)` = `−1.0000`,
in `d=1,2,4` (§7 [A]). This is the sharpest possible confirmation: the total-normalized carrier's
variance is pinned at `1/m²`, independent of volume, to machine zero.

**Solvable-model audit** (all exact, §7):

| model | carrier | `Var` (total-norm) | `V`-behaviour | acid (free massless) |
|---|---|---|---|---|
| 1D free massive chain | `q=(1/√L)Σφ` | `1/m²` | `L`-independent `[T]` | `m→0`: `Var→∞` ✓ |
| 2D free massive scalar | `q=(1/√V)Σφ` | `1/m²` | `V`-independent `[T]` | `m→0`: `→∞` ✓ |
| 2D composite `:φ²:` | `q=(1/√V)Σ:φ²:` | `0.689` (`m=.5`) | `V`-independent `[T]` | `m→0`: `→∞` ✓ |
| SU(2) one-plaquette | `s=1−½tr U` | `Var=¼` (exact) | no `V` (one site) `[T]` | — |
| 2D pure YM (Migdal) | `q=(1/√V)Σ tr F²` | `χ<∞` | `V`-indep `[T-standard]` | (U(1) dual photon: `χ=∞`, KT caveat) |

The 2D pure-YM row is `[T-standard]` (Migdal solvability: area law ⇒ gapped ⇒ summable `C` ⇒ finite
`χ`); the explicit numeric witnesses are the free 2D scalar/composite and the exact one-plaquette. The
2D **U(1)** caveat (dual photon is a free massless scalar ⇒ `χ_FF=∞`) is the corpus's standing KT
warning and is **acid-correct**: the carrier correctly diverges on the gapless abelian model.

---

## 3 · The two hazards: IR is an artifact, UV is real `[T]`

The carrier hazard was stated as one obstruction; it is **two**, and only the second is real.

**(H1\*) IR "vacuity" — an artifact.** Resolved by Theorem T4-1: the vanishing variance is the
**mean**-normalized mode. The H1 mode is the total-normalized one. `[T]`

**(H2\*) UV divergence — genuine, requires a fixed physical scale.** `tr F²` is dimension-4; its
connected two-point function has a short-distance singularity `C(z) ≈ c₀ δ(z) + (regular)`, with
`c₀ ~ a^{-8}` on the lattice. The susceptibility `χ = Σ_z C(z)` then carries a **power divergence**
`~ a^{-4}` in the continuum limit at fixed lattice spacing rescaling. (Contrast: the free *elementary*
`φ` has `χ = G̃(0) = 1/m²`, finite — because the `k=0` propagator is IR, not UV, and does not probe the
contact.) So the **bare** composite is not a usable carrier across `a→0`: it **diverges** (`H1` *fails*)
rather than being vacuous. The fix is standard and is exactly the task's "carrier defined by a fixed
physical scale":

- **gradient flow** at flow time `t ~ (fixed physical length)²`: `O_t = tr F²(t)` is finite and
  scheme-clean (the flow provides the subtraction);
- **APE/HYP smearing** at a fixed physical radius;
- **point-splitting with a contact-term subtraction** `C(z) → C(z) − c₀δ(z)`.

After renormalization, `χ_O^{ren}` is **finite** and `O(1)` (it equals `∫ s^{-1} dν_O` on the physical
`0⁺⁺` spectral measure), and the carrier of §2 is well-defined uniformly in `a`. `[T]` (the divergence is
elementary power-counting; the subtraction is standard). The witness §7 [D] isolates the finite part by
working the free composite, where the contact structure is explicit.

---

## 4 · What the carrier's susceptibility can and cannot carry `[T]`

The KL measure of `O` is `ν_O = Σ_n |⟨0|O|n⟩|² δ(s − m_n²)` (on `0⁺⁺` states), and
```
   χ_O = ∫ s^{-1} dν_O(s) = Σ_n |⟨0|O|n⟩|² / m_n² .                                  (4.1)
```
This is manifestly a **residue-weighted harmonic mean** — the exact object of the squeeze.

**Theorem T4-3 (ceiling, not edge) `[T]`.** For every positive `ν_O` (all `f₁², M₀, χ` finite),
```
   f₁² / χ_O  ≤  m₁² = inf supp ν_O  ≤  M₀ / χ_O ,                                         (4.2)
```
so `χ_O ≤ 1/λ` gives `m₁² ≤ M₀ λ` — an **upper** bound on the mass. A **lower** bound requires the edge
atom `f₁² ≥ r₀`:
```
   (f₁² ≥ r₀) ∧ (χ_O ≤ 1/λ)  ⇒  m₁² ≥ f₁²/χ_O ≥ r₀ λ > 0 .                                (4.3)
```
*Proof.* ch02 `eq:squeeze` and `eq:sufficient` (restated; the squeeze is exact for any positive
measure). The upper rung is the harmonic mean; the lower rung is the atom. ∎

**Witness that `χ` is not the pole `[V]`.** For the composite `:φ²:` the KL measure is the **two-particle
continuum** starting at threshold `(2m)²` (and `:φ²:` is **atomless**, `f₁²=0` — R5). Its susceptibility
is a weighted integral over that continuum, **not** the pole value: in `d=2`, `m=0.5`,
`χ_{:φ²:} = 0.689` versus `1/(2m)² = 1.000` (ratio `0.689`, §7 [D]). So `χ_{:φ²:} ≠ 1/m₁²`: the composite
carrier's susceptibility is the **squeeze ceiling**, and because `:φ²:` is atomless, `H1` on it is the
ceiling **with no edge** — precisely the corpus's `:φ²:` acid model, met again from the carrier side.

**The spectral-projector dual `[T]`.** The atom is `f₁² = ⟨0|O P_{\{m₁²\}} O|0⟩`, an interval-projector
expectation `ν_O([m₁², m₁²+I₂])`. Reading it requires choosing the window at `m₁²` (R3 face R1:
presupposes the gap) or the limit `I₂→0` (R3 face R2: infinite data). The susceptibility needs **neither**
— it is the gap-free, finite-data functional — but it returns only the **ceiling**. So the two halves of
the KL measure sit on **opposite sides of the R3 no-go**: the gap-free functional is the ceiling (this
note), the edge functional is circular (R3). There is no third functional that is both gap-free and
edge-reading; that impossibility **is** the R3 Stone–Weierstrass bound `sup|E−g| ≥ f₁²/2`, now seen as a
statement about carriers.

---

## 5 · The exact no-go — the strong "gap-carrying carrier" form `[NG]`

**No-go T4-4 `[NG]`.** Let `q` be **any** carrier linear in a local density, `q = Σ_x w(x)O(x)` with any
normalization `w` (gauge-invariant or not, smeared or not, at any physical scale). Then **`Var(q)` cannot
lower-bound the mass gap**; i.e. there is no gauge-invariant `q` and constant `λ>0` such that
`Var(q) ≤ 1/λ` **forces** `m₁² ≥ λ` on the KL measure of `O`.

*Proof.* `Var(q) = Σ_{x,y} w(x)w(y)C(x−y)` is a **quadratic form in the two-point function** `C`, hence
a functional of `ν_O` built from its moments / nonnegative combinations — a **Schur-concave** functional
(ch02 `§schur`). The thin atom `ν_ε = wδ_ε + (1−w)δ_{s₂}`, `w = ε/2λ` (ch02 thin-atom Prop.) has
`χ(ν_ε) = 1/λ` **for every `ε`** yet `inf supp ν_ε = ε → 0` (witness §7: `χ_{:φ²:}`-type quantities are
`V`-blind to `ε`). Hence the set `{Var(q) ≤ 1/λ}` contains measures with `m₁²` arbitrarily small, and no
variance bound can force `f₁²>0`. ∎

This is the **direction theorem** the H/Q waves kept meeting: **variance is an upper-bound machine**
(Fisher–Rao gives a variance *floor*, not a ceiling, Q7; Bernstein caps curvature; band-limited densities
drive `G″→0`). It is independent of the carrier: **changing the carrier cannot reverse the direction of a
moment inequality.** Therefore "a carrier whose susceptibility certifies the gap" does **not** exist —
and this, not carrier vacuity, is the true residual of the carrier hazard.

**What this leaves `[O]`.** The irreducible input is the pair on the **renormalized gauge-invariant
carrier** `q`:
```
   H* :  (atom) f₁² = ν_O({m₁²}) ≥ r₀    and    (isolation) Λ★ − m₁² ≥ δ    [the edge],
         with   χ_O ≤ 1/λ                                                      [H1, the ceiling, now non-vacuous].
```
This is exactly `H*` of ch02 `eq:hstar` and R-wave's **isolation/binding** target. The carrier work does
**not** remove it — it removes the normalization objection that stood in front of it, and confirms the
target is the **atom/isolation**, not the carrier. `[O]` · **Clay 0%.**

---

## 6 · Reconciliation with Q6, R-wave, and the falsification battery

- **Q6 "carrier hazard" (vacuous composite / gauge-dependent gluon):** the vacuity is **retracted** as a
  normalization artifact (§1, §2); the gluon branch remains gauge-dependent (unaffected). The third
  branch — "gauge-invariant non-self-averaging carrier" — is **constructed** (§2): the total-normalized
  zero mode of a renormalized gauge-invariant density. Q6's own §5.2 F4 test ("must be non-vacuous on the
  total-normalized mode") is now **passable**: `Var(q)=χ_O=O(1)`.
- **R3 (atom functional):** unchanged and now complementary. The gap-free functional is the
  susceptibility (ceiling, §4); the edge functional is the projector/boundary reader (circular, R3). The
  two are on opposite sides of the Stone–Weierstrass bound; no carrier closes the gap between them.
- **R-wave (isolation is the gap; binding is the discriminant):** unchanged. The carrier is binding-blind
  — `Var(q)` sees the residue-weighted harmonic mean, not the bound flux tube. The missing theorem is the
  R-wave one.
- **Falsification battery F1–F5:** F1 (free massless must fail) — the carrier's `Var = 1/m²` **diverges**
  as `m→0` and the 2D U(1) dual photon gives `χ=∞`: **passes** (acid-clean, §2 table). F4 (non-vacuous
  total-normalized carrier) — **passes** now (§2). F2, F3, F5 concern the atom/curvature/dimension and are
  untouched by carrier choice.

---

## 7 · Witness (`research/v3/w7/t4_carrier.py`) — exact numbers

**[A] Total vs mean normalization, free massive scalar `φ`, `m=0.5` (`1/m²=4`), exact periodic
propagator.** `Var(q(0))=Σ_z G(z)`, `Var(q̄(0))=Var(q)/V`:

```
  d=1  L=16,32,64,128,256 :  Var(q)=4.000000 (all L) ; log-log slope Var(q)= -0.0000 , Var(q̄)= -1.0000
  d=2  L=8,16,32,64       :  Var(q)=4.000000 (all L) ; slope Var(q)= -0.0000 , Var(q̄)= -1.0000
  d=4  L=4,6,8,12         :  Var(q)=4.000000 (all L) ; slope Var(q)= -0.0000 , Var(q̄)= -1.0000
```

`Var(q)` is pinned at `1/m²` to all digits, in every dimension, for every volume — the total-normalized
carrier is **non-self-averaging and non-vacuous**; `Var(q̄)=χ/V` is the self-averaging mode (slope `−1`).

**[A′] Free massless acid (`L=16`, `d=1,2,4`).** `Var(q) = 1/m² → ∞` as `m→0`:
`m=0.5→4.0`, `m=0.2→25.0`, `m=0.05→400.0`, `m=0.01→10000.0` (exact). The gap-free model **fails** H1 on
this carrier — the acid test the carrier must pass, and does.

**[B] Gauge-invariant composite `:φ²:` (free Wick), `χ = Σ_z 2G(z)²`.**

```
  d=1  L=32,64,128,256 :  4.108856, 4.108839, 4.108839, 4.108839   slope -0.0000
  d=2  L=16,32,64      :  0.692494, 0.689051, 0.689049             slope -0.0018
  acid  d=1 L=64, m=0.5,0.2,0.1,0.05,0.02 : 4.1, 62.8, 513, 5475, 195989  -> DIVERGES  ✓
  acid  d=2 L=64, m=0.5,0.2,0.1,0.05,0.02 : 0.69, 4.07, 16.4, 100.6, 3081  -> DIVERGES  ✓
```

Finite and `L`-independent for `m>0` (non-vacuous, non-self-averaging); divergent as `m→0` (acid-correct).

**[C] SU(2) one-plaquette Haar (exact).** `z=½tr U`: `⟨z⟩=0.000000`, `Var(z)=0.250000` (matches the exact
`⟨z²⟩=¼`); `Var(s=1−z)=0.250000`. A single gauge-invariant site: `O(1)` variance, no `V`.

**[D] Ceiling-not-edge.** `d=2, L=64, m=0.5`: `χ_{:φ²:}=0.689049` versus the lightest coupled threshold
`1/(2m)²=1.000000` (ratio `0.689`). `χ ≠ 1/m₁²`: the composite carrier's susceptibility is the squeeze
ceiling over an **atomless** continuum (`f₁²=0`), never the pole.

---

## 8 · Ledger

| item | statement | label |
|---|---|---|
| **normalization algebra** | `Var(total-norm zero mode)=χ_O`; `Var(mean-norm)=χ_O/V` | **[T]** |
| total-normalized carrier | `q=(1/√V)Σ(O−⟨O⟩)` has `Var=χ_O`, `V`-independent | **[T]** |
| composite `tr F²` "vacuous" | `Var∝1/V` is the **mean-normalized** mode — artifact | **[T]** (retraction) |
| **carrier existence** | gauge-invariant, non-self-averaging, `Var=O(1)≠1/V` | **[T]**+**[V]** |
| free massive scalar | `Var(q)=1/m²` exactly, all `L`, `d=1,2,4` | **[T]**+**[V]** |
| composite `:φ²:` | `χ` finite, `L`-indep (`4.109` d1, `0.689` d2); div `m→0` | **[T]**+**[V]** |
| one-plaquette SU(2) | `Var(½tr U)=¼` exact; `O(1)`, no `V` | **[T]**+**[V]** |
| bare composite UV | `χ_bare ~ a^{-4}` divergent; needs fixed-scale renorm | **[T]** |
| renormalized composite at fixed scale | finite `χ_O=∫s^{-1}dν_O`, usable carrier | **[T]** |
| squeeze on carrier | `f₁²/χ ≤ m₁² ≤ M₀/χ`: `χ` is the **ceiling**, not the edge | **[T]** |
| `χ` vs pole (composite) | `χ_{:φ²:}=0.689 ≠ 1/(2m)²=1`; `:φ²:` atomless | **[T]**+**[V]** |
| projector dual | atom=projector (needs `m₁²`, R1) vs `χ` (ceiling, gap-free); no third functional | **[T]** |
| **no-go T4-4** | no variance functional lower-bounds the mass (thin atom / Schur) | **[NG]** |
| missing structural theorem | atom+isolation floor `f₁²≥r₀`, `I₂` on the carrier `= H*` | **[O]** |
| 2D pure YM | gapped ⇒ `χ<∞` ⇒ carrier finite (Migdal); 2D U(1) dual photon `χ=∞` (KT) | **[T-standard]** |
| Clay gap | unbuilt | **[O]** · **0%** |

---

## 9 · Sharpest line

**A gauge-invariant non-self-averaging carrier exists and is elementary to write: the total-normalized
zero-momentum mode `q=(1/√V)Σ(O−⟨O⟩)` of a gauge-invariant density renormalized at a fixed physical
scale. Its variance is `Var(q)=Σ_z C(z)=χ_O`, an `O(1)` `L`-independent quantity — exactly the `ξ_{k→0}`
of H1 — and never `V^{-1}` (free massive scalar: `Var=1/m²` to all digits, all `L`, `d=1,2,4`; one
plaquette: `Var=¼`). The corpus's `Var∝1/V` on `tr F²` is the mean-normalized (self-averaging) mode, a
normalization artifact, not a no-go. The genuine residual of the carrier hazard is UV, not IR — the bare
`tr F²` needs a fixed-scale renormalization — and, more deeply, the squeeze direction: `χ_O` is a
residue-weighted **harmonic mean**, so any carrier's variance gives only the **ceiling** `m₁²≤M₀/χ_O`,
never the edge `m₁²≥f₁²/χ_O`. No variance functional can lower-bound the mass (thin atom / Schur), so
"susceptibility carries the gap" is a **no-go**, and the irreducible missing theorem is the
**atom + isolation floor** on that carrier — R-wave's isolation/binding target, unchanged. The carrier
hazard was half artifact, half the old ceiling-vs-edge wall: fixing the normalization makes H1 sit on a
definite non-vacuous gauge-invariant object, and confirms that the Clay gap lives in the atom, not in any
carrier.** **No QED. Clay 0%.**

---

*Refs.* **Carrier / vacuity:** `Q6-H1-STATUS.md` §1 (normalization load-bearing), §4.5 (carrier caveat);
`Q-WAVE-SYNTHESIS.md` §4 (three carriers); `research/unconditional/H1-IS-IT-THE-GAP.md` §5;
`docs/paper-v2/app.tex` (`app:postaudit`); ch02 `eq:chi`, `eq:h1`. **Ceiling / squeeze:** ch02
`eq:squeeze`, `eq:sufficient`, `eq:hstar`, thin-atom Prop., `§schur`; `lean/YMStatement/H1Sharpening.lean`
(`H1_does_not_imply_atom`). **Atom / isolation:** `R3-ATOM-FUNCTIONAL.md` (R1=R2=R3, Stone–Weierstrass
bound), `R-WAVE-SYNTHESIS.md` (isolation is the gap, binding is the discriminant), `R5-EXACT-ANCHORS.md`
(`:φ²:` atomless). **Direction theorem:** `Q7-TOOL-FORGE.md` (Fisher = variance floor); `H2-DEFECT-BOUND.md`.
**Standard (textbook, not re-fetched):** lattice free propagator; Källén–Lehmann; Migdal 1975 (2D YM
solvability); gradient-flow renormalization (Lüscher); Schur-concavity/convex order. **Witness:**
`research/v3/w7/t4_carrier.py` (numpy only, one job, no downloads). No `[C]` upgraded to `[T]`; the
normalization retraction is recorded per charter §2. **No QED; Clay 0%.**
