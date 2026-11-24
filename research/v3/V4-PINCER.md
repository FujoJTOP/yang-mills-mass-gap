# V4-PINCER — assembling the pincer on (★): what the strong-coupling jaw and the UV jaw actually cover of the tower, and what the head costs

> 2026-11-24 · Coverage / verdict note for the single live inequality
> ```
> (★)  Var_{μ(·|U′)}(⟨Pᵀ∇S_j, v⟩)  ≤  E_{μ(·|U′)}⟨PᵀH_jP v, v⟩ + W_j(v) − δ_W ,   δ_W>0,
> ```
> plus cross-shell decay `(★dec)`, uniform in shell `j`, boundary `U′`, cutoff `a`.
> Reads: `S6-DOBRUSHIN.md`, `S1-CONDITIONAL-POINCARE.md`, `S5-N-UNIFORM.md`, `C1-SZZ-NSCALING.md`,
> `C2-BALABAN-UNIFORM.md`, `STAR-WAVE-SYNTHESIS.md`, `G1-BLOCK-IDENTITY.md`, `G2`, `G3-PROJECTION-TERM.md`,
> `G4-GENERAL-BOUND.md`, `research/unconditional/G1-L3.md`, `G1-L3-TOWER-POINCARE.md`,
> `G1-L3-UV-FORCE-BOUND.md`, `PROOF-COARSE-POINCARE.md`, `PROOF-JACOBIAN-SAVES.md`.
> Witness: `research/v3/w7/v4_pincer.py` (shell-index arithmetic, output quoted §6).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step · **[V]** numerical witness · **[R]** reported/unrefereed.
>
> **Missing inputs, stated up front.** `V1-UV-NORMALIZATION.md`, `V2-TRUE-BOUND.md`,
> `V3-COUNTEREXAMPLE-SCOPE.md` **do not exist in the repo** (checked this session). Their intended
> content is carried by, respectively: the renormalised UV statement `L3 = sup_a Z_a⁻¹‖H_a‖_D < ∞`
> (`G1-L3.md` §1, `C2-BALABAN-UNIFORM.md` §1); the small-field bound `Corollary S1′`
> (`S1-CONDITIONAL-POINCARE.md` §3) together with Route-B cone-emptiness
> (`G1-L3-UV-FORCE-BOUND.md` §3); and the counterexample scope `G4-GENERAL-BOUND.md` §3–5 (which
> *supersedes* the S2 reduction — see §3.2). This note assembles the pincer from the files that exist.

---

## 0 · Verdict (no QED)

**The pincer is real but *unbalanced and conditional*, and it is not the theorem the task hoped for.**

1. **The strong-coupling jaw is a genuine `[T]` — and it grips zero tower shells.** S6's
   Dobrushin–Shlosman theorem gives `(★)∧(★dec)` for `β_j < β_Dobr = 1/72`; SZZ (via C1/Lemma SC)
   gives an `N`-uniform strong-coupling gap for `β < 1/48` (standard-Wilson `β_W ≈ 0.19`). But the
   tower's every shell sits at `β_j ∈ [2.5, ∞)`. The windows `(0, 1/72]`, `(0, 0.19]` and `[2.5, ∞)`
   are **disjoint**; the covered shell set is **empty** (`§2`, `[T]` arithmetic).

2. **The UV jaw is `[C]`, not `[T]`, and it is the only jaw that grips the tower.** Where the cone
   empties (`|F_j| = O(g_j)` and the block Laplacian floor beats the Savvidy negative term), the fibre
   is log-concave and Brascamp–Lieb gives the equipartition `Var ≤ E[H]`, reducing `(★)` to
   `W_j ≥ δ_W` (equivalently, on G4's corrected reading, `E[Σ] ≥ δ_W` with the Haar floor
   `Φ″ ⪰ N/12`). This covers **every shell `j ≥ j_UV` with `j_UV = O(1)`**, i.e. **all but `O(1)`
   shells of an arbitrarily deep tower**. But it is **[C]** on the missing gauge-invariant,
   `a`-uniform moment bootstrap `E|F_j| ≤ C g_j` (**[O]**, circular with the flow per
   `G1-L3-UV-FORCE-BOUND.md` §2), and only on the **small-field / excised** domain `(★_sf)` — not for
   all `U′` (`§3`).

3. **The head is not "one shell" — it is the `O(1)`-coupling IR bottom, and it is exactly the mass
   gap.** The uncovered set is a **finite, tower-depth-independent** band of `O(1)` shells at the
   bottom (physically `β ≈ 2.5`, `j ≲ 3` for `N=3` / `≲ 7` for `N=2`). There, `α(2.5) ≈ 180` (Dobrushin
   dead) **and** `PᵀH_jP` is Savvidy-indefinite (BL dead); the open object is
   `G″(v;a₀) = E[Σ] − D ≥ δ_W` at the RG-relevant coarse field, i.e. the **coarse spectral gap = H1**
   (`§4`, `[O]`). Moreover the strong "uniform in all `U′`" reading is **refuted** by G4's exact
   one-plaquette counterexample (`G″ < 0` at `β > 4.09`, maximal holonomy, double-well barrier) — the
   honest target is `(★_sf)`, not `(★)` on all of `SU(N)^{∂B}`.

4. **There is no "whole intermediate regime" left inside the tower.** The tower lives at `β ≥ 2.5`;
   the intermediate `β`-window `(0.19, 2.5)` is **not inhabited** by the physical trajectory (only by
   the one-loop extrapolation of the lowest, already-head shells). Quantitatively (`§6`), the
   uncovered band is `j = 1 … j_UV−1` with `j_UV ∈ {1,2,3}` — `1` shell at `J = 6,10,14` on the literal
   one-loop index — **independent of `J`**. So the union covers all but `O(1)` shells, not all but an
   interval.

5. **Therefore the honest label of the whole pincer is `[C]`/`[O]`, not `[T]`:**
   ```
   (★) on the tower  ⟸  [UV jaw: E|F_j| ≤ C g_j  +  (★dec),  all j ≥ j_UV, on (★_sf)]   [O]/[C]
                         ∧  [head: G″(v;a₀) ≥ δ_W  at j ≲ j_head,  i.e. H1]              [O]
   ```
   The strong jaw contributes **nothing on-tower**; the UV jaw is **conditional**; the head is **the
   Clay core**. Proving "all shells except the head" would be genuine progress on **fixed-volume
   existence (`G1`/OS0)**, and **none** on the mass gap: the head *is* the gap (`G1-L3.md` §5,
   `C2-BALABAN-UNIFORM.md` §4). Clay = 0 %. No QED.

---

## 1 · The object, the two `β`-conventions, and the shell index

The tower of conditionals is `ν_j = P_{j*}\,μ_{j+1}` over shells `j = 1,…,J`, `J = ⌊log_L(Λ/a)⌋`
(`C2-BALABAN-UNIFORM.md` §2.1; `G1-L3-TOWER-POINCARE.md` §1). Two bookkeepings coexist and must be
kept distinct:

| convention | coupling | head | deep UV | source |
|---|---|---|---|---|
| **physical / standard Wilson** | `β_j = 2N/g_j²`, `β_L ≈ 2.5` at the lattice bottom | `β ≈ 2.5`, `j = 1` (IR bottom) | `β_j → ∞`, `j → J` | S6 §1; G1-L3-TOWER §3 |
| **one-loop asymptotic** | `g_j² = 1/(2β₀ j)`, `β_j = 4Nβ₀ j`, `β₀ = 11N/(48π²)` | `β_j = 2.5` at `j_head` | `β_j ∝ j`, `j → J` | C2 §2.1; S5 §1 |

In the one-loop index `β_1 = 4Nβ₀ = 0.37` (`N=2`) / `0.84` (`N=3`); `β_j = 2.5` at
`j_head = 2.5/(4Nβ₀) = 6.7` (`N=2`) / `3.0` (`N=3`). The **physical** bottom sits at `β_L = 2.5` at `j=1`;
the one-loop expression is the UV asymptote, offset by a few shells from the physical bottom. **Both
bookkeepings agree on every conclusion below**; the shell indices differ only by this offset, which I
flag wherever it is used.

**Two layers of `(★)`.** `S6`/`S1`/`S5` write the literal form (L):
`G″_lit = E[H] − Var(G) + W_geo + X`, `W_j = W_geo + X`, `X = −2Cov(G,Φ′)`. `G4-GENERAL-BOUND.md`
shows that the object `(★)` is *about* — the coarse-mode marginal curvature — is (C):
`G″(v;a₀) = −∂²_a log ρ_v(a)|_{a₀} = E[Σ|a₀] − D`, `Σ = vᵀHess V v − ⟨K⁻¹c,c⟩`, and that on the
**unconditional** measure `G″_lit ≡ 0` (Lemma G2). The pincer's conclusions are robust to this
distinction (§3.2); where a file's reduction depends on (L) I say so.

---

## 2 · The strong-coupling jaw: `[T]`, off-tower, covers **∅**

### 2.1 The theorems, stated

- **S6 / Dobrushin [T].** On the shell-`j` block fibre, the Dobrushin coefficient is `α(β_j) = κ β_j`
  with `κ = O(d) ≈ 72` (`S6` §2, eq. (2.2)). `α < 1` at `β_j < β_Dobr = 1/72`, and then
  Dobrushin–Shlosman complete analyticity gives **(★) and (★dec)** verbatim, uniformly in `U′, a`
  (`S6` eq. (2.3)).
- **SZZ via C1 / Lemma SC [T], N-uniform.** SZZ's Bakry–Émery constant is `K_S = N(½ − 24|β|)`,
  positive for `|β| < 1/48` with **`N`-free radius** (`C1-SZZ-NSCALING.md` §§1–3); the per-colour
  constant `K_S/N` is `N`-independent. Pushed through the 1-Lipschitz block map `P` (`C_blk = O(1)`,
  Lemma SC), this gives the coarse Poincaré/LSI at `β < 1/48` (`PROOF-COARSE-POINCARE.md` (c)). Standard
  Wilson `β` for the same small coupling is `β_W ≈ 0.19` (`S6` §2 table, `ATOM-strong-coupling-extend`).

### 2.2 The shell range: empty

Every tower shell has `β_j ≥ β_1`:

| `N` | `β_1` (one-loop) | `< 1/72`? | `< 0.19`? | first shell `> 1/72` | first shell `> 0.19` |
|---|---|---|---|---|---|
| 2 | 0.3715 | no | no | `j=1` | `j=1` |
| 3 | 0.8359 | no | no | `j=1` | `j=1` |

So the number of integer shells `j ≥ 1` in the strong-coupling window is **0** for both `N` — the
thresholds are crossed at `j = 1/72·(4Nβ₀)⁻¹ < 1`. In the physical convention (`β_j ∈ [2.5,∞)`) it
is `∅` by inspection. **The Dobrushin/SZZ jaws prove `(★)∧(★dec)` for zero shells of the actual
tower** (`S6` §0 item 1, §4; arithmetic verified by `v4_pincer.py`, §6). This is a `[T]` statement
about an *off-tower* window; its relevance to the tower is `[O]` (none).

### 2.3 Why it cannot be pushed to the tower (`[NG]` for the method class)

`α(2.5) ≈ 180 ≫ 1`; the disagreement-percolation cluster is supercritical (`S6` §3); and, more deeply,
at `β_j ≥ O(1)` the fibre is not a small perturbation of a high-temperature measure but a sharp, curved
one whose correct linearisation is log-concavity (UV end) or Savvidy-indefinite (head) — neither of
which the influence/contraction toolbox reaches (`S6` §4). **`[NG]` for the local-influence /
Dobrushin / disagreement-percolation class at `β ≈ 2.5`.** No `1/48`-style threshold moves this: the
`1/48`–`2.5` distance is `≈ 13×` in Hooft `λ` (not the naive `120×`), and it is a *gap*, not a
constant to optimise (`S6` §2).

---

## 3 · The UV jaw: `[C]` on `[O]`, covers `j ≥ j_UV` on the small-field domain

*(This is the "V1 renormalised statement" ⊕ "V2 small-field bound" of the task, mapped onto the files
that exist.)*

### 3.1 The renormalised/UV-stable statement (`G1-L3`, `C2`)

`L3 = sup_a Z_a⁻¹‖H_a‖_D < ∞` on the smeared gauge-invariant core `D` (`G1-L3.md` §1;
`G1-UNIFORM-ZA` §1). Balaban's theorem is `a`-uniform **at fixed volume, at the level of the effective
action/densities, not the force Hessian** (`C2-BALABAN-UNIFORM.md` §1.1): Hessian control ⇒ action
control, not conversely, and the `R`-operation bounds are absolute-value. So the published
constructive-RG chain stops *at* L3; it does not supply the one-sided variance bound (`[T]` as
published; `[O]` for L3). **This is why the UV jaw is conditional, not a citation.**

### 3.2 The small-field / cone-empty bound (the mechanism that *would* close the UV shells)

Deep in the UV the field is `O(g_j) → 0`; the Savvidy cone empties
(`√(2g_j|F_j|) ∼ g_j ≪ Λ_{j+1}`, `G1-L3-UV-FORCE-BOUND.md` §3, Route B), so `H_j ⪰ c Λ_j² > 0` and the
fibre is log-concave. On the small-field boundary set `𝒮_t` the conditional modulus is explicit,
`ρ(U′) ≥ min{c_□ + N/12 − C₁βF_∂(U′), N/12 + R₀}` (**Lemma S1**, `[T]`; **Corollary S1′**, `[C]` on
the plaquette constant `C₁`), giving uniform `LSI/BL` and `D ≤ 0` (`S1` §3). Then Brascamp–Lieb gives the
exact equipartition
```
Var_{μ}(G_v) ≤ E[H_vv]      (equality iff Gaussian),     (S6 eq. (5.2); S1 §2(b))
```
and `(★)` collapses to the **Haar budget** `W_j ≥ δ_W`, which is **Lemma B** (`Φ″ ⪰ N/12`, hence
`W_geo ≥ 0` scale-free; `PROOF-JACOBIAN-SAVES.md` Lemma B; `S5` §3). On G4's corrected (C) reading the
same collapse is `G″ ≥ E[Σ] ≥ N/12 − Schur-loss` on the log-concave fibre (`G4` §3 "What *is*
provable"), i.e. the identical Haar-floor conclusion. **The two readings agree on the UV shells.**

**What the UV jaw covers, precisely.** All shells with `j ≥ j_UV`, where `j_UV` is the first shell at
which the Laplacian floor beats the Savvidy negative term (`§6`: `j_UV ∈ {1,2,3}` across the `O(1)`
constant scan; canonical `κ = c_L = 1` gives `j_UV = 2`). That is **`J − j_UV + 1` shells — all but
`O(1)` of an arbitrarily deep tower.**

**What it does not cover — the conditions.** (i) the moment bootstrap `E|F_j| ≤ C g_j` (**`[O]`**,
`G1-L3-UV-FORCE-BOUND.md` §6 item (1); circular with the flow); (ii) the equipartition defect /
`(★dec)` bound on `Var − E[H]` (**`[O]`**, §6 item (2); the `a⁻⁴` renormalon makes the bridge
*asymptotic*, not uniform); (iii) it holds on `(★_sf)` (small-field, Lemma-B-excised), **not** for all
`U′` — large-field boundaries make `K` indefinite and `G″` can be negative there (`S1` §4.1 `[NG]`;
`G4` §3(iv) counterexample `G″ < 0` at `β>4.09`, maximal holonomy, double-well barrier).

### 3.3 `N` and `a` are not the obstruction (`S5`, `[T]`/`[C]`)

`δ_W = c_w·N − O(1)` with `c_w ≈ 1/12–0.13` `N`-independent (adjoint dimension `N²−1` is multiplicity,
not magnitude: SZZ Lemma 4.1), and `δ_W = Θ(a⁰)` (`W_geo` coupling-free, bare Hessian growth `× g_j²`
cancelled). Conditioning adds no factor to the constant; the only exponential is the measure of the
large-field set, `e^{−cNj}` (`S5` §§3–5). **`[T]` mechanism / `[C]` synthesis — but this is coupling-
and `a`-bookkeeping, not coupling-window coverage.**

### 3.4 The cross-shell leg (`(★dec)`, `[O]`)

Even granting per-shell `(★)`, the tower min needs the off-diagonal covariance bound
`(★dec): |Cov(s_B,s_{B′})| ≤ C e^{−κ dist/a}` (`G1-L3-TOWER-POINCARE.md` §5; `S4-DECIMATION.md`).
`(★)` and `(★dec)` are the diagonal and off-diagonal of one covariance matrix and are **independent**
(`S4` `[T]`: `(★) ⇏ (★dec)` and vice versa). The bridge `(★) ⇒ (★dec)` needs uniform global fibre LSI
(`H2-global`, `C(β)<1`), which **fails at the head** (`C ≈ 15`; `S4` `[C]`). So the tower conjunction
`G″ ⪰ ½δ_W` (Gershgorin, `S4`) is `[O]` on the `(★dec)` half.

---

## 4 · The head: `β ≈ 2.5`, typical field — precisely what is open

**Configuration.** The IR bottom, `β ≈ 2.5`, `j ≲ j_head` (`≈3` for `N=3`, `≈7` for `N=2`), typical
field `g = O(1)` (not `O(g_j)`), RG-relevant coarse field `a₀` = the coarse **minimum** (not the
barrier).

**Why both jaws fail there — two independent, named reasons.**
- *Dobrushin:* `α(2.5) ≈ 180 ≫ 1` (`S6` §2) — local influence is supercritical.
- *Brascamp–Lieb:* `PᵀH_jP` has a Savvidy negative mode (min `≈ −0.3`, `MARGINAL`; `W6-BE`:
  `Hess S_W` min `−1.12β`, 21.7 % negative directions), so the fibre is **not** log-concave,
  `(Hess A)⁻¹ ⪯ H_j⁻¹` is false, and the equipartition `Var ≤ E[H]` fails.

**What is still open, stated once (this is exactly H1).** With the G4 corrected object (C),
```
G″(v;a₀) = E_{μ}[Σ|a₀] − D  ≥  δ_W ,     Σ = vᵀ Hess V v − ⟨K⁻¹c,c⟩ ,
```
i.e. the **coarse effective action `S_eff,j` is `δ_W`-strongly convex along every coarse mode `v` at
the RG-relevant coarse field**. By `S3-UV-FESHBACH`/`G4` Lemma G3, this is the conjunction
```
[fibre log-concavity enough that D ≤ E[Σ] − δ_W]   ∧   [E[Σ] ≥ δ_W].
        (small-field input, §3.2 [C])                  (coarse spectral gap = H1, [O])
```
The second factor is the **mass gap in curvature language** (`G″ ≈ 1/Var(a)`, i.e.
`Var(a) ≤ 1/δ ⟺ m_eff²(k→0) ≥ δ`; `PROOF-JACOBIAN-SAVES.md` §2, `PROOF-COARSE-POINCARE.md`).

**Measured, not proved.** `G″ = +0.27 … +0.51` at `β=2.5` (`W7`), `D/S ≈ 2 %` (`PROVE-Gpp`), with
`W_full ≈ +0.25–0.54` carrying the sign (`PROOF-JACOBIAN-SAVES` §3). This is **[V]**, not a theorem.

**The strong all-`U′` form is false, and the target must be corrected.** `S1` Prop U′-NG (`[T-standard]`+
`[V]`): any chromomagnetic boundary makes `K` indefinite on the soft mode (`λ_min: 0 → −1.43 → … → −8`).
`G4` §3(iv) (`[T]` in the one-plaquette model, `[V]` numbers): at maximal holonomy, longitudinal mode,
double-well barrier, `G″ = 0.206 − 0.0123β²` crosses `0` at `β* ≈ 4.09` and is **negative** for
`β > β*` (`β=5: −0.10`, `β=10: −1.02`). So `(★)` **uniform in all `U′` is disproved**; the tenable form
at the head is `(★_sf)` (uniform on `𝒮_t`, S1 Corollary) or `(★_avg)`, with `δ_W` allowed to vanish on
the excised large-field set (`STAR-WAVE-SYNTHESIS` §2). **This is a correction, not a defect: it says
the head is the gap *and* the gap is a strong-field phenomenon.**

**One further named head term (covariant projection).** A genuinely covariant (parallel-transport)
block spin adds `R_proj` to `G″`: exact identity `G″_P = 2W_geo + R_proj` (`G3` Theorem, `[T]`+`[V]`),
with `R_proj` sign-indefinite and `O(W_geo)` at the head (`|R_proj|/(2W_geo)` up to `0.76` in the toy),
while `O(g_j)`-small in the UV tower (`G3` §6, `[C]`). A genuine head bound on `|R_proj|` joins
`(★dec)`/H1 (`G3` §6, `[O]`).

---

## 5 · Coverage table

Tower shells `j = 1 (head) … J (deep UV)`; "phys" = physical `β∈[2.5,∞)` convention; "1-loop" = the
`β_j = 4Nβ₀j` index. `[T]` = theorem, `[C]` = conditional, `[O]` = open, `[NG]` = no-go for a class.

| shell range | theorem / tool | status | what it does **not** cover |
|---|---|---|---|
| **all `j`, `β_j < 1/72`** | S6 Dobrushin–Shlosman: `(★)∧(★dec)` | **[T]** (given DS) | **zero tower shells**: window `(0,1/72]` disjoint from `[2.5,∞)`; requires DS (imported) |
| **all `j`, SZZ window `β<1/48`** (std. Wilson `β_W ≈ 0.19`) | SZZ LSI (`K_S = N(½−24\|β\|)`) → Lemma SC | **[T]**, `N`-uniform | **zero tower shells** (`β_W≈0.19` vs `β_j≥0.37`); strong-coupling (high-temp) only |
| **`j` with `β_j < 0.19`** | (same, sharpened) | **[T]** off-tower | **vacuous on the tower**; `α(2.5)≈180` `[NG]` for the method class |
| **all `j`, any `β`** | `N`- and `a`-bookkeeping (`δ_W=c_wN`, `a⁰`) | **[T]** mechanism / **[C]** synthesis | says nothing about the β-window; the defect itself is `[O]` |
| **small-field `𝒮_t`, all `j`** | Lemma S1 / Cor S1′: `ρ(U′)` modulus, `D≤0` | **[T]** / **[C]** (on `C₁`) | gives the **fibre half `D≤0`** only, not the coarse gap `E[Σ]≥δ_W` |
| **deep UV `j ≥ j_UV = O(1)`** | BL equipartition `Var≤E[H]` ⇒ `W_j≥δ_W` (≡ `E[Σ]≥δ_W`) | **[C]** on `E\|F_j\|≤Cg_j` `[O]` | all but `O(1)` shells, but **conditional**; only on `(★_sf)`/cone-empty domain, **not all `U′`** |
| **cross-shell, all `j`** | `(★dec)` band decay / Gershgorin conjunction | **[O]** (`H2-global` fails at head) | the tower `min`; independent of `(★)` (`S4` `[T]`) |
| **head `j=1 … j_head`** | — | **[O]** = **H1** | the gap itself; all-`U′` form `[NG]`-refuted (`G4`: `G″<0` at `β>4.09`) |
| **intermediate `β∈(0.19, 2.5)`** | — | **not inhabited by the tower** | if a shell lands here (1-loop low index) it is head territory |
| **covariant `P`, `j≥j_UV` / head** | `G″_P = 2W_geo + R_proj` | **[C]** (`O(g_j)`) / **[O]** | `R_proj` unbounded by `W_geo` at head (`G3`) |

**Union.** `0 ∪ [j_UV, J] = all shells except `j = 1 … j_UV−1``. The union is **not** a union of
theorems: the strong jaw contributes `∅`, and the UV jaw is `[C]` on an `[O]` bootstrap.

---

## 6 · Gap analysis (quantified)

Witness `research/v3/w7/v4_pincer.py` (deterministic; `py -3.12`; output quoted):

```
N = 2 : beta_j = 0.3715 j      beta_1 = 0.3715
        shells j>=1 with beta_j < 1/72 = 0 ;  < 0.19 = 0
        beta_j = 2.5 at j = 6.73 ;  cone-empty j_UV = 1..3 (canonical kappa=cL=1 -> 2)
        tower J=6,10,14 : |SZZ shells|=0, |UV shells|=J-1, uncovered=1 (j=1)
N = 3 : beta_j = 0.8359 j      beta_1 = 0.8359
        shells j>=1 with beta_j < 1/72 = 0 ;  < 0.19 = 0
        beta_j = 2.5 at j = 2.99 ;  cone-empty j_UV = 1..2 (canonical -> 2)
        tower J=6,10,14 : |SZZ shells|=0, |UV shells|=J-1, uncovered=1 (j=1)
```

**Read-off, in order of what the task asks.**

1. **Covered shell range — strong coupling:** **∅** for both `N`, both thresholds. The thresholds are
   crossed at `j < 1`; even `β_1 = 0.37–0.84 > 0.19 > 1/72`. `[T]` arithmetic.

2. **Covered shell range — UV:** `j ≥ j_UV`, `j_UV ∈ {1,2,3}` across an `O(1)` scan of the two cone
   constants (`κ, c_L`); canonical value `j_UV = 2`. This is `J − j_UV + 1` shells. `[C]` on
   `E|F_j|≤Cg_j` (`[O]`) and on `(★_sf)`.

3. **The gap between the strong-coupling and UV ranges (in `β`).** The uncovered `β`-interval is
   `(β_SZZ, β_{j_UV}) ≈ (0.19, O(1)·a₀^{...})`. In the **physical** convention the tower has **no**
   shells at `β < 2.5`, so this interval is **outside the tower** entirely: it is a gap in `β`-space,
   **not** an uncovered band of tower shells. In the **one-loop** convention the lowest shells do
   extrapolate into `β∈[0.37,2.5]`; those shells are precisely the head, so the `β`-gap and the head
   **coincide**. Either way the uncovered tower set is `j = 1 … j_UV − 1`, a **finite,
   `J`-independent** band (`1` shell on the literal index; `j ≲ j_head ≈ 3` (`N=3`) / `7` (`N=2`) in the
   physical-head convention — the offset between the two counts is the one-loop-vs-physical offset
   flagged in §1).

4. **Is the union "all but one shell / one region"?** Yes, structurally: **all but `O(1)` shells**, and
   the remainder is **one region** — the head. **But the union is conditional, and the head is the
   Clay core.** So "all but the head" is the correct *geometric* description of the pincer, while the
   *logical* statement is `[C]`-on-UV `∧` `[O]`-at-head, not a theorem.

5. **Does this reduce Clay to the head?** **For existence, yes; for the gap, no.** Closing the pincer
   gives `G1` fixed-volume existence (OS0/1/3, subsequential `μ_L`) and **nothing** about the spectrum
   (`G1-L3.md` §5; `C2-BALABAN-UNIFORM.md` §4). The head *is* `H1 = G2`. So the pincer does **not**
   narrow Clay; it narrows *existence* to `(★_sf)∧(★dec)` and isolates that the surviving obstruction
   is the gap read in curvature language.

---

## 7 · Honest verdict — labels

| claim | label | basis |
|---|---|---|
| S6 Dobrushin `(★)∧(★dec)` at `β_j<1/72` | **[T]** (given DS) | `S6` §2 |
| SZZ `N`-uniform gap/LSI at `β<1/48` | **[T]** | `C1` §§1–3; `PROOF-COARSE-POINCARE` (c) |
| strong-coupling window `(0,1/72]` / `(0,0.19]` vs tower `[2.5,∞)` disjoint | **[T]** (arithmetic) | `S6` §4; `v4_pincer.py` |
| covered strong-coupling tower shells | **∅** | `v4_pincer.py` |
| UV BL equipartition on cone-empty shells | **[C]** on `E\|F_j\|≤Cg_j` `[O]` | `S6` §5; `G1-L3-UV-FORCE-BOUND` §3,6 |
| UV jaw covers `j ≥ j_UV`, `j_UV=O(1)` | **[V]** (index) / **[C]** (inequality) | `v4_pincer.py`; `S1` Cor S1′ |
| `N`/`a` uniformity, `δ_W=cN`, `Θ(a⁰)` | **[T]** mechanism / **[C]** synthesis | `S5` §§3–5 |
| small-field log-concavity `D≤0` on `𝒮_t` | **[T]** / **[C]** (on `C₁`) | `S1` §3 |
| head `G″≥δ_W` (coarse, at RG min) | **[O]** = **H1** | `S1` §6; `S3`; `G4` §3; `PROOF-JACOBIAN-SAVES` §2 |
| uniform-in-all-`U′` `(★)` | **[NG]** (disproved) | `S1` Prop U′-NG; `G4` §3(iv) |
| `(★)` and `(★dec)` independent | **[T]** | `S4` |
| covariant projection term `R_proj` | `O(g_j)` `[C]` UV / `O(W_geo)` `[O]` head | `G3` §6 |
| **"(★) holds on all tower shells except the head", as a theorem** | **NO — `[C]`/`[O]`** | this note §0, §5, §6 |
| **"closing the pincer proves the gap"** | **NO — it proves existence only** | `G1-L3` §5; `C2` §4 |

**The one-line verdict.** The pincer is a genuine **conditional** reduction, not a theorem: the
strong-coupling jaw is a theorem that misses the tower, the UV jaw reaches all but `O(1)` shells but
only on the small-field/excised domain and only modulo the `[O]` moment bootstrap and `(★dec)`, and the
head is a finite, tower-depth-independent band that *is* `H1`. There is **no whole intermediate regime**
inside the tower — but there is also **no theorem on "all shells except the head"**. The precise
statement the campaign can carry is: `(★)` on the tower reduces to `[UV: E|F_j|≤Cg_j + (★dec) on
(★_sf)]` `∧` `[head: coarse-gap H1]`, with the strong-coupling legs `[T]` but off-tower. **No QED.
Clay 0 %.**

---

## 8 · Witness (`research/v3/w7/v4_pincer.py`)

Deterministic, numpy-only, `py -3.12`, exit 0. It evaluates the campaign's one-loop tower
parametrisation (`g_j²=1/(2β₀j)`, `β_j=4Nβ₀j`, `β₀=11N/(48π²)`), the S6/SZZ thresholds (`1/72`,
`0.19`), and the cone-empty crossover (Savvidy negative term `~κ g_j²` vs box floor `~c_L/a_j²` with
`a_j=a₀L^{-j}`), scanning the two `O(1)` cone constants over `{0.5,1,2}²`. Output quoted in §6. It is
**bookkeeping, not a proof**: the shell counts are exact given the parametrisation; the thresholds and
`κ, c_L` are the quoted physical inputs; the conclusion (strong jaw `∅`, UV jaw `j≥j_UV=O(1)`, head
`= j<j_UV` finite and `J`-independent) is invariant across the scan. The inequality content is
`§2` (S6), `§3.2` (BL/S1), `§4` (H1); the witness locates, not proves.

---

## 9 · Verification ledger

| item | checked | status |
|---|---|---|
| Dobrushin–Shlosman complete analyticity; van den Berg–Maes; Brascamp–Lieb (JFA 22 (1976) 366); Prékopa; Holley–Stroock | standard tools, carried by `S6`/`S1`/`S5` verification ledgers | **[standard]** (not downloaded) |
| SZZ CMP **400** (2023) 805, arXiv:2204.12737; `K_S=N(½−24\|β\|)`, radius `1/48`, `N`-free | `C1-SZZ-NSCALING.md` §7 (verified against e-print there) | **[V]** (carried) |
| `1/48 ↔ β_W≈0.19 ↔ Hooft λ≈96`; honest `≈13×` | `S6` §2 table; `ATOM-strong-coupling-extend` §1 | **[V]** (carried) |
| Balaban CMP 109/116/119/122; action-level `a`-uniform, not Hessian | `C2-BALABAN-UNIFORM.md` §7 (zbMATH-verified there) | **[V]** (carried) |
| `W_geo≥0` (Prékopa), `Φ″≥N/12`, Lemma B | `PROOF-JACOBIAN-SAVES` §2; `S5` §3 | **[T]** (carried) |
| `G″=2W_geo` (constant `v`), `G″=2W_geo+Δ` (varying `v`), `R_proj` | `G1-BLOCK-IDENTITY.md` §3–4; `G3` §3–4; residuals `≤1e-7` | **[T]**+**[V]** (carried) |
| `G″<0` counterexample `G″=0.206−0.0123β²`, `β*>4.09` | `G4-GENERAL-BOUND.md` §3(iv), witness `g4_general_bound.py` | **[T]** model + **[V]** (carried) |
| shell-index counts (`β_1`, `j_head`, `j_UV`, `∅` strong coupling, `J−1` UV) | this session, `v4_pincer.py` | **[V]** |

**Notes.** (i) `V1/V2/V3-*.md` do not exist; their roles are mapped in the header and §3. (ii) The
one-loop index and the physical `β_L≈2.5` bottom are offset by a few shells (§1); every conclusion is
stated so that it holds in both. (iii) `G2-TWO-PLAQUETTE.md` has unfilled numeric placeholders
(`PLACEHOLDER_*`) and was not used for any load-bearing claim. (iv) `G4` supersedes the S2 reduction
`(★)⇔2W_geo≥δ_W` for the *strong* form; the pincer's UV conclusion (`E[Σ]≥δ_W` with Haar floor on the
log-concave set) is unchanged by that correction, as noted in §3.2.

*Refs.* `S6-DOBRUSHIN.md` §0–§5; `S1-CONDITIONAL-POINCARE.md` §1,§3,§4.1,§6; `S5-N-UNIFORM.md` §2–§5;
`C1-SZZ-NSCALING.md` §1–§3; `C2-BALABAN-UNIFORM.md` §1,§2.1–2.4,§4; `STAR-WAVE-SYNTHESIS.md` §1–§4;
`G1-BLOCK-IDENTITY.md`; `G3-PROJECTION-TERM.md` §4–§6; `G4-GENERAL-BOUND.md` §2–§5;
`research/unconditional/G1-L3.md` §1–§2,§5; `G1-L3-TOWER-POINCARE.md` §1–§5;
`G1-L3-UV-FORCE-BOUND.md` §1–§6; `PROOF-COARSE-POINCARE.md`; `PROOF-JACOBIAN-SAVES.md` §1–§3;
`S3-UV-FESHBACH.md`; `S4-DECIMATION.md`; witness `research/v3/w7/v4_pincer.py`. No third-party sources
downloaded.
