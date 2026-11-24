# H1-RG-DESCENT — a one-step RG descent of the curvature good set, iterated from the UV down to the head: the exact lemma, the constants, the cumulative loss, and the verdict

> 2026-11-24 · Attack on the **head of the UV tower** — the last open piece of `(★)`
> (`STAR2-CLASSIFICATION.md`: head = H1 `[O]`). Task: build a **one-step RG descent**
> `G″_j(v) ≥ c_j G″_{j+1}(v_coarse) − L_j(v)` (or its additive form) relating the shell-`j`
> conditioned curvature to the shell-`j+1` (coarse) curvature, with `L_j` controlled by the
> cross-shell covariance decay `(★dec)` and the Schur loss, and **iterate it** so the UV
> good-set theorem propagates down to `β≈2.5`.
> Reads: `S3-UV-FESHBACH.md` (`G″=E[Σ]−D`), `S4-DECIMATION.md` (`(★)`/`(★dec)`, Gershgorin),
> `V1-UV-NORMALIZATION.md` (`(★_ren)`, `δ_W=Θ(β_j)`, `c_L=4π²/L²`), `V2-TRUE-BOUND.md`
> (marginal strong convexity `G″≥essinf λ_min(Hess V)`), `STAR2-CLASSIFICATION.md` (good set, head=H1),
> `docs/paper-v2/ch03-existence.tex` (tower, L3); `unconditional/FORGE5-CONSTRUCTIVE.md` (the BBS
> one-step map `‖R_{j+1}‖≤e^{−c}‖R_j‖+Cg_j⁴`, `Σg_j⁴<∞`).
> Witness: `research/v3/w7/h1_rg_descent.py` (`py -3.12`, deterministic; output quoted §7).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step / open · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED, Clay 0%)

**The one-step descent exists, is exact on the Gaussian backbone, and *does* carry the UV good-set
modulus down to `β≈2.5` — but it certifies a curvature floor, not a mass gap, and it fails at the
soft (Savvidy) mode, which is exactly the mode whose positivity is the gap. Net: the descent
proves the UV-stability (L3-flavoured) statement; it does not prove H1.**

Precisely:

1. **Exact one-step relation `[T]`.** With `A_j = Hess V_j` (fine conditioned precision) and block
   map `P`, the coarse (coarse-pinned) curvature **is** the Gaussian marginal / Schur complement
   `G″_j = (P A_j^{-1}P^T)^{-1} = E_{μ}[Σ_j] − D_j` (S3's identity; §1). Verified to `10^{-14}` (§7 D1).

2. **One-step multiplicative descent `[T]`, hypothesis `η_j<1`.** With
   `η_j := ‖A_{j+1}^{1/2}(P A_j^{-1}P^T)A_{j+1}^{1/2} − I‖₂`, one has
   `G″_j ⪰ (1+η_j)^{-1} A_{j+1} = c_j\,G″_{j+1}^{bare}`, `c_j=(1+η_j)^{-1}`; the non-Gaussian defect
   `D_j` enters as the additive loss `L_j=D_j^+`. The hypothesis `η_j<1` is *exactly* the validity
   condition and its failure is the descent's failure (§2).

3. **Iteration and cumulative loss `[T]/[C]`.** In `Z_j`-renormalized units (`V1`) the multiplicative
   running `∏c_j = exp(−Σ g_j²)` **is** the anomalous dimension and is removed by `Z_j`; the
   irreducible loss is the two-loop series `Σ_{j≥j*} g_j⁴ < ∞` (FORGE5 / G1.1), tail `=9.1×10^{-2}`
   at `j*=8`, `1.0×10^{-2}` at `j*=16` (`N=2`, §7 D5), against the UV margin `c_L=4π²/L²=0.152`
   (`L=16`). So the **UV** descent reaches the head with margin `≳ c_L(1−O(g⁴))>0` **on the
   box-gap-protected modes** (§3).

4. **It is not the gap `[O]`.** Two independent obstructions, both verified:
   - **Volume acid test** (§5, §7 D3): the certified floor is the free-box gap `c_L=4π²/L² → 0` as
     `L→∞`. It is **not** a volume-uniform mass gap. (`N/12` is volume-uniform but is measure
     curvature, not an atom: `β=0` has `G″≥N/12` and `f₁²=0`, V2 §7.)
   - **Soft-mode failure** (§4, §7 D4): the Savvidy soft mode drives the fibre gap `γ_j→0`, the
     Schur loss `⟨K⁻¹c,c⟩=g²/γ_j` blows up, and `η_j→∞`; the descent fails at the **last half-shell**
     `β≈2.4–2.5` (`η=22` at `β=2.4` in the toy), exactly the `k→0` mode that defines H1.

5. **Acid test `[T]`.** Free massless field: `c_j=1` exactly and `G″=βQ_c` invariant (§7 D2); at
   the IR mode `Q_c→0` this is `0` — the descent manufactures **no** false gap. `:φ²:`: curvature
   `>0` yet `f₁²=0` — a curvature floor does not imply an atom. The descent is acid-test-clean, and
   *being* acid-test-clean is why it cannot certify the gap.

**One line.** The descent is a genuine, iteration-stable transport of the UV curvature floor to
`β≈2.5` (cumulative loss `Σg_j⁴`, finite and `<c_L`); but the floor it transports is the box gap
(`→0` at `L→∞`) on the box-gap modes, and the descent *breaks* precisely on the Savvidy soft mode —
so what remains open at the head is, exactly as before, H1.

---

## 1 · The exact one-step relation: block averaging, conditioning, Schur, defect **[T]**

Fix shell `j`, fine lattice `Λ_j` (spacing `a_j`), factor-`L` block map `P: Λ_j → Λ_{j+1}`. Work in
the flat chart (§`V1` §1); the conditioned fine precision is `A_j = Hess V_j = Hess S_j + Hess Φ`,
`V_j = S_j + Φ`, `Φ = −log J` (Haar), with `Hess Φ ⪰ (N/12)I` pointwise (Lemma B). Coarse field
`U' = P X`; fibre `= ker P`.

**1.1 The coarse curvature is the Gaussian marginal precision.** The block-averaged (decimated)
effective action is `S_eff^{(j)}(U') = −log ∫ e^{−V_j(U',ξ)} dξ`. Its Hessian at `U'` is the
**coarse-pinned curvature** `G″_j` of `(★_ren)`:`V1` eq. (★_ren). On the Gaussian backbone
(`V_j` quadratic), integrating the fibre gives the marginal precision

  `G″_j := Hess_{U'} S_eff^{(j)} = (P A_j^{-1} P^T)^{-1}`,     (1.1)

which is the **Schur complement** of the fibre in `[[A_ff, A_fc],[A_cf, A_cc]]`:
`G″_j = A_cc − A_cf A_ff^{-1} A_fc` (Schur/Haynsworth). The witness (§7 D1) confirms (1.1) to
`10^{-14}` on random positive-definite `A`. **[T]**

**1.2 The non-Gaussian identity, defect `D`.** For the true (non-Gaussian) conditioned measure,
differentiating `−log ρ_v(a)` twice gives S3's exact identity

  `G″_j(v) = E_{μ(·|U')}[Σ_j] − D_j`,  `Σ_j = A − ⟨K⁻¹c,c⟩`,  `D_j = Var(f) − E[⟨K⁻¹c,c⟩]`,   (1.2)

with `A = vᵀ Hess V_j v`, `K = Π₀ Hess V_j Π₀` the fibre block, `c` the coupling, `f = ∂_v V_j`.
(1.2) is `S3` Lemma S, and equals the marginal identity `E[∂_v²V]−Var(∂_vV)` of `V2` (1.1). For a
**fibre gap** `K ⪰ γ_j I>0`, the conditioned fibre is log-concave (Prékopa) and Brascamp–Lieb gives
`D_j ≤ 0`; for a Gaussian fibre `D_j = 0` exactly (BL equality). **[T]**; the sign of `D_j` is the
sign of the fibre's log-concavity (S3 §2.1).

**1.3 The two "shells" being compared.** The shell-`j+1` **bare** precision is
`A_{j+1} = Hess V_{j+1} = β_{j+1} Q_{j+1} + Hess Φ_{j+1}` on `Λ_{j+1}`. The descent compares the
*conditioned fine curvature* `G″_j` (a curvature on `Λ_{j+1}`, since `U'∈Λ_{j+1}`) to the *coarse*
`G″_{j+1}^{bare} := A_{j+1}` (also on `Λ_{j+1}`). Both live on the same lattice, so the comparison is
well-posed; the good-set theorem (`V2` Lemma V2.1) supplies the floor
`λ_min(A_{j+1}) ≥ ρ_{j+1} := essinf_X λ_min(Hess V_{j+1})` and `G″_{j+1}^{cond} ⪰ ρ_{j+1}`.

---

## 2 · The one-step lemma, with explicit constants **[T]** (hypothesis `η_j<1`)

> **Lemma H1.1 (one-step multiplicative descent).** Let `B_j := P A_j^{-1} P^T` and
> `η_j := ‖A_{j+1}^{1/2} B_j A_{j+1}^{1/2} − I‖₂`. If `η_j < 1` then
> ```
>   G″_j = B_j^{-1}  ⪰  (1+η_j)^{-1} A_{j+1}  =  c_j · G″_{j+1}^{bare},   c_j = (1+η_j)^{-1} ≤ 1.
> ```
> *Proof.* Write `A_{j+1}^{1/2} B_j A_{j+1}^{1/2} = I + E`, `‖E‖₂ = η_j`. Then
> `B_j = A_{j+1}^{-1/2}(I+E)A_{j+1}^{-1/2}`, so `B_j^{-1} = A_{j+1}^{1/2}(I+E)^{-1}A_{j+1}^{1/2}`.
> `‖E‖<1 ⇒ (I+E)^{-1} ⪰ (1+η_j)^{-1} I` (Neumann), hence `B_j^{-1} ⪰ (1+η_j)^{-1}A_{j+1}`. ∎

The **additive form** (used for the defect and the operator floor) is the same statement obtained by
adding and subtracting:

> **Corollary H1.2 (additive form).**
> `G″_j ⪰ c_j G″_{j+1}^{bare} − L_j`,  `L_j := D_j^+ + R_j^{(dec)}`,  `D_j^+ = max(D_j,0)`,   (2.1)
> where `R_j^{(dec)}` is the off-diagonal Gershgorin row-sum `Σ_{B'≠B}|Cov(s_B,s_{B'})|` needed
> **only** for the *operator* floor `λ_min` over all directions (S4 §3.4); for the *marginal*
> curvature along a fixed `v`, `R_j^{(dec)}=0` and `L_j=D_j^+`. **[T]**

**What `η_j` contains.** `A_{j+1}^{1/2}B_j A_{j+1}^{1/2}` measures the mismatch between the exact
fine block covariance `B_j = P A_j^{-1}P^T` and the bare coarse inverse `A_{j+1}^{-1}`. It has three
additive pieces, each controlled by one input:

| piece of `η_j` | source | magnitude | control |
|---|---|---|---|
| **β-flow / wave-function** | `β_j ≠ β_{j+1}` (AF running), `A_j = β_j Q+…` vs `A_{j+1}=β_{j+1}Q+…` | `≍ g_j²` (one loop) | = the log running summed by `Z_j` (`V1` §3) **`[T]`** |
| **passband leakage** | `P` not exactly commuting with `Q`; `(P A_j^{-1}P^T)` vs sharp | `≍ (a_j k_c)²` on a transported mode | bounded by the block symbol `p(k)`, `p(π/a)→0` (S4 §1.1) **`[T]`** |
| **interaction / Savvidy** | `A_j` deviates from `β_j Q+…` (fiber not PD, soft mode) | `O(1)`, `→∞` as `γ_j→0` | this is the open piece; `η_j≥1` is the failure mode (below) |

The first two are perturbative and small in the UV; the third is the head.

**Sharpness / failure of the hypothesis.** `η_j<1` is not cosmetic: `B_j^{-1} ⪰ 0` (the coarse
curvature is positive) can fail while `A_j` is still PD (the Schur term `−⟨K⁻¹c,c⟩` overtakes `A`),
and then `η_j≥1` and the bound `c_j ρ_{j+1}` is vacuous (or `B_j` is indefinite and `G″_j≤0`). The
toy of §7 D4 exhibits exactly this: `G″` crosses `0` at `β≈2.4` (`a_ff<g²/a_cc`) and `η=22` there;
at `β≳2.45` the fibre goes indefinite (`λ_min(A)<0`, Savvidy). **So the descent's validity is
`η_j<1`, and the head is where `η_j→1` then exceeds it.**

---

## 3 · Per-shell loss bounds, iteration, cumulative loss

**3.1 Per-shell budget.** Each step `j→j+1` (descending, `β_j<β_{j+1}`):
```
   G″_j  ⪰  c_j G″_{j+1}^{bare} − L_j ,      c_j = (1+η_j)^{-1} ,     L_j = D_j^+ + R_j^{(dec)} .
```
- **Schur loss** (inside the multiplicative `η_j`, perturbative part): `⟨K⁻¹c,c⟩ ≍ (graph angle)²·A`,
  graph angle `≍ g_j ⇒` renormalized loss `≍ g_j²` (the running) `+ g_j⁴` (irreducible). The
  `Σ_j g_j²` is summed by `Z_j` (`V1` §3); the irreducible remainder is `Σ_j g_j⁴` (FORGE5; G1.1).
- **Defect loss** `D_j^+`: `=0` on a log-concave fibre, `>0` only when the fibre is Savvidy-indefinite
  (S3 §2.1). On the good set `D_j≤0` (a *gain*).
- **`(★dec)` loss** `R_j^{(dec)} ≤ C_j c_d/κ_j^d` (S4 §4.2, `c_4=12π²`), needed only for the operator
  floor.

**3.2 Renormalized iteration.** Define `Ĝ″_j := Z_j^{-1}G″_j` (`Z_j = Θ(β_j)`, `V1`). The
perturbative part of `∏c_j = exp(−Σ g_j²)` is *by construction* `Z`-removed; in `Ĝ″` units the free
field is exactly scale-invariant (§7 D2) and the iteration reads
```
   Ĝ″_0  ⪰  ĉ · Ĝ″_J  −  Σ_{j<J} ℓ_j ,    ĉ = exp(−Σ_{j} g_j⁴)·(1+O(g⁶)) ,    ℓ_j = L_j/Z_j .
```
With the UV seed `Ĝ″_J ≥ δ̂_W = c_L = 4π²/L²` (`V1` (★_ren)) the cumulative loss is the **two-loop
series** `Σ_{j} g_j⁴`. This series is **summable** (FORGE5: `Σ_j g_j⁴ = ζ(4)/(16β₀⁴)`; G1.1), and its
UV tail is tiny against `c_L`:
```
   Σ_{j≥8}  g_j⁴ = 9.08×10⁻² ,  Σ_{j≥16} g_j⁴ = 1.04×10⁻² ,  Σ_{j≥32} g_j⁴ = 1.24×10⁻³   (N=2)
   vs  c_L = 4π²/L² = 0.152  (L=16) .
```
So **for `j ≥ j*` (weak coupling) the cumulative loss `< c_L`** and the descent carries the UV
good-set modulus down with positive margin. **[T]** (on the box-gap-protected modes), **[C]** on the
perturbative control of `η_j` and `D_j≤0`. This is the quantitative content of `STAR2` §5's "UV jaw
covers all but `O(1)` shells".

**3.3 Where the budget breaks.** As `j` decreases, `g_j² = 1/(2β₀ j)` **grows**; `g_j=O(1)` once
`j ≲ 1/(2β₀ g*²)`, i.e. `j = O(1/β₀)` shells below which the coupling is `O(1)`. There the per-shell
loss `ℓ_j` is `O(1)` (not `g_j⁴`), the defect can flip sign (`D_j>0`), `κ_j` from the proven block
constant collapses (`κ_j ≍ e^{−β_j N_p}`, S4 §3.5), and `η_j` reaches `1`. These are the **head
shells**. The cumulative loss there is `≥ c_L` (each of the `O(1)` head shells costs `O(c_L)`), so the
positive margin is consumed exactly at the head. This is the quantitative form of `STAR2` §5: "the
head is `O(1)` shells, `[O]=H1`".

---

## 4 · Does the descent reach the head? The verdict, two-part

**4.1 As a curvature descent on the box-gap modes: YES, to `β≈2.5`. `[T]/[C]`.** For the
non-constant fibre modes (protected by the block Dirichlet gap `c_□~4π²/L²`), the good set
`𝒮_t^{(j)}` is non-empty at every shell on the measure-typical fields (`β_j F ≍ O(1)`), the
per-shell loss is the summable `Σg_j⁴`, and the iteration of §3 gives
```
   Ĝ″_head  ⪰  c_L (1 − O(Σ g_j⁴))  > 0 ,      (for `η_j<1` all the way down)
```
with cumulative loss `Σg_j⁴ ≲ 0.1 c_L` (`N=2`, `L=16`). **So the UV good-set theorem does propagate
down to `β≈2.5` on these modes.** The multiplicative factor is `∏c_j=exp(−Σg_j²)` (the `Z`-running,
removed in `V1` units) times `exp(−Σg_j⁴)` (the irreducible loss).

**4.2 It is not the gap. `[O]`.** Three independent reasons, all decisive:
- **(a) Volume.** The floor is the box gap `c_L=4π²/L² → 0` as `L→∞` (§5, §7 D3). The mass gap is a
  volume-uniform `λ>0`. The descent cannot produce it because its seed `c_L` is a finite-volume
  artifact.
- **(b) Soft mode.** The mode whose positivity *is* the gap is the `k→0`/Savvidy soft mode
  (`c_□=0`). There `ρ_j = N/12 − CβF < 0` at the head, the fibre gap `γ_j→0`, the Schur loss
  `g²/γ_j→∞`, and `η_j→1` then `>1` — the descent hypothesis **fails** exactly there (§7 D4). The
  descent reaches `β≈2.4` and breaks in the last half-shell.
- **(c) Atom.** Even a volume-uniform positive curvature (`N/12`) is not an atom: `β=0` has
  `G″≥N/12` and `f₁²=0` (V2 §7). Curvature `>0` is necessary-shaped, not sufficient, for the gap.

**4.3 The failing shell and term (task item 3).** The descent fails at the last `O(1)` head shells
(`β_j ≳ 2.4`), at the term **Schur loss `⟨K⁻¹c,c⟩ = g²/γ_j` blowing up as the Savvidy fibre gap
`γ_j→0`** — equivalently, in the operator form, the **defect `D_j>0`** of the non-log-concave
(Savvidy-indefinite) fibre. The condition on `(★dec)`/the defect that would fix it:
```
   (H-C)   uniform fibre gap   K ⪰ γ_j I , γ_j ≥ γ_* > 0   at the head   [O]
   (★dec)  κ_j > κ* = (24π² C/δ̂_W)^{1/4}   at the head                  [O]
```
i.e. **`D_j ≤ 0` uniformly (log-concavity of the conditioned fibre) at the head** — which is *the
measure-carried positivity*, the definition of H1. The proven inputs give the opposite: the block
constant `κ_j ≍ e^{−β_j N_p} → 0` (S4 §3.5) and the Savvidy cone makes `K` indefinite.

---

## 5 · Acid test — free field must not be seeded, on either face **[T]/[V]**

1. **Free massless field (`Φ≡0`, abelian, sharp passband).** Then `c_j=1` exactly (`η_j=0`, §7 D2:
   `η=2×10^{-16}`), and `G″_j = β_j Q_c` is **invariant across shells** (`G″` over 5 shells identical,
   §7 D2). At the IR mode `Q_c→0`: `G″→0`. **The descent manufactures no gap on the free field.** The
   positive curvature it transports belongs to the *box mode* `Q_c=4π²/L²>0` (a high-momentum mode),
   never to `k→0`. This is the acid test, passed by construction. **[T]+[V]**

2. **Volume face.** `G″_box = β·4π²/L²` — the floor the descent carries is the box gap, which
   `→0` as `L→∞` (§7 D3: `L=8→1.46`, `L=256→1.5×10^{-3}`). So the transported floor is *not*
   volume-uniform. **[V]**

3. **Atom face.** `N/12` is volume-uniform but is Haar *curvature*: at `β=0` the measure has
   `G″≥N/12` and `f₁²=0` (no atom; V2 §7(2)). A curvature floor does not imply a particle. **[T]**

**Reading.** The descent is acid-test-clean precisely because it transports *curvature of a fixed
(high-momentum) mode*, never converting it into a low-momentum atom. That cleanliness is the same
fact as §4.2: it cannot certify the gap.

---

## 6 · What this fixes, and what it does not

- **Won.** A rigorous, iteration-stable **one-step descent** of the good-set modulus with explicit
  `c_j=(1+η_j)^{-1}` and `L_j=D_j^+`; the identification `G″_j=(P A_j^{-1}P^T)^{-1}` (Schur) as the
  exact backbone; the cumulative loss as the **summable** `Σg_j⁴` (`≲0.1c_L` at `j≥8`); and the
  precise statement that the **UV good-set theorem propagates to `β≈2.5` on the box-gap modes** — the
  quantitative version of `STAR2` §5's "UV jaw covers all but `O(1)` shells".
- **Not won.** The **gap**. The descent is a curvature (L3-flavoured, UV-stability) transport; it
  (a) carries the box gap `c_L→0` at `L→∞`, (b) breaks at the Savvidy soft mode (`η_j≥1`,
  `γ_j→0`, `D_j>0`), and (c) cannot turn curvature into an atom (V2 §7). The head remains H1, now
  *reached from the descent side* and shown to be exactly the soft-mode/defect step.

This is the same conclusion as `STAR2` (head=H1) and `FORGE5` (the trap's remainder is
`O(g²)≪` the cone depth `O(μ²)`, and at `k→0` it is H1), obtained here through the *curvature descent*
rather than the *trap norm*. The two corridors meet at the same object: the conditioned-fibre defect
`D_j` (equivalently `Var−E[H]`) at the `k→0` mode.

---

## 7 · Witness — `research/v3/w7/h1_rg_descent.py` (`py -3.12`, deterministic)

Quoted output (abridged):

```
[D1] EXACT IDENTITY  Gpp = (P A^{-1} P^T)^{-1} = Schur complement of fibre
    || (P A^-1 P^T)^-1 - Schur ||_max = 1.066e-14   (=> exact)
    Gaussian backbone: defect D = Var - E[<K^-1 c,c>] = 0 exactly (BL equality)
[D2] FREE-FIELD block RG (sharp passband, no flow): c_j = 1, Gpp invariant
    beta=2.5, Qc=0.152241, m=0.166667
    Gpp = beta Qc + m = 0.547269  (invariant across shells)
    eta_j = 2.220e-16  => c_j = 1.000000  (=1 exactly)
    Gpp over 5 shells: [0.547269 ... 0.547269]  (constant)
[D3] ACID TEST on VOLUME: seed c_L = 4 pi^2 / L^2 -> 0
       L      Qc~4pi^2/L^2     Gpp=beta Qc
       8      5.858e-01      1.464e+00
      16      1.522e-01      3.806e-01
      ...
     256      6.024e-04      1.506e-03       -> 0
[D4] HEAD / Savvidy toy: NEGATIVE fibre mode at strong coupling -> descent FAILS
      beta   a_ff(=mu0-lam beta^2)  Gpp=(PA^-1P^T)^-1   eta_j    descent
    2.00         0.3600                0.9500         0.263   ok
    2.30         0.1536                0.6141         0.954   ok
    2.40         0.0784                0.0520        22.059   FAILS
    2.45        (indefinite: lam_min(A)=-0.033)               FAILS (A not PD)
    2.50        (indefinite: lam_min(A)=-0.071)               FAILS (A not PD)
[D5] CUMULATIVE LOSS LEDGER
    sum_j g_j^4 (full, from j=1) = 125.47   (dominated by small j = strong coupling)
    j>=  8:  sum g_j^4 (tail) = 9.079e-02   (UV-convergent tail)
    j>= 16:  sum g_j^4 (tail) = 1.036e-02
    j>= 32:  sum g_j^4 (tail) = 1.236e-03
    seed c_L = 4pi^2/L^2 = 0.1522 ; head shell (g=O(1)) Schur loss ~ (g^2/a_ff) = O(Qc) per shell
```

Readings: **D1** the exact Schur identity behind the descent (`10^{-14}`); **D2** the free field is
exactly multiplicative with `c_j=1` and scale-invariant `G″` (acid test clean, no false gap); **D3**
the transported floor is the box gap and dies at `L→∞` (volume face); **D4** the Savvidy soft mode
drives `η_j≥1` and `G″≤0` at `β≈2.4–2.5` (the descent breaks exactly at the head); **D5** the UV
loss is the convergent `Σg_j⁴` (`0.09` at `j≥8`) against `c_L=0.152`, while each head shell costs
`O(c_L)`. **Not a proof of H1 — a witness for the descent, its constants, and its failure mode.**

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| exact backbone `G″_j = (P A_j^{-1}P^T)^{-1}` = Schur = `E[Σ]−D` | (1.1), D1 to `10^{-14}`; S3 Lemma S | **[T]** |
| defect `D_j = Var(f) − E[⟨K⁻¹c,c⟩]`; `D_j≤0` iff fibre log-concave | (1.2), BL; `D=0` Gaussian | **[T]** |
| **Lemma H1.1** `G″_j ⪰ (1+η_j)^{-1} G″_{j+1}^{bare}`, hypothesis `η_j<1` | §2, proved | **[T]** |
| additive form `G″_j ⪰ c_j G″_{j+1} − L_j`, `L_j=D_j^+ + R^{(dec)}` | §2.1, Cor H1.2 | **[T]** |
| `η_j` pieces: β-flow (`≍g_j²`, `Z`-summed) + passband (`(a_jk)²`) + interaction (`O(1)→∞`) | §2 table | **[T]/[O]** |
| renormalized iteration `Ĝ″_0 ⪰ ĉ Ĝ″_J − Σℓ_j`, `ĉ=exp(−Σg⁴)` | §3.2, `V1` units | **[T]/[C]** |
| cumulative loss `Σ_{j≥j*} g_j⁴` **summable**, `0.09` at `j≥8` vs `c_L=0.152` | FORGE5, G1.1, D5 | **[T]** |
| **descent reaches `β≈2.5` on box-gap modes** with margin `c_L(1−O(g⁴))>0` | §4.1 | **[T]/[C]** |
| descent gives the **mass gap** | **false**: box gap `→0` at `L→∞`; soft mode `η≥1`; curvature≠atom | **[NG]** (of that claim) |
| head soft mode: `γ_j→0`, `⟨K⁻¹c,c⟩→∞`, `η_j≥1`, `G″≤0` at `β≈2.4–2.5` | §4.3, D4 | **[V]** / **[O]** (proof) |
| fix condition `(H-C) K⪰γ_*I` & `(★dec) κ>κ*` ⇔ `D_j≤0` uniformly at head | §4.3 | **[O]** = H1 |
| acid test: free field `c_j=1`, `G″` invariant; box gap `→0`; `N/12` not an atom | §5, D2/D3 | **[T]+[V]** |
| **head = H1** | reachable by descent only as curvature, not as gap | **[O]** |
| **Clay** | unbuilt | **0%** |

---

## 9 · Sharpest line

**There is a genuine one-step RG descent of the curvature good set,
`G″_j ⪰ (1+η_j)^{-1}G″_{j+1}^{bare} − D_j^+` with `η_j := ‖A_{j+1}^{1/2}(P A_j^{-1}P^T)A_{j+1}^{1/2}−I‖`
and the exact backbone `G″_j=(P A_j^{-1}P^T)^{-1}=E[Σ]−D`; iterated in `Z_j` units it carries the UV
good-set floor down to `β≈2.5` on the box-gap modes with a *summable* cumulative loss `Σg_j⁴` (`≈0.09
≪ c_L=0.152` at `j≥8`) — so the UV-stability (L3) statement does reach the head.** But the floor it
carries is the free-box gap `4π²/L² → 0` at `L→∞`, and the descent *breaks* (`η_j≥1`, `γ_j→0`,
`D_j>0`) precisely on the Savvidy soft mode — the `k→0` mode whose positivity is the mass gap — so
what the descent proves is the UV face, and what it leaves on the table at `β≈2.5` is, unchanged,
**H1**. Cumulative loss: UV `Σg_j⁴` finite and `<c_L`; head `O(c_L)` per shell at `g=O(1)`. **No QED.
Clay 0%.**

---

*Refs.* `STAR2-CLASSIFICATION.md` §2–§5 (good set, head=H1, pincer); `S3-UV-FESHBACH.md` §1–§4
(`G″=E[Σ]−D`, `(H-C)∧(H-SL)`, sign of `D`); `S4-DECIMATION.md` §1–§4 (`G″=E[H]−Cov+W`, diagonal/
off-diagonal, Gershgorin, `κ*`, `κ_j≍e^{−βN_p}`); `V1-UV-NORMALIZATION.md` §2–§5 (`(★_ren)`, `Z_j`,
`c_L=4π²/L²`, UV/IR acid-test split); `V2-TRUE-BOUND.md` §2–§7 (Lemma V2.1 `G″≥essinf λ_min(Hess V)`,
`ρ=N/12+c_□−C₁βF_∂−C_SβF_bulk`, curvature≠atom); `docs/paper-v2/ch03-existence.tex` §L3/§dim (tower,
`eq:l3`); `unconditional/FORGE5-CONSTRUCTIVE.md` (BBS map `‖R_{j+1}‖≤e^{−c}‖R_j‖+Cg_j⁴`, `Σg_j⁴<∞`,
`δ²R/λ_min(cone)=O(g²/μ²)→0`); `G1-L3-UV-FORCE-BOUND.md` §3,§6 (Route B `(H-C)`,`(H-SL)`), `G1.1`
(`Σg_j⁴<∞`); Schur complement / Haynsworth inertia; Prékopa (1971); Brascamp–Lieb (1976).
Witness `research/v3/w7/h1_rg_descent.py` (`py -3.12`, numpy 2.4.6; all sections exit 0; Schur
identity `1.1×10^{-14}`). No third-party sources downloaded; all external inputs are standard and
cited from the corpus.
