# SM3-ANISOTROPIC-RG — a two-step / anisotropic RG descent that keeps the Savvidy soft mode out of the fibre: the one-step inequality, the constants, `η<1` on the soft mode, the fixed-point obstruction, and the soft-mode floor

> 2026-11-24 · Attack on the **Savvidy soft mode** — the one mode on which the UV-face RG descent
> (`H1-RG-DESCENT.md` §4.2b) breaks. Task: **modify the descent so the soft mode stays controlled** —
> anisotropic blocks (long in the soft direction), a separate soft-mode variable, or a two-step
> descent with an intermediate coarse-graining — state the exact one-step inequality with constants,
> decide whether `η<1` holds on the soft mode at `β=2.5`, and if not prove the obstruction.
> Reads: `H1-RG-DESCENT.md` (Lemma H1.1, `η_j`, the `β≈2.4` break, D4), `S3-UV-FESHBACH.md`
> (`G″=E[Σ]−D`, fibre gap ⇒ `D≤0`), `S4-DECIMATION.md` (`(★dec)`, `κ*`, block constant), `docs/paper-v3/main.tex`
> §6 (`eq:onestep`, the soft-mode break, the Haar floor `N/12`, the flat-chart caveat §6.6).
> Witness: `research/v3/w7/sm3_anisotropic_rg.py` (`py -3.12`, deterministic; output quoted §7).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step / open · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED, Clay 0%)

**The modified descent DOES control the soft mode in the `η<1` sense — by retaining the soft mode in the
coarse field so the fibre holds only hard modes — but it proves a *curvature* floor on the good set, not the
gap; and a fixed-point obstruction shows that no blocking (isotropic, anisotropic, or two-step) can lift the
soft-mode curvature beyond its bare value, so what is left at `β=2.5` is the Savvidy/measure sign, i.e. `H1`
unchanged.** Precisely:

1. **Why isotropic blocks fail [T].** The soft mode is the lowest (longest-wavelength) mode. An isotropic
   block map `P` leaks a fraction `O((a k_soft)²)` of it into the fibre `ker P`; that leaked component still
   carries the *soft* curvature `ε`, so the fibre gap `γ_j ≤ ε → 0`, the Schur loss `g²/γ_j → ∞`, and
   `η_j ≥ 1` — the `H1` break (§1, `H1-RG-DESCENT` D4: `η=22` at `β=2.4`).

2. **The scheme [T]/[C].** Anisotropic blocks — *long in the soft direction*, i.e. the soft mode is retained
   at full resolution and only the hard directions are decimated (equivalently: a two-step descent that first
   cascades the hard sublattice, then coarsens the soft direction; equivalently: an explicit soft-mode
   variable integrated against the gapped hard complement). The fibre now contains **only** hard modes, with a
   gap `γ_hard ≥ c_□^{blk} = (π/L_blk)² + 2g|B| > 0`. The one-step lemma is unchanged,
   ```
   G″_j  ⪰  (1+η^{an}_j)^{-1} A_{j+1},      η^{an}_j = ‖A_{j+1}^{1/2}(P_an A_j^{-1}P_an^T)A_{j+1}^{1/2} − I‖,
   ```
   but the interaction piece of `η^{an}_j` is now `(graph angle)²·A_j/γ_hard` — **finite** — instead of the
   divergent `g²/γ_soft` (§2–§3, Lemma SM3.1).

3. **`η<1` on the soft mode at `β=2.5`: YES, conditionally [T]/[C].** With a hard fibre, `η^{an}_soft ≈ g²/γ_hard`.
   At the head `g²≍1–1.6`, `γ_hard ≳ (π/L_blk)² ≳ 2.5` (lattice units), so `g²/γ_hard ≲ 0.6 < 1` (§4, W1/W5).
   On the `H1` toy numbers (no Haar) the anisotropic scheme moves the `η<1` window from `β≈2.35` (isotropic)
   out to `β≈2.45`; the last step to `β=2.5` needs the **Haar floor** `Φ″ ⪰ N/12` to keep the bare soft
   curvature `ε ≥ N/12 > g²/γ_hard` (W5: `η^{an}_soft = 0.28` at `β=2.5`). So the soft-mode `η≥1` failure was a
   *bad-blocking artifact*; removing it exposes the real obstruction (point 4).

4. **The real obstruction — the soft mode is a fixed point of the blocking [T]/[NG].** If the soft mode `v`
   is retained (`P_an P_an^T v = v`), then for **any** precision `A_j`,
   `v^T(P_an A_j^{-1}P_an^T)v = v^T A_j^{-1}v`, so `G″_j(v) = A_j(v)` exactly (`η=0`; W2 to `7×10^{-15}`). The
   descent is, on the retained soft mode, an **identity**: it *transports* the soft curvature, it does not
   *create* positivity. Eliminating the gapped (positive) hard fibre is convention (A) and can only *lower*
   `G″` (Schur loss `≥0`, `S3` §2.4). Hence **no blocking scheme can push `G″_soft` above the fine bare soft
   curvature `A_j(v_soft)`** — `[NG]` for "an RG blocking lifts the soft-mode curvature." The only positivity
   that survives is what is already in `A_j(v_soft)`.

5. **What it proves — the soft-mode floor [T]/[C], not the gap [T]/[O].** On the good set the bare soft
   curvature is `ε_good = βQ(k_soft) + Φ″_soft + Sav_soft ≈ β c_L + N/12` (`Sav≈0` on small field). The scheme
   transports the floor
   ```
   G″_soft  ⪰  (1+η^{an})^{-1}[ β c_L + N/12 ]  −  g²(β c_L + N/12)/γ_hard  > 0     (good set, §5, W4)
   ```
   — a **β-dependent, good-set, volume-uniform (via the Haar `N/12`) soft-mode curvature floor** (`≈0.40` at
   `β=2.5`, W4). But (a) it is `N/12`-protected *curvature*, and `β=0` has `G″≥N/12` with `f₁²=0` — **not an
   atom** (`V2` §7 acid test); (b) the **unconditional essinf floor is negative at `β=2.5`**: on the Savvidy
   tail `Sav→−1.1β`, so `ε_good − 1.1β < 0` for `β > β* ≈ 0.18` (W4). The good-set↔essinf gap is exactly the
   measure-carried positivity, `H1` (§5). The floor is the flat-chart mode's curvature (gauge-dependent,
   `main.tex` §6.6), not the physical gauge-invariant gap.

**One line.** The anisotropic/two-step descent *fixes the `η≥1` symptom* (hard fibre instead of a soft one),
but the soft mode is a **fixed point** of every blocking, so the descent only *transports* the soft-mode
curvature `ε` — which on the good set is the Haar floor `N/12` (curvature, not atom) and on the Savvidy tail
is negative. So the modified descent proves a good-set curvature floor and *nothing more*: `H1` (the uniform,
i.e. essinf, positivity of the soft mode at `β=2.5`) is untouched, now shown from the blocking side to be
exactly the Savvidy sign.

---

## 1 · Why isotropic blocks fail on the soft mode — the leakage mechanism **[T]**

Recall (`H1-RG-DESCENT` §2): `A_j = Hess V_j` (fine conditioned precision), `P` the factor-`L` block map,
`G″_j = (P A_j^{-1}P^T)^{-1} = E[Σ] − D`, and
`η_j = ‖A_{j+1}^{1/2}(P A_j^{-1}P^T)A_{j+1}^{1/2} − I‖`; the one-step lemma holds iff `η_j<1`.

The soft mode `v_soft` is the lowest (longest-wavelength) eigenvector of `A_j`, eigenvalue `ε := A_j(v_soft)`,
possibly `→0` (Savvidy softness). An isotropic low-pass `P` with symbol `p(k)` (`p(0)=1`, `p(π/a)→0`,
`S4` §1.1) satisfies `p(k_soft) = 1 − O((a k_soft)²)`, so a fraction `ℓ := 1−|p(k_soft)|² = O((a k_soft)²)` of
`v_soft` is projected into the fibre `ker P`. That leaked vector is still an approximate eigenvector of `A_j`
at eigenvalue `≈ ε`, so

```
    γ_j := λ_min( K ),   K = Π₀ A_j Π₀  (Π₀ = projector onto ker P)      ⇒   γ_j ≤ ε + O(ℓ·spread)  →  0 ,
    Schur loss  ⟨K^{-1}c,c⟩  ≳  (leaked weight)²·ε / γ_j   →  ∞ ,       η_j ≥ 1 .
```

In the `H1` 2-mode caricature (`H1-RG-DESCENT` D4, `a_ff = μ0−λβ²`) the fibre coordinate *is* the soft mode:
`G″ = a_cc − g²/a_ff`, which crosses `0` while `A` is still PD, and `η → ∞`. **The failure is the divergent
fibre loss `g²/γ_soft` with `γ_soft→0` — an artifact of putting the soft mode (or a piece of it) into the
fibre.** This is the fact the modified scheme exploits.

---

## 2 · The scheme: anisotropic two-step descent, soft mode retained

Two equivalent presentations; we use the first for the proof and the second for orientation.

**(a) Anisotropic blocks, long in the soft direction.** Partition the shell-`j` lattice into blocks that are
`L_∥` cells *along the soft-mode polarization/momentum direction* and `2` cells across it, with `L_∥ ≫ 2`
(long enough that `p(k_soft)=1` to machine precision on the soft mode). The block map `P_an = P_an(s_soft)`
then **retains the soft mode exactly** (`P_an P_an^T v_soft = v_soft`) and decimates only modes with
`|k| ≥ k_c := π/L_∥` across / `π/2` along. The fibre is `ker P_an = span{modes with |k| ≳ k_c}`, i.e. **only
hard modes**.

**(b) Two-step descent (equivalent).** (i) *Hard cascade:* iterate ordinary isotropic blocking **on the hard
sublattice** (project out the soft mode first) until the hard modes reach the cutoff — every fibre at every
sub-step is hard (`γ_hard > 0`), so each sub-step has `η^{hard}<1` with the standard constants
(`H1-RG-DESCENT` §3). (ii) *Soft coarsening:* once the hard modes are gapped and integrated, coarsen the soft
direction; the soft mode is the single surviving low mode, and its curvature is the transported effective
curvature. This is the standard adiabatic/Born–Oppenheimer split: integrate the gapped fast modes at fixed
soft-mode amplitude, then let the slow amplitude flow.

**(c) Separate soft-mode variable (equivalent, used for the floor).** Write `X = a·v_soft + ξ`, `ξ ⊥ v_soft`,
and integrate `ξ` exactly at fixed `a` (well-posed on the good set since the `ξ`-Hessian `K_ξξ` has gap
`γ_hard`). The soft-mode curvature is then the marginal curvature `G″_soft(a) = ∂_a²(−log ρ_v)(a)`, and the
`ξ`-fibre is gapped so `D≤0` (Prékopa + Brascamp–Lieb, `S3` §2.1).

Throughout, the **hard-fibre gap** is the block-scale Laplacian plus the Savvidy hard gap,
`γ_hard ≥ c_□^{blk} + 2g|B| = (π/L_blk)² + 2g|B| > 0` — *independent of the soft curvature `ε`*, because the
soft mode is no longer in `ker P_an`.

---

## 3 · The one-step inequality, with constants **[T]** (hypothesis `η^{an}<1`)

> **Lemma SM3.1 (anisotropic one-step descent).** Let `P := P_an` retain the soft mode and let
> `A_{j+1} := P A_j P^T` be the Galerkin coarse bare precision. With
> `η^{an}_j := ‖A_{j+1}^{1/2}(P A_j^{-1}P^T)A_{j+1}^{1/2} − I‖`, if `η^{an}_j < 1` then
> ```
>   G″_j  ⪰  (1+η^{an}_j)^{-1} A_{j+1}  =  c^{an}_j G″_{j+1}^{bare},     c^{an}_j = (1+η^{an}_j)^{-1} ≤ 1 ;
> ```
> additively `G″_j ⪰ c^{an}_j A_{j+1} − D_j^{+,hard}`, with `D_j^{hard} ≤ 0` (the fibre is hard ⇒ log-concave,
> Prékopa/BL, `S3` §2.1), so the defect is a *gain* on the good set.

*Proof.* Lemma H1.1 of `H1-RG-DESCENT` applied to `P=P_an`; `‖E‖<1 ⇒ (I+E)^{-1} ⪰ (1+η^{an})^{-1}I`. ∎

**The `η^{an}_j` budget — three pieces, the third now finite.**

| piece | source | magnitude | note |
|---|---|---|---|
| `η^{flow}` | `β_j ≠ β_{j+1}` (AF running) | `≍ g_j²` | removed by `Z_j` (`V1` §3) **[T]** |
| `η^{pb}` | passband leakage of the *hard* modes | `≍ (a_j k_c)²` | small; `p(π/a)→0` **[T]** |
| `η^{int}` | interaction / Savvidy, **hard** fibre | `ϑ²/γ_hard`, `ϑ` = soft–hard graph coupling | **finite** (was `g²/γ_soft→∞`) **[T]/[C]** |

The whole point: **`η^{int}` is `ϑ²/γ_hard` with a *hard* `γ_hard`, not `g²/γ_soft` with `γ_soft→0`.** Because
`γ_hard ≳ (π/L_blk)² + O(g|B|)` is bounded below independently of the soft curvature `ε`, the interaction piece
stays finite as `ε→0`. Two calibrations of the coupling appear below and agree to leading order: the **raw**
`ϑ=g` (H1-toy; loss `g²/γ_hard`, used in W1/W5) and the **graph-angle** `ϑ²≍g²A_j` (loss `g²A_j/γ_hard`, used
in W4); in both `η^{int}_soft ≈ g²/γ_hard` and `<1` for a hard fibre.

**Fixed-point sharpening on the soft mode itself.**

> **Proposition SM3.2 (soft mode is a fixed point of the blocking) [T].** Let `v` be retained by `P_an`
> (`P_an P_an^T v = v`, e.g. `v` a long-wavelength mode in the passband). Then for every PD `A_j`,
> ```
>   v^T (P_an A_j^{-1}P_an^T) v  =  v^T A_j^{-1} v ,      hence (v an eigenvector, A_j v = ε v)   G″_j(v) = ε = A_j(v).
> ```
> If the flow is turned off on the soft mode (`A_{j+1}(v) = ε`), then `η^{an}_soft = 0` and Lemma SM3.1 is an
> **equality**: `G″_j(v) = A_{j+1}(v)`.

*Proof.* `P_an P_an^T v = v ⇒ P_an^T v` is the lift, and `v^T P_an A_j^{-1}P_an^T v = (P_an^T v)^T A_j^{-1}(P_an^T v)`
`= v^T A_j^{-1} v`. For an eigenvector `A_j v = ε v`: `v^T A_j^{-1}v = |v|²/ε`, so `G″_j(v) = ε`. ∎

**Consequence `[NG]`.** On the retained soft mode the descent is an *identity* — it transports `ε`, it does not
manufacture positivity. And since the gapped hard fibre is a **positive** complement, its Schur elimination is
convention (A) (`S3` §2.4) and only *lowers* `G″`: `G″_j(v_soft) = ε − (Schur loss) ≤ ε = A_j(v_soft)`. Hence
**any** blocking — isotropic, anisotropic, or two-step — gives `G″_soft ≤ A_j(v_soft)`: no reblocking can push
the soft-mode curvature above its fine bare value. **The obstruction to a positive soft-mode floor is not `η`
(it is fixed to `0`/small); it is that the floor *is* the bare soft curvature, and its sign is the Savvidy
sign.**

---

## 4 · Does `η<1` hold on the soft mode at `β=2.5`? — YES, with the Haar floor **[T]/[C]/[V]**

Because `γ_hard ≳ 2.5` (lattice units) and `g²(β=2.5) ≍ 1–1.6`, the anisotropic interaction piece is
`η^{int}_soft ≍ g²/γ_hard ≲ 0.6`; with `η^{flow}` `Z`-removed and `η^{pb} = O((a/L_∥)²)≪1`, so
`η^{an}_soft < 1`. The witness makes the contrast with the isotropic break explicit, on the **same** soft
curvature (W1, W5):

```
  β      ε (no Haar)   η_ISO     η_AN(no Haar)   ε (+Haar N/12)   η_AN(+Haar)
 2.30    +0.1536        0.95        0.31            +0.3203           0.13
 2.40    +0.0784       22.06        0.85            +0.2451           0.17
 2.45    +0.0396         inf        10.00           +0.2063           0.21
 2.50    +0.0000         inf         inf            +0.1667           0.28
```

Readings: (i) at `β=2.4` the isotropic descent has `η_ISO = 22` (the `H1` break) while the anisotropic one
has `η_AN = 0.85 < 1` — *same* soft curvature, different blocking; the fix is real. (ii) On the Haar-free
toy, `ε→0` at `β=2.5` and even the anisotropic scheme fails (`η→inf`) — the pure-Savvidy soft mode with no
measure floor is genuinely gapless. (iii) **With the Haar floor** `Φ″ ⪰ N/12` (`Lemma B`), the bare soft
curvature is protected, `ε ≥ N/12 = 0.167 > g²/γ_hard`, and `η^{an}_soft = 0.28 < 1` all the way through
`β=2.5`. **The cure at `β=2.5` is carried by the measure (Haar) floor, not by the blocking** — foreshadowing
that the floor it proves is `N/12` (curvature), not the dynamical gap.

So the answer to the task's part (2) is: **`η<1` on the soft mode at `β=2.5` is achieved** (conditionally on
the hard-fibre gap, and at `β=2.5` on the Haar floor). The obstruction is *not* `η`; it is positivity
(Prop SM3.2 / `[NG]`), addressed in §5.

---

## 5 · The soft-mode floor: good-set curvature, and the essinf gap = `H1` **[T]/[C]/[O]**

By Lemma SM3.1 + Prop SM3.2 the soft mode obeys
```
   G″_soft  ⪰  (1+η^{an})^{-1} A_{j+1}(v_soft)  ≈  A_{j+1}(v_soft) − (hard Schur loss) ,     A_{j+1}(v_soft) = β_{j+1}Q(k_soft) + Φ″_soft + Sav_soft .
```
Two evaluations of the bare soft curvature `A_{j+1}(v_soft)`:

- **Good set (small field, `Sav≈0`).** `A_{j+1}(v_soft) ≈ β c_L + N/12` (`c_L = 4π²/L² = Q(k_soft)` to
  leading order; `Φ″_soft ⪰ N/12`). The transported floor is
  ```
     G″_soft  ⪰  (1+η^{an})^{-1}[β c_L + N/12]  −  g²(β c_L + N/12)/γ_hard     (W4, positive at β=2.5: ≈0.40, η≈0.28)
  ```
  Positive for `g²<γ_hard` (a hard fibre), across `β∈[0.5,5]` (W4). It is **volume-uniform via the Haar
  `N/12`** in the `L→∞` limit (`β c_L→0`). This is the positive content: a **β-dependent soft-mode curvature
  floor on the good set**.

- **Unconditional essinf (Savvidy tail).** On the Savvidy-indefinite tail `Sav → −1.1β` (`S4` §3.2;
  `PROVE-Gpp-positive` §3), so
  ```
     ε_essinf  =  β c_L + N/12 − 1.1β  <  0     for  β > β* ≈ 0.18     (W4; −2.35 at β=2.5)
  ```
  and the essinf floor `ε_essinf − Schur loss` is **negative at `β=2.5`**. The `[T]` fact: a positive
  *uniform* (essinf) floor would need `Sav_soft < β c_L + N/12 − δ_W ≈ 0.5` uniformly on the tail, which the
  proven `Sav` floor `−1.1β = −2.75` violates. **Proving the essinf floor positive at `β=2.5` is exactly
  `H1`** (the measure-carried positivity); the block descent cannot supply it (Prop SM3.2 / `[NG]`).

**What the floor proves, and what it does not.** The scheme proves the good-set floor, a **β-dependent,
volume-uniform-via-`N/12` soft-mode curvature bound** — the task's part (3). It does **not** prove the mass
gap, for three independent reasons (`V2` §7, `H1-RG-DESCENT` §4.2): (a) the volume-uniform part is the Haar
`N/12 = g_μ(0)`-type *measure curvature*, and `β=0` has `G″≥N/12` with `f₁²=0` — curvature ≠ atom; (b) the
unconditional (essinf) floor is negative at `β=2.5` (Savvidy), the good-set↔essinf gap being the measure's
concentration on the good set (`H1`); (c) the flat-chart soft mode is gauge-dependent (`main.tex` §6.6) — the
physical gauge-invariant mode's floor is not established here.

---

## 6 · Acid test — free field must not be seeded **[T]/[V]**

1. **Free massless field.** The soft mode is the lowest box mode; the sharp/anisotropic passband *retains* it
   (`p(k_soft)=1`), so `η^{an}=0` exactly (W3: `η≤1.1×10^{-16}`), and `G″_soft = A_{j+1}(v_soft) = βQ(k_soft) + m²`.
   At `k_soft→0` (`m=0`) this is `0` (W3: `L=8→1.46, L=64→2.4×10^{-2}, L=256→1.5×10^{-3} → 0`). **The scheme
   manufactures no gap on the free field** — it returns the physical soft curvature `βQ(k_soft)`, which `→0`.
   `[T]+[V]`
2. **The fixed-point identity (Prop SM3.2) IS the acid test.** `G″_j(v_soft) = A_j(v_soft)` for the retained
   soft mode (W2 to `7×10^{-15}`) is exactly the statement that the descent transports, never creates — so it
   *cannot* manufacture positivity on the free field, by construction. `[T]+[V]`
3. **Atom face.** `N/12` is volume-uniform but is Haar curvature (`V2` §7): the floor the scheme proves is
   `N/12 + βc_L − loss`, and on the free field (`Φ≡0`) this reduces to `βQ(k_soft)→0`. The descent is
   acid-test-clean precisely because it carries *curvature of the flat-chart soft mode*, never converting it
   into an atom. `[T]`

---

## 7 · Witness — `research/v3/w7/sm3_anisotropic_rg.py` (`py -3.12`, deterministic)

Quoted output (abridged):

```
[W1] eta on the SOFT mode: isotropic (soft fibre) vs anisotropic (hard fibre)   g=0.30, Lam=2.5
      eps      loss_ISO=g^2/eps  Gpp_ISO   eta_ISO   loss_AN=g^2/Lam  Gpp_AN   eta_AN
     1e-01       9.000e-01      +0.3000    3.00e+00      3.600e-02  +0.0640   5.62e-01
     1e-02       9.000e+00      -7.8000       inf        3.600e-02  -0.0260      inf
     1e-06       9.000e+04      -9.0e+04      inf        3.600e-02  -0.0360      inf
[W5] head sweep beta=2.3..2.5 : ISO vs ANISO eta, with / without the Haar floor
      beta   eps(no Haar)  eta_ISO  eta_AN(noH)  eps(+Haar)  eta_AN(+H)
      2.40   +0.0784       22.06      0.85        +0.2451      0.17
      2.50   +0.0000         inf       inf        +0.1667      0.28
[W2] fixed-point identity: P retains the soft mode  =>  Gpp(v_soft) = A(v_soft)
    Gpp(v_soft) = 6.0772209563 ; A(v_soft) = 6.0772209563 ; diff = 7.11e-15
[W3] free-field acid test: eta=0, Gpp = beta Q(k) + m^2 ; massless k->0 -> 0
    L=   8 m=0.0 Gpp=1.46e+00 ;  L= 256 m=0.0 Gpp=1.51e-03 ;  massive -> m^2=0.1667
[W4] beta-floor : GOOD-SET Gpp_good>0 across beta (at 2.5: 0.396, eta=0.282) ;
     UNCONDITIONAL essinf floor <0 for beta>beta*=0.176 (at 2.5: -2.354).
```

Readings: **W1** the isotropic fibre loss `g²/ε` diverges while the anisotropic loss `g²/Λ` is finite — the
scheme replaces the divergent loss by a bounded hard one; **W5** at `β=2.4` `η_ISO=22` (break) becomes
`η_AN=0.85` (pass), and `β=2.5` needs the Haar floor to keep `ε≥N/12`; **W2** the fixed-point identity is exact
(`7×10^{-15}`) — transport, not creation; **W3** free massless `→0` (no false gap); **W4** the good-set floor
is positive across `β` while the essinf floor is negative above `β*≈0.18` (head = `H1`). *Not a proof — a
witness for the scheme, its constants, and its obstruction.*

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| isotropic soft-mode leakage: `γ_j ≤ ε → 0`, loss `g²/γ_j → ∞`, `η≥1` | §1; `H1-RG-DESCENT` D4 (`η=22` at `β=2.4`) | **[T]/[V]** |
| anisotropic block / two-step / separate-soft-variable scheme, hard fibre | §2, `γ_hard ≥ (π/L_blk)²+2g|B|` | **[T]** |
| **Lemma SM3.1** `G″_j ⪰ (1+η^{an})^{-1}A_{j+1}`, `η^{an}=η^{flow}+η^{pb}+g²A_j/γ_hard` | §3, proved | **[T]**, `η^{an}<1` **[C]** |
| interaction piece **finite** (`g²A_j/γ_hard`) where isotropic had `g²/γ_soft→∞` | §3, W1 | **[T]/[V]** |
| **Prop SM3.2** retained soft mode ⇒ `G″_j(v)=A_j(v)` exactly (fixed point, `η_soft=0`) | §3, W2 (`7×10^{-15}`) | **[T]** |
| no blocking lifts `G″_soft` above `A_j(v_soft)` (convention A: hard-fibre Schur only lowers) | §3 | **[NG]** for "RG lifts soft-mode curvature" |
| **`η<1` on the soft mode at `β=2.5`** | §4, W1/W5, with Haar floor | **[T]/[C]/[V]** |
| good-set soft-mode floor `(1+η)^{-1}(βc_L+N/12) − g²(·)/γ_hard > 0` (≈0.40 at β=2.5) | §5, W4 | **[T]/[C]/[V]** |
| volume-uniform part = Haar `N/12`; `β=0` has `G″≥N/12`, `f₁²=0` — curvature ≠ atom | §5–§6; `V2` §7 | **[T]** |
| unconditional (essinf) soft-mode floor **negative** at `β=2.5` (`Sav→−1.1β`, `β*>0.18`) | §5, W4 | **[T]/[V]**, positivity = **[O]=H1** |
| modified descent gives the **mass gap** | **false**: floor is (good-set) curvature, `N/12` not atom; essinf negative; flat-chart mode | **[NG]** (of that claim) |
| acid test: free field `η=0`, `G″=βQ(k_soft)→0`; fixed point = no manufacture | §6, W2/W3 | **[T]+[V]** |
| **H1 (uniform soft-mode positivity at `β=2.5`)** | untouched; the blocking route is closed by Prop SM3.2 | **[O]** |
| **Clay** | unbuilt | **0%** |

---

## 9 · Sharpest line

**The Savvidy soft mode is a fixed point of every renormalization-group blocking. Retaining it in the coarse
field with anisotropic blocks (fibre gap `γ_hard ≥ (π/L_blk)²+2g|B|` instead of `γ_soft→0`) removes the
`η≥1` symptom — `η^{an}_soft ≈ g²/γ_hard ≲ 0.6`, and `<1` at `β=2.5` once the Haar floor keeps the bare soft
curvature `ε≥N/12` — so `G″_soft ⪰ (1+η^{an})^{-1}[βc_L+N/12] − (hard Schur loss) ≈ 0.40 > 0` on the good
set. But because `v^T(P A^{-1}P^T)v = v^T A^{-1}v` for the retained soft mode (exact, `7×10^{-15}`), the
descent only *transports* the soft curvature `ε`; the gapped fibre is a positive complement whose Schur
elimination can only *lower* it. Hence no blocking creates soft-mode positivity — the good-set floor is the
Haar `N/12` (curvature, not an atom; `β=0` has `f₁²=0`), and the unconditional essinf floor is negative at
`β=2.5` (Savvidy `−1.1β`). The modified descent therefore **controls the soft mode and proves a β-dependent
good-set curvature floor, not the gap**; and it shows from the blocking side that the residual is, exactly and
only, the Savvidy/measure sign — `H1`. **No QED. Clay 0%.**

---

*Refs.* `H1-RG-DESCENT.md` §1–§5 (Lemma H1.1, `η_j` budget, D4 break `η=22` at `β=2.4`, bound `(P A^{-1}P^T)^{-1}`);
`S3-UV-FESHBACH.md` §2.1–§2.4 (fibre gap ⇒ `D≤0` via Prékopa/BL; convention (A) positive complement hurts);
`S4-DECIMATION.md` §1.1 (block symbol `p(k)`), §3.2 (`E[H]≥−1.1β`), §4.2 (`κ*`); `V2-TRUE-BOUND.md` §7
(`N/12` is curvature not atom); `V1-UV-NORMALIZATION.md` §2–§5 (`c_L=4π²/L²`, `Z_j`, `δ_W`);
`docs/paper-v3/main.tex` §6 (`eq:onestep`, `eq:fisher`, §6.5 spectral no-go, §6.6 flat-chart gauge-dependence);
`PROOF-JACOBIAN-SAVES.md` §2 (Lemma B, `Φ″ ⪰ N/12`); Prékopa (1971); Brascamp–Lieb (1976); Schur complement /
Haynsworth. Witness `research/v3/w7/sm3_anisotropic_rg.py` (`py -3.12`, numpy 2.4.6; fixed-point identity
`7.1×10^{-15}`; all sections exit 0). No third-party sources downloaded; all external inputs standard and
cited from the corpus.
