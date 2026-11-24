# SLOT-BS-FIXPOINT — the Birman–Schwinger / eigenvalue fixpoint `κ=F(κ)` for the `0⁺⁺` channel, made **self-anchoring** by self-referencing its threshold: it is controller-free and acid-clean, but `F(0)=0` (the gapless fixed point survives), and the sole exclusion is the absent `k→0` `0⁺⁺` sign-definite instability

> 2026-11-24 · Fills the ONE open slot of `PROFILE.md` (P4: a `k→0` sign-definite self-consistency `R(𝔉,H)` with no gapless fixed point, uniform, acid-clean) by the **Birman–Schwinger / eigenvalue-fixpoint** route. Witness `research/unconditional/w7/slot_bs_fixpoint.py` (`py -3.12`, single job, toy model — topology only).
> Labels (`CLAUDE.md §2`): **[T]** true theorem / banked · **[C]** conditional · **[O]** exact missing step · **[NG]** no-go for a *named class* · **[V-here]** new this round (incl. verified in the witness) · **[J]** judgment · **[R]** external/inherited · **[RT]** rigorous toy (proved in the model, not transferred).
> **This note does NOT prove Clay and does NOT prove no BS fixpoint can exist (`∀class ≠ ∀possible`). No QED.**

---

## VERDICT (up front)

> **NOT CONSTRUCTED.** The Birman–Schwinger eigenfixpoint `κ=F(κ)` on the `0⁺⁺` two-gluon channel **is** controller-free (the threshold is self-referenced to `(2m₁)²=4κ`, i.e. the fixpoint's own edge — the "self-anchor"), non-positivity, nonabelian, and acid-clean. **But `F(0)=0`:** on the *stable* (subcritical) vacuum the massless `0⁺⁺` channel has no bound state below its own threshold, so `κ=0` is a fixed point and the gapless solution survives. **[NG]**
>
> The only exclusion is the **negative fixed point / sign-definite inversion** `F(0)<0`, equivalently the BS determinant at threshold `D_{\rm thr}<0`, equivalently `K_{BS}(\theta)>1` — i.e. a `0⁺⁺` **threshold bound state at zero binding**, which for massless constituents is a **`0⁺⁺` tachyon** = a `k→0` sign-definite instability of the vacuum. **This is absent at β≈2.5:** the coarse `k→0` mode is *measured positive* (`G″≈+0.4…+0.51`, W7), the rigorous bound is only `G″≳−0.08>-\text{const}` (Cauchy–Schwarz), and the negatives are 90–97% high-momentum (W6). The one genuine nonabelian negative mode (Savvidy, `λ_min=−g|B|`) is **`1⁻⁻`/adjoint**, and handling the channel (projecting to the `0⁺⁺` color-singlet) **removes the wrong sector but does not transfer its sign to `k→0`** — the singlet kernel is sign-indefinite (repulsive seagull + Oehme–Zimmermann ghost sign) and subcritical. So the "channel fix" trades *wrong sector* for *no sign*: net no exclusion. **[NG]/[O]**
>
> **The object is A3.** The BS determinant `det(1-K_{BS})` **is** the Fredholm determinant `D_O(z)=det(1-zT_O)` (`META-CROSSDOMAIN-PROTOTYPE §1.4.3`); its zero is the atom, its derivative the residue. Locating that zero uniformly as `a→0` is the marginal zero-finding = transport = `I2=H1`. The self-referenced edge makes the map non-contractive across β and the `α_s~O(1)` ladder is IR-uncontrolled. **Non-uniform. [O]**
>
> **Exact failing step = the gapless fixed point persisting; the escape = the k→0 sign-definite `0⁺⁺` instability (the same `H1`); the object = `A3` (READ, not FORCE).** Controller-free but **circular = A3**. Clay 0%.

---

## §0 · Setting: operator, Hilbert space, nonabelian input

**Transfer matrix.** `T=e^{-aH}` is self-adjoint, positive on the physical Hilbert space `H_phys=ker Q_BRST⊂L²(𝒜)` (`Gauss`-invariant, positive-norm). `HΩ=0`, unique vacuum; on `(vac)^⊥` bottom `E₁=m₁` (`0⁺⁺`), next `E₂`. The `0⁺⁺` interpolator `O=¼ tr F²`. Källén–Lehmann measure `ν≥0` on `s=m²`, `C(τ)=⟨O(τ)O(0)⟩=∫e^{-√sτ}dν`, atom `f₁²=ν({m₁²})`, support floor `m₁=inf supp ν` (`H1-IS-IT-THE-GAP`, `[T]`).

**The `0⁺⁺` two-gluon subspace.** `P_{0⁺⁺}` = projection onto `J^{PC}=0⁺⁺`, **color-singlet** sector of `adjoint⊗adjoint`, total momentum `P→0` (rest frame). On this subspace the interacting two-body operator is
```
     𝓗 = 𝓗₀ + 𝓥 ,        𝓗₀ = the free two-gluon kinetic operator ( S⁻¹⊗1 + 1⊗S⁻¹ ) ,
     𝓥  = the two-gluon-irreducible BS kernel (one-gluon exchange + seagull + self-energy/ghost).
```
`S` is the single-gluon propagator; `𝓖₀(s)=(𝓗₀-s)^{-1}` the free two-body resolvent. The homogeneous BS equation is `(1-𝓖₀(s)𝓥)Φ=0`.

**Nonabelian input (two exact places `[T]`).**
1. **Color factor.** For two adjoints coupling to the singlet `R=1` (`C_R=0`), one-gluon exchange has eigenvalue `½(C_R-2C_A)=½(0-2N)=-C_A=-N` (`N≥2`) — *sign-definite attractive*, stronger than `q\bar q`'s `-4/3` (`BS-KERNEL-SIGN §1`). For pure `U(1)` the adjoint is trivial (`C_A=0`, no charged vector) ⇒ factor `0` ⇒ **IL-ii satisfied structurally** (`KEY-DIRECT-DERIVATION §1`).
2. **The Savvidy negative mode.** `Hess=-D²+2\,ad F_B`, transverse-magnetic `λ_min=-g|B|<0` (`PROOF-SAVVIDY-EXACT-SIGN`) — an **adjoint/`1⁻⁻`** object, **not** the `0⁺⁺` singlet. It enters the `0⁺⁺` kernel only at second order (`1⁻⁻⊗1⁻⁻⊃0⁺⁺`) through high-momentum loop modes (`BS-KERNEL-SIGN §5`). **Channel handled**: projecting to the `0⁺⁺` singlet removes it from the external `1`-set.

---

## §1 · The functional `F` and the fixpoint equation (Task 1)

**Birman–Schwinger eigenvalue.** For general self-adjoint `𝓥` (sign-indefinite — we do **not** assume `V≥0`), the bound-state content is read by the BS operator / determinant
```
   K_{BS}(s) := 𝓥^{1/2}(𝓗₀-s)^{-1}𝓥^{1/2}  (form sense on the 0⁺⁺ channel),
   D(s) := det(1 - 𝓖₀(s)𝓥) = det(1 - K_{BS}(s))   (Fredholm det, trace-class kernel).   (1)
```
An isolated `0⁺⁺` state sits at the largest `s=κ` below the continuum with `D(κ)=0` (eigenvalue `1` of `K_{BS}(κ)`). The threshold is `θ` (`θ=0` massless; `θ=(2m_g)²` with constituent mass; physically the multi-glueball/`2m₁` scale, `BS-THRESHOLD §1`).

**The self-anchoring fixpoint.** Let `κ:=m₁²`. Define the **self-consistent** spectral functional *with every input dressed at scale `κ` and the continuum edge self-referenced to the fixpoint's own scale*:
```
   F(κ) := the location s* of the lowest zero of D_κ(s)  in the 0⁺⁺ channel,
           where S=S_κ (self-energy dressed at κ), 𝓥=𝓥_κ, and the edge θ = θ(κ) = 4κ
           (the two-glueball threshold (2m₁)², Single-Scale Covariance).              (2)
```
The **fixpoint equation** is
```
        κ = F(κ),      κ = m₁² ≥ 0.                                                (★)
```
Two structural identifications, both exact:
- **`θ(κ)=4κ` is the self-anchor.** In confinement there is no asymptotic gluon; the only physical threshold above the lightest `0⁺⁺` is the multi-glueball scale `θ_c≥(2m₁)²=4κ` (`BS-THRESHOLD §1`, `[T]`). Setting the edge to `4κ` makes `(★)` a genuine self-consistency rather than a one-shot criterion with an external edge. **This *is* Single-Scale Covariance `E₂=ρ_r m₁`** (`G2D-ISOLATION-I2 §2`, `SSC`, `[T]`): the isolation margin `θ-m₁²=(4-1)κ` is a fixed multiple of the gap. So the "self-anchor" is the statement that the edge is the unknown — **it is not an independent anchor; it is the gap appearing on the right side of `(★)`.** **[V-here]** (the BS-coordinate reading of SSC).
- **`D(s)` is the Fredholm determinant `D_O(z)`** under `z=e^{-a√s}`, `T_O=T|_{0⁺⁺}` (`META-CROSSDOMAIN-PROTOTYPE §1.4.3`; `docs/new-tools-registry.md §A3`). `(★)` in the `z` variable is the Fredholm zero `z₁=e^{-am₁}`.

**Witness `[RT]`.** The separable–kernel toy (`slot_bs_fixpoint.py`) realizes `(★)`'s topology exactly: `D_{\rm thr}=1-g/g_c(m_g)`, `g_c(m_g)=2π²m_g/\intχ²`, bound state `m_B=2m_g-B`; see §2.

---

## §2 · Gapless-fixed-point test (Task 2)

### 2.1 `F(0)=0` — the gapless solution survives. `[NG]`

At `κ=0` the input is the free massless field (`Σ_0=0`), edge `θ(0)=0`; the `0⁺⁺` continuum is the free `trF²` cut from `0` (`BS-THRESHOLD §1`, T1). A stable vacuum has **no `0⁺⁺` state below its own (zero) threshold** — "bound below `0`" is a negative-energy state = a tachyon, which is absent. Hence `K_{BS}(0)` is **subcritical** (`spec<1`), `D_0(s)>0` for `s≥0`, no zero; the lightest object is the continuum edge, so
```
        F(0) = 0 .                                                                  (3)
```
`κ=0` is a fixed point of `(★)`. **The intrinsic relation admits the gapless solution** — the same obstruction `PROFILE-FACET-ANCHOR Theorem SA` proves for every structural relation (`[T-schema]`), here in the BS coordinate.

**β≈2.5 signal `[T]+[V-here]`.** That `K_{BS}(0)` is subcritical is not merely unproven — the *measured* massless vacuum is **locally stable** at the coarse mode: `G″(k→0)≈+0.4…+0.51` (`W7-Gpp0`), the rigorous bound only `G″≳W_geo-2√(VarG\,VarΦ′)=-0.08` (Cauchy–Schwarz, sign not certified; `B4-VARIANCE-MAG`, `ABC-FINAL-LEDGER`), and the negative Hessian directions are 90–97% high-momentum (`W6-BE-at-2.5`). A stable coarse mode ⇒ `K_{BS}(0)<1` ⇒ `D(0)>0` ⇒ `F(0)=0`. **The measured sign points *away* from the exclusion.** `[J]`

### 2.2 The only exclusion is the sign-definite inversion `F(0)<0` — and it is the absent `k→0` `0⁺⁺` instability. `[NG]/[O]`

By the task's own criterion, since `F(0)=0` one must exhibit a **negative fixed point / sign-definite inversion**. In BS form, excluding `κ=0` requires `D_{\rm thr}<0`, i.e.
```
   D_{\rm thr} = det(1-K_{BS}(θ)) < 0   ⟺   K_{BS}(θ) has an eigenvalue >1
   ⟺   a 0⁺⁺ bound state at ZERO binding (at the massless threshold) = a 0⁺⁺ TACHYON.
```
This is **verbatim the `k→0` sign-definite instability** the corpus proves absent:
- **Gap equation:** `F(0)<0` ⟺ `F′(0)>1` pitchfork ⟺ a `k→0` tachyon (`MS-GAP-EQUATION §2`, `[T]`). The BS fixpoint and the gap equation are **the same fixpoint** `κ=m₁²`, `F` increasing iff the coarse mode is unstable. `[V-here]`
- **Banach:** `T(0)=0` and the only contractively-certified fixed point near `0` is the gapless one; a nonzero branch needs `DT` crossing `1` = the order parameter (`BANACH-FIXED-POINT §5`, `[T]`).

At β≈2.5 the instability is **not in the coarse `0⁺⁺` mode** (coarse mode positive/measured, §2.1). The one genuine negativity is the Savvidy `1⁻⁻` `λ_min=-g|B|`, which is (i) **wrong sector** (adjoint, `1⁻⁻`) and (ii) **high-momentum** (`D-BBS-LOWRANK-KATO §(c)`: the band bottoms at `k_∥=0` but its 3-4-momentum lives at `k_∥≤√(gB)` and it is 90–97% orthogonal to the coarse slice).

### 2.3 The "channel fix" trades *wrong sector* for *no sign* — net no exclusion. `[NG]/[V-here]`

The task's instruction — *"the negative/Savvidy mode is adjoint, wrong channel; bind the `0⁺⁺` as a 2-gluon below the cut"* — is **implementable** and we implement it: project `𝓥` onto the `0⁺⁺` color singlet. Consequences, all `[T]`/`[V-here]`:
- **Acid gain:** the projected kernel no longer certifies the Savvidy tachyon (it is `1⁻⁻`), so the false positive `‖K‖≥1`-detecting-an-instability is gone. The channel is now the *right* one (singlet, `J^{PC}=0⁺⁺`).
- **Sign loss:** with the wrong-sector negativity removed, the *remaining* `0⁺⁺` kernel is **sign-indefinite and subcritical**: colour is attractive (`-N`) but the full `J=0` projection is indefinite — repulsive 4-gluon seagull + Oehme–Zimmermann ghost/self-energy sign change (`BS-KERNEL-SIGN §2`, `[T]`). And even a *net-attractive* kernel **cannot lift a massless-threshold state off `0`** (attraction binds *below* a threshold it cannot raise; `BS-KERNEL-SIGN §3`, `[T]`). So `D_{\rm thr}>0`, `F(0)=0`.
- **The strongest self-anchoring attempt fails.** Using the `1⁻⁻` negativity as a *seed* to bind the singlet — Feshbach/Schur elimination of the tachyon pair from two-gluon space — does `C≺0 ⇒ -BC^{-1}Bᵀ⪰0`, i.e. it **only adds positivity to a kernel already indefinite**, so `Σ≻0⇔K≻0 = H1` (`FORGE-FESHBACH-ANGLE`, `BS-KERNEL-SIGN §5`, `[T]`). The second-order feed lands on the high-momentum modes the coarse mode does not couple to (`PROOF-IR-H1 §4`).

**§2 verdict.** `F(0)=0` holds; the sign-definite inversion `F(0)<0` requires a `k→0` `0⁺⁺` tachyon that is **absent** at the coarse mode at β≈2.5; channel-cleaning is correct but removes the sign along with the wrong sector. **The gapless fixed point is not excluded.** `[NG]`

**Witness `[RT]`** (`slot_bs_fixpoint.py`): `D_{\rm thr}=1-g/g_c` is `<0` **iff** `g>g_c` (supercritical/tachyon) — the sign-definite inversion is exactly supercriticality; the subcritical stable vacuum has no bound state below threshold (`F=θ`, gap `0`); a positive gap appears only for `g>g_c` (tachyon onset) or with a lifted threshold `m_g>0` (controller).

---

## §3 · Uniform-attainment test: is `(★)` attained uniformly across β≈2.5? (Task 3)

**No — and for two independent, exact reasons.** `[O]`

1. **Zero-finding is the Fredholm zero = transport = `I2`.** Even granting a positive zero `κ*>0`, its uniform location is the zero `z₁=e^{-am₁}` of `D_O(z)` (`§1`). To bracket `z₁` and bound `D_O′(z₁)` (the residue `f₁²∝N_O/D_O′`) **uniformly as `a→0`** is the marginal BPHZ zero-finding wall (`META-CROSSDOMAIN-PROTOTYPE §1.4.3`; `LOC-DETERMINANT §5`; `new-tools-registry §A3`). The consistency condition is exactly the isolation `S₂=o(λ₁)=I2`. This is **`READ ≠ FORCE`** (`META-GAPFORCING-IDENTITY §3.5`): the Fredholm determinant *reads* the atom when given the spectrum but cannot *force* it uniformly.
2. **The self-referenced edge is non-contractive; the ladder is IR-uncontrolled.** Setting `θ(κ)=4κ` puts the unknown on both sides of `(★)`; the map is not a contraction across β (no fixed spatial reference — the "anchor" moves with the answer; `SSC`, `G2D-ISOLATION-I2 §2`). At the binding scale `g²N/4π=O(1)` (`BS-KERNEL-SIGN §4`): cross-ladder/vertex corrections are same order and the truncation breaks the STI (scheme/gauge-dependent). **No small parameter ⇒ IL-iii fails.**

So the fixpoint is **non-uniform** even where a positive branch might numerically exist; uniformizing it is `I2/H1`.

---

## §4 · Acid test + the three iron laws (Task 4)

**Acid test (free massless must NOT be certified). `[T]+[RT]`** The mechanism passes it **by construction and non-vacuously in form**: `g=0 ⇒ 𝓥_eff=0 ⇒ K_{BS}≡0 ⇒ D≡1>0`, no zero, `F(0)=0`, `f₁²=0` — the free field is **not** certified (`F1`). The three adversaries:
| adversary | `𝔉_{BS}` | reason |
|---|---|---|
| `μ₀` free massless | `0` | `g=0⇒K_{BS}=0`, no bound state |
| `μ_{U(1)}` pure | `0` | adjoint trivial ⇒ color factor `0`, no Savvidy charged vector (F2) |
| `μ_m=:φ²:` free massive (gapped but **atom-free**, the *discriminating* acid) | `0` | no extra gauge binding below its own `2m` threshold ⇒ no atom |
**Acid PASSES on all three** — but this is exactly the trap `BANACH-FIXED-POINT §5` names: the fixpoint is acid-clean **because it certifies nothing positive** (it returns the gapless fixed point).

**The three iron laws (`KEY-DIRECT-DERIVATION §0`). `[T]`**
- **(IL-i) non-positivity.** The BS criterion is *not* a positivity/convexity/LSI inequality — it is driven by the sign-indefinite kernel and negative part. **PASS in form.** *But* the output collapses to `F(0)=0` on the free field, so it is **vacuous in output** — it never produces a strict bound. `[NG]`
- **(IL-ii) uses the nonabelian interaction.** PASS structurally: colour factor `-N` (`½(C_R-2C_A)`, `=0` for `U(1)`); the Savvidy negativity exists only for `N≥2` (`KEY-DIRECT-DERIVATION §1`, `[T]`). `[T]`
- **(IL-iii) no small parameter.** **FAIL:** the positive branch (if any) needs either the absent `k→0` tachyon or an `O(1)`-coupling ladder that is IR-uncontrolled at `α_s~O(1)`, β≈2.5 (`BS-KERNEL-SIGN §4`). `[NG]`

So the fixpoint satisfies IL-i (in form) and IL-ii, fails IL-iii, and is **acid-clean vacuously** — it passes the letter and fails the spirit (it never leaves the gapless branch).

---

## §5 · Verdict: constructed OR the exact failing step (Task 5)

**NOT CONSTRUCTED.** The BS eigenfixpoint `(★)`, made self-anchoring by the self-referenced edge `θ(κ)=4κ`, is:
- **(a) no gapless fixed point?** **NO.** `F(0)=0` (§2.1). The gapless solution survives — the *primary* failing step. The only exclusion is the sign-definite inversion `F(0)<0` ⟺ a `k→0` `0⁺⁺` tachyon, **absent** at the coarse mode (§2.2); channel-cleaning removes the wrong sector but not its absence (§2.3).
- **(b) uniform attainment?** **NO** (§3). Zero-finding = Fredholm = `I2/H1`; edge non-contractive; ladder `O(1)`.
- **(c) non-positivity, nonabelian, acid-clean?** Non-positivity PASS-in-form (IL-i), nonabelian PASS (IL-ii), acid PASS — but **vacuously** (returns the gapless branch).

**The exact failing step is one object in three coordinates (§2.2):** the `k→0` sign-definite `0⁺⁺` instability — `F(0)<0` (BS) ⟺ `F′(0)>1` pitchfork (gap equation, `MS-GAP-EQUATION §2`) ⟺ `DT` crossing `1` (Banach, `BANACH-FIXED-POINT §5`) ⟺ a `0⁺⁺` tachyon = the single-gluon IR mass gap `m_g>0` in the physical sector (`BS-KERNEL-SIGN §3`). This is `H1` in the BS alphabet.

**Verdict label.** `[NG]` for the **Birman–Schwinger eigenfixpoint class** (it is controller-free, genuine, nonabelian, acid-clean, and *defeated by the gapless fixed point*, exactly as `PROFILE-FACET-ANCHOR Theorem SA` predicts), and `[O]` for the exact missing step (the `k→0` `0⁺⁺` sign-definite exclusion = `H1`). The object is **`A3`** (BS det = Fredholm `D_O(z)`, the unique floor-sensitive survivor, **READ not FORCE**). So the route is **circular = A3**: the "self-anchor" `θ=4κ` is `SSC/I2` = the gap on the RHS; excluding `0` needs the gap; reading the positive branch uniformly is transport.

**What would flip it (`[T]`-criterion).** A proof of **`D_{\rm thr}<0`** at β≈2.5 by non-edge data — i.e. a `0⁺⁺` threshold bound state at zero binding, equivalently a proven `k→0` `0⁺⁺` sign-definite instability (not the `1⁻⁻` Savvidy tachyon, not high-momentum). By `PROOF-IR-H1 §4` a fast-mode (high-momentum) negativity cannot gap the slow mode; by `KEY-DIRECT-DERIVATION §1` the only signed nonabelian structure is the `1⁻¹`/adjoint instability. **None is known; exhibiting it is `H1` itself.** No QED. **Clay 0%.**

---

## §6 · Relation to the corpus (honest placement)

- **Confirms `PROFILE.md` P4 and the "representational freedom is illusory" claim.** The BS fixpoint is **the same fixpoint** as the gap equation (`MS-GAP-EQUATION`), the Banach contraction (`BANACH-FIXED-POINT`), and the Fredholm zero (`LOC-DETERMINANT`): all have `F(0)=0` and all need the identical `k→0` sign-definite exclusion. **New here:** the BS-coordinate *identification* of the self-anchor `θ(κ)=4κ` with **SSC/I2** — the "controller-free" threshold is the gap itself, so the BS route is *A3-circular*, not merely "controller-needing" (`[V-here]`).
- **Sharpens the prior "BS needs a controller" finding (`META-CROSSDOMAIN-PROTOTYPE §1.4`, `BS-THRESHOLD`, `BS-KERNEL-SIGN`, `BS-VARIATIONAL-LOWER`) into a self-anchoring test:** making the threshold self-referenced (the *only* available self-anchor) replaces the external controller `(2m_g)²` by `(2m₁)²` — which does **not** supply an independent anchor, it reproduces `SSC` and inherits the gapless fixed point. Channel-handling (the task's promotion) is **correct and necessary** but **insufficient**: it removes the wrong-sector sign without creating the `k→0` sign.
- **Boundary.** `∀class ≠ ∀possible`: this is `[NG]` for the BS eigenfixpoint with a self-referenced edge (the natural self-anchoring). A BS variant with an independent, non-edge, `k→0` sign-definite seed would escape — but supplying that seed *is* `H1`.

---

## References

**Internal (absolute paths; labels carried):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall=Savvidy concavity / key=`G″(0)`=H1); `docs/BRIEF-FOR-MODEL.md` (§1 atom, §2 five bars, §4 transport≡confinement); `docs/new-tools-registry.md §A3` (Fredholm determinant).
- `research/unconditional/PROFILE.md` + `PROFILE-FACET-ANCHOR.md` (the one open slot = `k→0` sign-definite self-consistency; Theorem SA, SSC-adjacent; SA-1..SA-5); `SPECIES-EXISTENCE.md` (`F3∧F4=∅`, equivalence/circularity, READ≠FORCE); `META-GAPFORCING-IDENTITY.md` (§3.1 identity-forcing no-go; §3.5 READ≠FORCE).
- `META-CROSSDOMAIN-PROTOTYPE.md §1` (BS exact structure; §1.4 import: topology inverted, detects Savvidy tachyon, residue=A3); `BS-THRESHOLD.md` (T1/T2 thresholds, `θ_c≥(2m₁)²`); `BS-KERNEL-SIGN.md` (colour `-N` attractive, kernel indefinite, **threshold crux §3**, Feshbach §5); `BS-VARIATIONAL-LOWER.md` (Temple needs kinematic edge `ρ=4m_g²`; `0⁺⁺*` is the next level, not the continuum); `BS-BLIND.md` (Temple/BS no lower bound); `D-BBS-LOWRANK-KATO.md` (band soft at `k=0`, `ε=O(1)`, Kato circular at `k→0`).
- `MS-GAP-EQUATION.md` (§2 `κ=0` always, `F′(0)>1` pitchfork; tachyon UV); `BANACH-FIXED-POINT.md` (§5 `T(0)=0`, bifurcation=order parameter); `G2D-ISOLATION-I2.md §2` (SSC `E₂=ρ_r m₁`); `H1-IS-IT-THE-GAP.md` (squeeze, `G″=1/χ`, `:φ²:` gapped atom-free); `LOC-DETERMINANT.md` (`D_O` zeros, `S₂=o(λ₁)=I2`); `PROOF-IR-H1.md §4` (fast-mode cannot gap slow-mode).
- `KEY-DIRECT-DERIVATION.md §0` (**three iron laws**), §1 (Savvidy nonabelian-only, `λ_min=-g|B|`); `KEY-DIRECT-OPTIMIZE.md §A.4`; `B4-VARIANCE-MAG.md` / `ABC-FINAL-LEDGER.md` (CS `G″≳-0.08`, measured `+0.4`); `W6-BE-at-2.5.md` (negatives 90–97% high-momentum); `W7-Gpp0` (`G″≈+0.51`); `PROOF-SAVVIDY-EXACT-SIGN.md` (`λ_{n,s}` ladder).
- Witness: `research/unconditional/w7/slot_bs_fixpoint.py` (this work).

**External `[R]` (standard; not re-verified this session):** Birman, *Mat. Sb.* **55(97)** (1961) 125; Schwinger, *PNAS* **47** (1961) 122; Simon, *Quantum Mechanics for Hamiltonians Defined as Quadratic Forms* (1971); Reed–Simon IV §XIII; Salpeter–Bethe, *Phys. Rev.* **84** (1951) 1232; Källén (1952)/Lehmann (1954); Savvidy, *Phys. Lett. B* **71** (1977) 133; Nielsen–Olesen (1978); Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Kato, *Perturbation Theory for Linear Operators* (1966); Prékopa (1971); Brascamp–Lieb, *Adv. Math.* **22** (1976).

---

*Status: **SLOT-BS-FIXPOINT.** Formulates the Birman–Schwinger eigenfixpoint `κ=F(κ)` of the `0⁺⁺` two-gluon channel with a **self-referenced threshold** `θ(κ)=4κ` (the SSC self-anchor). It is controller-free, nonabelian, acid-clean — but `F(0)=0` (gapless fixed point survives; stable subcritical coarse mode, measured `G″≈+0.4`). The sole exclusion `F(0)<0` (sign-definite inversion) is a `k→0` `0⁺⁺` tachyon, absent at β≈2.5; channel-cleaning removes the wrong-sector Savvidy `1⁻⁻` but not its absence. Non-uniform (Fredholm zero-finding = `I2`; edge non-contractive; `O(1)` ladder). Object = BS det = Fredholm `D_O(z)` = `A3` = READ not FORCE. **Verdict: `[NG]` for the BS-eigenfixpoint class, `[O]` exact step = `k→0` `0⁺⁺` sign-definite instability = `H1`; circular = A3.** `∀class≠∀possible`. Clay NOT proved; 0%. No QED.*
