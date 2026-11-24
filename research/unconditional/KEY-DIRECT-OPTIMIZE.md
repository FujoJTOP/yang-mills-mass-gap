# KEY-DIRECT-OPTIMIZE — optimizing the Savvidy signed minorant (KEY-DIRECT attack I): can a positive-at-`B=0` functional `Φ_total` cross the floor, or does it die on the same wall from a sharper angle?

> 2026-11-24 · **Focused optimization** of the best attack in `KEY-DIRECT-DERIVATION.md` (Attack I: the Savvidy signed minorant). The task: build `Φ_total(B) = [signed Savvidy part] + [W_geo / dynamical part]` that is **strictly positive at `B=0`**, attempt `Φ_total ≤ m₁`, and decide whether **any** optimization gives a genuine `0⁺⁺` atom floor `f₁² ≥ c>0` (or `m₁² ≥ λ>0`) at the crossover `β≈2.5`, uniform `(L,a)`. Witness: `research/unconditional/w7/key_direct_optimize.py` (`py -3.12`, single small job; exact spectra + 1D/2D quadrature + tiny toys; no lattice MC).
>
> Labels (charter §2): **[T]** true theorem · **[C]** conditional · **[NG]** no-go for a named class · **[O]** exact missing step · **[V-here]** computed this round · **[R]** external. **No QED on anything unproven. Clay 0% is the expected honest answer; the point is to try genuinely and then pin the exact structural death.**

---

## VERDICT (up front)

> ### **No optimization yields a genuine `0⁺⁺` atom floor. The literal task — a functional strictly positive at `B=0` — is achieved (`Φ_total(0)=W_geo>0`, [V-here]), but it is a `G″` (curvature = `1/χ`, inverse susceptibility), which by the exact squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` is an UPPER bound on `m₁²`, not a floor; and it is blind to a thin/absent edge atom (witness `dν=√s ds`: `χ<∞`, `G″>0`, yet `inf supp ν=0`, verified [V-here]). So the object gives exactly H1 (`G″>0`), not `H*`.**
>
> **Separately, the channel is wrong (§B): the Savvidy negative mode is a colored (adjoint), transverse-magnetic, `s_z=−1` single gluon — a 1-form (`J≥1`) object with **zero** `0⁺⁺` projection at one-particle level (color singlet⊕adjoint = 0 by Schur; helicity ±1 excludes `J=0`). The `0⁺⁺` atom `f₁²` is a color-singlet scalar residue; the negative mode enters `0⁺⁺` only as a two-gluon CUT above `2m_eff=2√(gB)>0`, and **does not exist at all at the physical vacuum `⟨B⟩=0`** (operator `−D²+2adF_B = −D²⪰0`). So even a hypothetical `Φ_total ≤ m₁` would bound the wrong channel.**
>
> **Precise structural death = a CONJUNCTION: (1) H1-not-`H*` (the positive piece `W_geo` is a curvature — upper-bounds `m₁²`, blind to the atom); (2) wrong-channel (the signed piece is a colored 1-form magnetic mode, not `0⁺⁺`); (3) width-not-floor (the Savvidy scale `√(2gB)→0` precisely where the gap is largest); (4) no small parameter (crossover at `α_s~O(1)`, no dynamical lever). §C shows the H1→`H*` identity upgrade fails (READ ≠ FORCE). This is the Clay wall, seen from the sharpest angle yet — transparent, measured, unmoved. Clay 0%. No QED.**

---

## §0 Setup, the construct, and what "optimize" can and cannot change

**The construct (from `KEY-DIRECT-DERIVATION.md` §1).** The gluon fluctuation Hessian in a constant background `B` is `−D² + 2 ad F_B`. Its only genuinely *signed* (indefinite) structure — unavailable to U(1) — is the spin–field coupling `2 ad F_B`, which yields a **negative mode** `E² = −gB` whose sign is fixed by `C₂(adj)=N>0`. Take `Φ_Sav(B) := √(2gB)`. **[T]** It passes IL-(i) (signed, not a positivity), IL-(ii) (U(1): 0 negative modes, SU(2): ≥1), and the acid test (`→0` at `B=0`, free/U(1)/massless). **Its death is `Φ_Sav ≤ m₁`:** `√(2gB)` is monotone increasing, `→0` as `B→0`; the physical vacuum has `⟨B⟩=0` yet `m₁>0`, so `Φ_Sav(0)=0≪m₁`; it is the curvature of an unstable saddle (a width), not a spectral floor.

**The optimization instruction.** Combine the signed Savvidy structure with the **proved-positive Haar curvature** `W_geo≥0` (Prékopa; `PROOF-JACOBIAN-SAVES.md`) and/or the dynamical curvature, to make
```
Φ_total(B) := [ signed Savvidy part √(2gB) ]  +  [ W_geo (≥0) ] ,
```
**strictly positive at `B=0`**, and try `Φ_total ≤ m₁`. The named missing inequality is `W_geo > D − Δ_dyn − X` (measured margin `D/S≈2%`).

**What "optimize" can change, and what it cannot.** An optimization can make the *free-field / acid / small-`B*`* properties cleaner (this it does, §A.4). It **cannot change the TYPE** of the object: any combination of a *background-field scale* (`√(2gB)`) and *curvatures* (`W_geo`, `Δ_dyn`) is still a `G″`/scale object. The floor `m₁²≥λ` is a statement about `inf supp ν`, which requires the **edge residue** `f₁²` (`H*`), a different functional the optimization never supplies. That type fact is the whole verdict; §A–§C make it exact.

---

## §A The positive-at-`B=0` functional `Φ_total` — derivation, numeric, and H1-vs-floor

### A.1 It IS strictly positive at `B=0` (the literal ask is achieved) [V-here]

The exact split of the coarse-marginal curvature (`PROOF-JACOBIAN-SAVES.md` §1) is
```
G″(0) = Δ_dyn + W_geo + X ,   Δ_dyn = E[H]−Var(G) ,   W_geo = E[Φ″]−Var(Φ′) ,   X = −2Cov(G,Φ′) ,
```
with `Φ = −log J` the (nonabelian) Haar Jacobian. Two facts fix the pieces:
- **`W_geo ≥ 0`, β-independent, B-independent** (Prékopa: `Φ` convex ⇒ `e^{−Φ}dX` log-concave ⇒ the pushed linear marginal is log-concave ⇒ `V_eff` convex). At `β=0` it is the whole curvature, `W_geo = G″ ≥ N/12` (`[T]` PROOF-JACOBIAN-SAVES; SU2 `1/6`, SU3 `1/4`). At `β=2.5` measured `W_full = W_geo+X ≈ +0.22–0.31` (`[R]` W7-RERUN, PROOF-JACOBIAN-SAVES §3).
- **`W_geo` is nonabelian-only.** For U(1) the link measure is flat (`J=1`, `Φ=0`), so `W_geo=0`. For SU(N) the Vandermonde makes `Φ≠0` and `W_geo>0`. (§A.5 verifies this on a single link.)

Hence
```
Φ_total(B) = √(2gB) + W_geo  ⟹  Φ_total(0) = 0 + W_geo = W_geo > 0   (SU(N)) ,  = 0 (U(1)) .
```
**This satisfies the task's literal request: a functional built from the signed Savvidy negativity plus the proved-positive Haar curvature that is strictly positive at `B=0`.** **[V-here]** Numerically §A.5 gives `W_geo^{SU2}≈0.205` on a single link (convention-dependent; sign & U(1)-vanishing are the invariants), consistent with the corpus's `~0.22–0.31` at `β=2.5`.

### A.2 The decisive point — it is a CURVATURE, so it gives H1, NOT a floor [T]

`Φ_total` is a sum of a background scale and curvatures; evaluated at `B=0` it is a `G″`. But `G″` is, exactly, the **inverse zero-momentum susceptibility** (`H1-IS-IT-THE-GAP.md`):
```
G″(k=0) = 1/χ = 1/∫ s⁻¹ dν ,    χ = ∫ s⁻¹ dν = C̃(0).
```
The exact squeeze on the KL measure `ν` (`ν({m₁²})=f₁²`, `∫dν=M₀`) is
```
f₁²/χ  ≤  m₁²  ≤  M₀/χ  =  M₀ · G″ .            (真定理, H1-IS-IT-THE-GAP §2.2)
```
So `G″ ≥ m₁²/M₀` — **`G″` is a residue-weighted harmonic mean, which is ≥ min**, i.e. `Φ_total` is an **UPPER bound on `m₁²`**, the *opposite* of a floor. In particular `Φ_total ≤ m₁` is **not a valid minorant inequality**: the true relation runs the other way (`m₁² ≤ M₀·Φ_total`).

**Two witnesses that `G″>0` does not imply `m₁>0` [T, verified [V-here] §A.5]:**
- **Task's witness** `dν = √s ds` on `[0,Λ]`: `χ = ∫₀^Λ s^{−1/2}ds = 2√Λ < ∞`, `G″ = 1/(2√Λ) > 0`, but it is purely a.c. ⇒ `f₁² = 0` and `inf supp ν = 0`. (For `Λ=4`: `G″=0.2501`, `m₁=0`.) A positive curvature coexists with a **gapless** edge.
- **Thin atom** `ν = w δ_ε + (1−w) δ_S`, `w/ε` bounded: `G″>0` while `m₁²=ε→0`.

**Therefore: `Φ_total` certifies `G″ > 0` = H1 (= `χ<∞`), and NOTHING about a spectrum floor.** The optimization succeeds at "positive at `B=0`" and still gives **H1, not `H*`.** This is the "width-not-floor / H1-not-`H*`" death, now in the cleanest curvature language.

### A.3 The named missing inequality, derived and shown to BE `G″>0` = H1 [T]

Impose `G″>0` on the exact split: `Δ_dyn + W_geo + X > 0`, i.e.
```
W_geo  >  −Δ_dyn − X .                             (◆')
```
Equivalently, using the fiber-defect split `G″ = S − D` (`INVENT-VARIANCE-INEQ.md`; `S` Gaussian-Schur, `D>0` the Savvidy fiber defect, `D/S≈2%`), `G″>0 ⟺ S>D`. The KEY-DIRECT §6(d) form `W_geo > D − Δ_dyn − X` is the same condition with the Savvidy defect `−D` made explicit (`D>0` makes it marginally stronger than (◆')). **All of these are literally `G″>0`** — an identity rewritten, not an inequality. Proving any of them uniformly is `H1` (`Var_μ(∂_nV|0) < E_μ[h_n|0]`). The `K⁻¹` Schur term cancels between `S` and `D` (`KEY-DIRECT-DERIVATION.md` §2), so the nonabelian indefinite fiber gives **no** separating lever. **[T]** The margin `D/S≈2%` is a *measurement of Gaussian-exactness*, not a proof.

### A.4 Iron laws + acid test on `Φ_total` (passes in form; fails only the minorant)

| test | `Φ_total = √(2gB) + W_geo` | verdict |
|---|---|---|
| **IL-(i) non-positivity / not free-field-satisfiable** | signed Savvidy part (instability sign) + Haar positivity; **both = 0 for free U(1)** ⇒ `Φ_total(free)=0`; free **massive composite** `:φ²:`: no Haar, no gauge negative mode ⇒ `0` (gapped, `f₁²=0`, not falsely certified) | **PASS** (acid-clean: `0` on both free cases) |
| **IL-(ii) nonabelian** | Haar non-flat only for SU(N) ⇒ `W_geo>0` only nonabelian; negative mode nonabelian-only | **PASS** |
| **IL-(iii) no small parameter** | defined directly at `β=2.5`; `W_geo` measured there | **PASS** |
| **minorant `Φ_total ≤ m₁`** | `Φ_total` is a `G″`; squeeze gives `m₁² ≤ M₀·G″` (wrong direction) and `G″>0` is compatible with `m₁=0` (witness `√s ds`) | **FAIL** — this is the exact death |

So the optimization **repairs the free-field / acid properties** (which the bare `Φ_Sav` already had) but **does not cross the floor**: the one step it must pass — `Φ_total ≤ m₁` as a genuine spectrum minorant — is unprovable and, by the squeeze direction, false in general.

### A.5 Numerics on SU(2)/U(1) small patches [V-here]

**`W_geo` (single link, pure Haar marginal of `a=⟨v,X⟩`).** Flat chart `X∈su(2)`, `|X|=φ∈[0,2π]`, Haar density `J(φ)∝(sin(φ/2)/(φ/2))²`; marginal `p(a)∝∫_{|a|}^{2π} sin²(φ/2)/φ dφ`; `W_geo=V_eff″(0)`.
- SU(2): `p(0)=1.2188`, `W_geo = 0.2051` (analytic cross-check `(1/4)/p(0)=0.2051`, exact match). **`>0`.**
- U(1): flat Haar ⇒ `p(a)=const` ⇒ `W_geo = 0`. **`=0`.**

This is the concrete **nonabelian-only positivity**: the Haar-Jacobian curvature is strictly positive for SU(N), identically zero for U(1). It is the cleanest positive invariant the construct can build — and it is a curvature.

**Signed Savvidy spectrum.** U(1): `min E²=0`, **0** negative modes. SU(2)@`gB=1`: `min E²=−1.000`, **2** negative modes (`E²=−gB`). `Φ_Sav(B)=√(2gB)`: `1.414` at `B=1`, `0.141` at `B=0.01`, `0.000` at `B=0` — monotone `→0`.

**Squeeze witnesses.** `dν=√s ds` (`Λ=4`): `χ=3.999`, `G″=0.250`, `f₁²=0`, `m₁=0`. Thin atom `wδ_ε+(1−w)δ_S`: `G″∈{0.80,0.19,0.05}>0` with `m₁²=ε∈{0.02,0.02,0.001}`. **`G″>0` with `m₁→0` throughout.**

---

## §B The channel — `J^{PC}` of the Savvidy negative mode: 1-form magnetic (`J=1`), NOT `0⁺⁺`

The task's decisive check: is the negative mode a **1-form / magnetic (`J=1`)** excitation, or does it carry a nonzero **`0⁺⁺` (scalar, A1g)** projection? The glueball atom `f₁²` is the `0⁺⁺` residue. Four exact facts settle it.

**B.1 Color — adjoint, zero `0⁺⁺` overlap [T].** The negative mode is an eigenvector of the Hessian in the **adjoint** color representation: `ad(T₃)` on `su(2)` has eigenvalues `{0,±i}`, so the unstable gluons carry charge `q=±1` (chromo-charged). The atom operator `O=¼trF²` is a **color singlet**. By Schur/color conservation `⟨singlet|adjoint⟩ = 0`: **a single adjoint gluon has identically zero projection onto any color-singlet operator**, hence onto `f₁²`. This is the strongest possible channel mismatch and it is gauge-independent. **[T, V-here]** (`K0-OVERLAP-AUDIT.md` §b independently records the same decoupling.)

**B.2 Spin / helicity — `J≥1`, no `J=0` [T].** The unstable level is the lowest Landau level with `s_z=−1` (spin **anti-parallel** to `B`), one of the **two physical transverse** polarizations of a massless spin-1 field. A massless spin-1 has helicities `λ=±1`; a one-particle state therefore has `J ≥ |λ| = 1`. **There is no `J=0` one-gluon state** — exact, by the little group. The `0⁺⁺` (A1g, `J=0`, scalar) channel is unreachable by one negative-mode gluon.

**B.3 The LLL band spans all `J_z`, but `J_z=0 ≠ J=0` [V-here].** In symmetric gauge the LLL is infinitely degenerate (guiding center): `ψ_m ∝ z^m e^{−|z|²/(2ℓ²)}`, with total `J_z = m + s_z`. With `s_z=−1` the band is `J_z = m−1`, `m=0,1,2,…` = `{−1,0,1,2,…}` — so it **does** contain a `J_z=0` member (`m=1`). But each member is still a single **adjoint transverse vector** with physical `J≥1`; a `J_z=0` state with `J≥1` is *not* a `J=0` state. (The band's `J_z` spread is the Landau degeneracy, not scalar content.)

| `m` | `J_z=m−1` | is it `J=0`? |
|---|---|---|
| 0 | −1 | no (J≥1) |
| 1 | 0 | no (J_z=0, but J≥1) |
| 2 | +1 | no (J≥1) |
| 3 | +2 | no (J≥1) |

**B.4 Two-gluon `0⁺⁺` is a CUT, not an atom; and the mode is absent at `⟨B⟩=0` [T, V-here].** The minimal **color-singlet `J=0`** object is a **two-gluon** state. Its spectral support is a **continuum** starting at the two-particle threshold `2m_eff = 2√(gB) > 0`:
```
gB=0.5 → 2m_eff=1.414 ;   gB=1 → 2.000 ;   gB=2 → 2.828   (0⁺⁺ support ⊂ [2m_eff, ∞): a CUT).
```
An **isolated `0⁺⁺` atom** — the glueball *pole below the two-gluon cut* (`PROOF-NO-DISSOLUTION-KINEMATIC`) — is **not** produced by the negative mode; the negative mode contributes only to the `0⁺⁺` continuum **above** threshold. Worse: at the **physical vacuum `⟨B⟩=0`** the fluctuation operator is `−D²+2adF_B = −D² ⪰ 0` — **the negative mode does not exist.** The "channel bridge" evaporates exactly where the gap is measured.

**§B verdict.** The Savvidy negative mode is a **colored (adjoint), transverse-magnetic, `s_z=−1`, 1-form (`J≥1`) single-gluon excitation**. Its continuum `J^P` is that of a gluon (`1⁻`); after color-neutralization its *magnetic* partner is a **vector** (1-form) glueball, `J^{PC}=1^{+-}`-type (a magnetic vector), never `0⁺⁺`. The `0⁺⁺`/A1g projection is **zero** at one-particle level (color + helicity), and only a threshold **cut** at two-particle level. **[T]** So even granting a hypothetical `Φ_total ≤ m₁`, it would **bound the wrong channel** — a colored magnetic 1-form, not the `0⁺⁺` atom. **Structural channel mismatch.**

---

## §C The H1 → `H*` upgrade test — can an exact identity FORCE the atom `f₁²>0`?

The best `Φ_total` gives only `G″>0` (H1). Can combining it with an exact identity — A3 Fredholm determinant `f₁²=−λ₁ Res_{z=z₁}G(z)`, or B6 complete-monotone — **upgrade** H1 to the atom `f₁²>0`? Prior evidence (`META-GAPFORCING-IDENTITY.md`) says identities READ but do not FORCE. We test and give the exact reason.

**C.1 `G″`/`χ` cannot determine the atom — same `χ`, different `f₁²` [T, V-here].** Construct two KL measures with **identical** `χ` (hence identical `G″=1/χ`):
- (a) atom + continuum: `f₁² = 0.719`, `m₁² = 1`, `χ = 0.750`.
- (b) pure a.c. `√s ds` on `[0,0.141]`: `f₁² = 0`, `m₁² = 0`, `χ = 0.750`.

**Same `χ`, same `G″`; atom `f₁² = 0.719` vs `0`; gap `m₁² = 1` vs `0`.** A functional of `G″` alone — which is all `Φ_total` (a curvature) is — **provably cannot output** `f₁²` or `m₁`. H1 is a **moment** (`χ = ∫s⁻¹dν`); the atom is an **edge point evaluation**; a moment is blind to a point (linear-data / moment-blindness, `META-GAPFORCING-IDENTITY.md` §3.3c). **[T]**

**C.2 Fredholm determinant reads, does not force [NG].** For a discrete transfer operator `T_O` (spectrum in `(0,1)`):
- pure cut: `D_O(z)=det(1−zT_O)` has **no zero** for `z∈(0,1)` (eigenvalues only on the cut);
- atom `μ=0.2` + cut: `D_O(z)` has an **isolated zero at `z₁=1/μ=5`**, i.e. `z₁ = e^{+a m₁}`.

So `D_O` **reads** the atom (an isolated zero appears — exactly the tool's value, `INVENT-NEW-IDENTITY.md`). **But its location `z₁` IS the mass `m₁`** — locating it uniformly in `a` (the minimal zero in the continuum) is precisely the transport / H1 step (`LOC-DETERMINANT`; `META-GAPFORCING-IDENTITY.md` §4–5). The identity does **not** force `f₁²>0`; it **relocates** the unknown. Complete-monotone (B6) is the same: `r(t)=−(1/t)log(C(t)/M₀)` is a decreasing **upper** bound on `m₁²`, and the atom is an *entire* function with no singularity (blind to analytic continuation, `TOOL-COMPLETELY-MONOTONE`). **[NG]** for the identity class.

**§C verdict.** No upgrade. H1 + exact identity does **not** give the atom. Precise reason: H1 is a moment (blind to the edge point), and the only atom-reading identity (Fredholm) requires the zero location `z₁=e^{am₁}`, whose uniform evaluation *is* the transport/H1 problem. **READ ≠ FORCE.** So the optimization remains at H1.

---

## §D Precise verdict and the exact structural death

**Question.** Does ANY optimization of the Savvidy construct give a genuine `0⁺⁺` atom floor `f₁² ≥ c>0` (or `m₁² ≥ λ>0`, not merely `G″>0`) at `β≈2.5`, uniform `(L,a)`?

**Answer. NO.** The literal sub-goal — a functional strictly positive at `B=0` — is **achieved** (`Φ_total(0)=W_geo>0`, [V-here], nonabelian-only). But it cannot be a floor. The death is a **conjunction of four precise structural facts**, each independently fatal to the floor and jointly exhaustive of what the construct can build:

1. **H1-not-`H*` (the positive piece `W_geo` is a curvature).** `Φ_total` at `B=0` is `G″=1/χ`, the inverse susceptibility. The exact squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` makes `G″` an **upper** bound on `m₁²` (harmonic mean ≥ min) — the wrong direction — and blind to the atom: `dν=√s ds` has `G″>0`, `f₁²=0`, `m₁=0`. §A. **This is the single most important line.** A curvature is not a mass; the floor needs the residue `H*` (`f₁²≥c`), which no curvature supplies.
2. **Wrong channel (the signed piece is a colored 1-form magnetic mode).** The Savvidy negative mode is adjoint (chromo-charged), transverse-magnetic, `s_z=−1`, `J≥1`. Zero `0⁺⁺`/A1g projection at one-particle level (color singlet⊕adjoint = 0; helicity ±1 excludes `J=0`); only a two-gluon **cut** above `2m_eff` at two-particle level; and absent at `⟨B⟩=0`. So the signed negativity lives in a different (`1`-form/vector, adjoint) channel than the `0⁺⁺` atom `f₁²`. §B. **Even `Φ_total ≤ m₁` would bound the wrong channel.**
3. **Width-not-floor (the Savvidy scale vanishes where the gap is largest).** `√(2gB)` is monotone `→0` at `B=0`; the physical vacuum has `⟨B⟩=0`, `m₁>0`, so `Φ_Sav(0)=0≪m₁`. It is the curvature/width of an **unstable saddle**, not a ground-state spectral floor. Converting it = Savvidy-transport = H1/(CP). §A/§0.
4. **No small parameter.** At `β≈2.5` (`α_s~O(1)`) there is no small-`g` expansion and the free `H₀` is gapless: the only dynamical lever (`‖V‖ < gap(H₀)`) is absent (`KEY-DIRECT-DERIVATION.md` §3). The optimization changes nothing here.

**Which is "the" reason?** For *this specific construct* the sharpest, most constructive statement is **the pairing (1)+(2)**: the construct's **positive** ingredient (`W_geo`) is a *curvature* (H1 — cannot see an atom or a floor), and its **signed** ingredient (the Savvidy negative mode) is in the *wrong channel* (colored 1-form magnetic, not the `0⁺⁺` scalar). Optimizing cannot fix either: making the functional positive at `B=0` turns it into a `G″` (kills the floor by (1)); making it sensitive to the signed negativity ties it to an adjoint 1-form mode that (a) is not `0⁺⁺` and (b) disappears at `⟨B⟩=0`. These are **type mismatches**, not estimates to be sharpened. **[T]**

**What genuinely survives (bank honestly).**
- **[V-here] `Φ_total(B)=√(2gB)+W_geo` is strictly positive at `B=0`** with a **provably nonabelian-only** positive part (`W_geo^{SU2}=0.205>0`, `W_geo^{U1}=0` on a single link; Prékopa `≥N/12` at `β=0`). A clean acid-passing functional.
- **[T] The channel of the Savvidy negative mode is definitively `1`-form magnetic (adjoint, `s_z=−1`, `J≥1`), zero `0⁺⁺` at one-particle level.** A sharp, reusable discriminator (with `K0-OVERLAP-AUDIT`).
- **[T]/[V-here] Two crisp counterexamples that curvature-positivity ≠ gap:** `dν=√s ds` and the thin atom — `G″>0`, `m₁→0`.
- **[NG] The H1→`H*` identity upgrade is closed for the identity class** (same `χ`, different atom; Fredholm reads-not-forces).

**Honest bottom line.** The optimization does everything it can *within its type*: it builds a strictly-positive-at-`B=0`, nonabelian-only, acid-clean functional, and it certifies `G″>0` (H1) with clean ingredients. It **cannot** produce a `0⁺⁺` atom floor, because (1) its positive part is a curvature (H1, an upper bound on `m₁²`, blind to the atom) and (2) its signed part is a colored 1-form magnetic mode (not `0⁺⁺`, absent at `⟨B⟩=0`). The remaining gap is still exactly `H*` (`f₁²≥c` + isolation) **and** the channel bridge to `0⁺⁺` — i.e. the same nonperturbative transport machine the whole campaign needs, which **does not exist**. **Clay remains open; 0%. No QED on anything unproven.**

---

## §E References

**Witness.** `research/unconditional/w7/key_direct_optimize.py` (`py -3.12`, single small job): A1 single-link marginal curvature SU(2)/U(1); A2 exact Savvidy spectrum; A3 squeeze witnesses (`√s ds`, thin atom); B channel (`ad(T₃)` charges, `s_z`/helicity, LLL `J_z` band, two-gluon threshold); C same-`χ`-different-`f₁²` and Fredholm read-not-force toy.

**Internal (absolute paths).** `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall=Savvidy concavity, key=`G″(0)`/H1); `D:\Dev\math\research\unconditional\KEY-DIRECT-DERIVATION.md` (Attack I; the death; §2 `S−D` cancellation; §3 variational; §6(d) named missing inequality); `docs\new-tools-registry.md` §B8; `docs\YANG-MILLS-MASTER.md` §5–6, §9; `PROOF-JACOBIAN-SAVES.md` (`G″=Δ_dyn+W_geo+X`, `W_geo≥0` Prékopa, `Φ″≥N/12`, measured `W_full≈+0.22–0.31`, LB); `INVENT-VARIANCE-INEQ.md` (`G″=S−D`, `D>0⇔K indefinite`, `D/S≈2%`); `H1-IS-IT-THE-GAP.md` (squeeze `f₁²/χ≤m₁²≤M₀/χ`, `G″=1/χ` inverse susceptibility, thin-atom `wδ_ε+(1−w)δ_S`, `√s`-type blindness); `K0-OVERLAP-AUDIT.md` (§b Savvidy negatives = transverse 1-forms adjoint, `0⁺⁺` decouples by color/C/composite); `META-GAPFORCING-IDENTITY.md` (READ≠FORCE; every identity satisfied by gapless free field; Fredholm reads atom but zero-location = transport); `INVENT-NEW-IDENTITY.md` (Fredholm `f₁²=−λ₁Res`, pole-below-cut); `TOOL-COMPLETELY-MONOTONE.md` (B6 `r(t)↓m₁²` upper bound, atom entire ⇒ analytic-blind); `PROOF-NO-DISSOLUTION-KINEMATIC.md` (pole-below-cut structure); `SCHUR-CONVEX.md` (positivity ⇒ upper bounds); `LOC-DETERMINANT.md` (zero-location = transport).

**External `[R]`.** Savvidy, *Phys. Lett. B* **71** (1977) 133; Nielsen–Olesen, *Nucl. Phys. B* **144** (1978) 376; Prékopa (1971); Brascamp–Lieb, *Adv. Math.* **22** (1976) 366; Källén (1952)/Lehmann (1954); Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Shen–Zhu–Zhu, *CMP* **400** (2023) 805.

*Status: **Optimization attempted and carried to its exact death. The literal goal — `Φ_total` strictly positive at `B=0` — is achieved [V-here]; but `Φ_total(0)=W_geo` is a curvature (`G″=1/χ`) that UPPER-bounds `m₁²` and is blind to the atom (`√s ds`: `G″>0`, `m₁=0`), and the Savvidy signed part is a colored 1-form magnetic mode (`J≥1`, adjoint, `s_z=−1`) with zero `0⁺⁺` projection. Structural death = H1-not-`H*` AND wrong-channel AND width-not-floor AND no-small-parameter. H1→`H*` identity upgrade fails (READ≠FORCE). Clay 0%. No QED.***
