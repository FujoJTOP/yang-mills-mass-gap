# SLOT-DSE-FIXPOINT — the Dyson–Schwinger gap equation at `k→0`: the **scheme-independent** piece is the *seagull identity that certifies the gapless solution*; the only sign-definite exclusion of `m²=0` is the vertex-pole residue `C(r²)`, which is the unprotected order parameter (=`H1`)

> 2026-11-24 · Route: **Dyson–Schwinger / 2PI gap equation** for the `0⁺⁺` sector, `k→0`. Target: the ONE open slot of `PROFILE.md` §3 — a **`k→0` sign-definite self-consistency with no gapless fixed point**, uniform at the crossover, non-positivity, nonabelian, acid-clean. This note attacks it through the DSE.
>
> **判词 (Verdict).** The exact gluon DSE is a hierarchy that never closes; the PT-BFM/PT **scheme-independent** piece of it (the STI-fixed / background-quantum Ward content) is precisely the **seagull identity**, and for vertices regular at `q=0` the seagull identity **forces `Δ⁻¹(0)=0`, i.e. `m²=0`** (Aguilar–Binosi–Figueiredo–Papavassiliou 2016; verified). So the one piece of the DSE that is *both* scheme-independent *and* free-field-acid-clean is, under the acid test, the piece that **certifies the gapless solution** — it is an identity satisfied by the free field (`META-GAPFORCING §3.1`). The *only* way the DSE produces `m²(0)=c>0` is to relax the regularity assumption and insert a **longitudinally-coupled `1/q²` vertex pole** (the Schwinger mechanism, Cornwall 1982 / ABP), whose residue `C(r²)` is **not fixed by the STI** (no anomalous slot, `W3-marginal-BPHZ`) and is determined by a *homogeneous* Bethe–Salpeter equation whose *scale* is fixed by a separate inhomogeneous vertex SDE — i.e. it **is** the dynamical order parameter, `=H1`. Hence: **there is no scheme-independent, sign-definite DSE piece that excludes `m²=0`; the two requirements are in exact complement** (Theorem D below). **[NG]** for the "scheme-independent DSE exclusion" class; **[O]** the exact missing step = the pole residue `C(r²)`. **No 5/5. Clay 0%. No QED.**

---

## §0 · Setting, the two objects, labels, and the iron laws

**Label key.** **[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for a *named class* · **[V-here]** new this round · **[R]** inherited/banked.

**The target (Clay, `PROFILE.md`).** `H` self-adjoint transfer Hamiltonian on `H_phys=ker Gauss`, `HΩ=0`, bottom of `(vac)^⊥` `=m₁`, KL measure `ν`, atom `f₁²=ν({m₁²})`. The gap `= H*` (edge atom + isolation) `=` the `k→0` sign-definite self-consistency without a gapless fixed point (`PROFILE §3`). The task asks whether the DSE supplies that self-consistency.

**Two DSE objects — and an exact channel mismatch to respect.**
- **(i) The gluon DSE** `D⁻¹_{μν}=D₀⁻¹_{μν}+Π_{μν}` is the **`1⁻⁻`** (vector) self-energy. Its IR value defines the *gauge-fixed dynamical gluon mass* `m_g²=Π_T(0)`. The gluon is a **vector**, `J^{PC}=1⁻⁻`.
- **(ii) The physical `0⁺⁺` channel** is the scalar composite `O=¼ tr F²`; its mass `m_{0⁺⁺}²` is the pole/edge of `⟨O O⟩_c`. The two are linked only through `1⁻⁻⊗1⁻⁻⊃0⁺⁺` and the kernel (`BS-KERNEL-SIGN`).
- **Mismatch [T].** The task's phrase "the DSE for the `0⁺⁺` channel … `m²=Σ(m²)`" is the **gluon-DSE language**, which is the **`1⁻⁻`** object. `m_g` is (a) gauge-fixed, (b) not a positive-norm asymptotic state (`A∉H_phys`), (c) not the `0⁺⁺` mass. **The bridge `m_g>0 ⇒ m_{0⁺⁺}>0` is exactly the transfer wall** (`PROOF-FRG-DSE-GAP-EQUATION`, `PROOF-INDUCED-MASS-GAUGE-INVARIANT`). I attack both objects and the bridge; the verdict is the same on both sides.

**The three iron laws (`KEY-DIRECT-DERIVATION §0`).** Any lower-bound derivation must be **(IL-i) non-positivity** (the free massless field must *fail* it; identities/positivity are Schur-concave ⇒ only upper-bound the atom), **(IL-ii) uses the nonabelian interaction** (else the U(1) no-go), **(IL-iii) no small parameter** (must hold at `α_s~O(1)`, `β≈2.5`). **Free-field acid test:** clean iff `Φ(free massless) ≤ 0` (not `>0`), and *discriminating* iff it also gives `≤0` on the gapped-but-atom-free `:φ²:`.

---

## §1 · The DSE written concretely (gluon + ghost + the exact 2PI kernel), and its scheme-independent piece

### 1.1 The gluon DSE, transverse part, full vertices `[T, standard]`

In Landau gauge, with the transverse projector `P_{μν}(k)=δ_{μν}−k_μk_ν/k²`, the exact one-particle-irreducible gluon self-energy is `Π_{μν}=Π_T P_{μν}` and
```
D⁻¹_T(k) = k² + Π_T(k) ,                                          (1)
Π_T(k) = ½ g² C_A ∫ d⁴q/(2π)⁴ Γ^{(0)}_{μαβ} D_{αα'}(q) D_{ββ'}(k−q) Γ_{να'β'}(q,k)
        −  g² C_A ∫ d⁴q/(2π)⁴ q_μ (k−q)_ν D_G(q) D_G(k−q) Γ^{gh}_{α}(q,k) Γ^{gh}_α(k,q)   (2)
        +  (4-gluon tadpole / seagull) ,
```
where `D_{αα'}` is the full gluon propagator, `D_G` the full ghost, `Γ_{να'β'}` the **full** three-gluon vertex, `Γ^{gh}` the full ghost-gluon vertex (`Roberts–Williams 1994, hep-ph/9403224`; `Cornwall 1982`). The **`0⁺⁺`-channel analogue** is the DSE/BSE for the composite `O=¼trF²`:
```
G_O(k) = ⟨O(k)O(−k)⟩_c ,   G_O⁻¹(k) = (bare) + Λ^{(2)}_O[ D , D_G , Γ , Γ^{(4)} ] ,   (3)
```
`Λ^{(2)}_O` the **exact 2-particle-irreducible kernel** of the `O`-channel (the amputated, irreducible 4-point function).

### 1.2 The "exact 2PI kernel" is not a closure `[T]`

The 2PI effective action `Γ_2PI[φ,G]` gives `G⁻¹=G₀⁻¹+Π` with `Π=½ Λ^{(2)}`. Stationarity `δΓ_2PI/δG=0` yields an equation for the **full** `Λ^{(2)}` that itself sources the **irreducible 6-point** kernel, and so on. Equivalently, in (2) the full `Γ, Γ^{gh}, Γ^{(4)}` are unknown; no Ward identity fixes their **transverse** parts. This is the lattice **Makeenko–Migdal** open hierarchy (`cpb-attack/10-dse`, `[T]`): the 2PI formulation **trades a vertex-truncation for a kernel-truncation**, and the kernel is another infinite series. **The DSE is exact only as an infinite hierarchy; any solve is a truncation.** `[T]`.

### 1.3 The scheme-independent piece — identified, and it is the **seagull identity** `[T/R]`

Contract the background leg of the three-gluon vertex with `q_μ`: the **Abelian STI / background-quantum identity (BQI)** gives a QED-like Ward identity **provided no `1/q²` poles are present** (`Binosi–Papavassiliou, Phys. Rept. 479 (2009) 1, arXiv:0909.2536`; `Aguilar–Binosi–Papavassiliou, PRD 78 (2008) 025010, arXiv:0802.1870`). Feeding these into (2) produces an exact diagrammatic cancellation whose **residual is the seagull identity** — a total-derivative/translation-invariant integral that **vanishes under any invariant regularization** (`Aguilar–Binosi–Figueiredo–Papavassiliou, PRD 94 (2016) 045002, arXiv:1604.08456`, verified):

> **Seagull identity [T/R].** For vertices **regular at `q=0`**, the STI+BQI-fixed content of the PT-BFM gluon DSE reduces to an identity that forces
> ```
>     Δ⁻¹(0) = 0      i.e.   m²(0) = 0 .
> ```
> It "is valid independently for each set of diagrams" (`ABP 2016`, `Binosi funcQCD-2018`); the only escape is a **vertex `1/q²` pole**.

This is the key structural fact. The scheme-independent (STI-fixed), regularization-independent, gauge-invariant piece of the DSE is **precisely the piece that certifies `m²=0`.** The PT mass formula `m²(k²)=∫…C(y)…Δ²(y)` of `Aguilar–Binosi–Papavassiliou (PRD 78 (2008) 025010; Front. Phys. 11 (2016) 111203, arXiv:1511.08361)` is obtained **only after** inserting the vertex pole; it is gauge-invariant in PT-BFM but its **existence** is the pole, not the STI.

### 1.4 The candidate "sign-definite scheme-independent pieces" — each examined `[T]`

| Candidate | scheme-independent? | sign-definite / excludes `m²=0`? | verdict |
|---|---|---|---|
| **Longitudinal Ward / BQI** `q_μΠ^{μν}=0` | **yes** (`[T]`, BQI) | **no** — it is an *identity*; free field satisfies it; seagull ⇒ `m²(0)=0` | **acid-vacuous; certifies the gapless solution** |
| **Pinch-technique self-energy** `Π^{PT}` | **yes** (gauge-invariant, PT-BFM, `[T/R]`) | **no** — transversality + seagull give `Π^{PT}(0)=0` for regular vertices | **vacuous without the pole** |
| **PT effective charge** `α_s^{PT}(k)` | **yes** (RG-invariant to 2 loops, `[T/R]`) | **no** — a spacelike observable; IR *freezing* is a finite coupling, **not** a pole/atom; no `m` | **freezing ≠ gapping `[J]`** |
| **Trace anomaly / RG-invariant** `⟨T^μ_μ⟩=(β/2g)⟨trF²⟩` | **yes** (`[T]`) | **no** — aggregate UV/moment data; **upper-bounds** the atom (`GAP-G2-edge`) | **wrong direction (IL-i)** |
| **Vertex-pole residue** `C(r²)` (Schwinger) | **NO** — unprotected transverse vertex, no anomalous slot (`W3-marginal-BPHZ`) | **YES** — `m²(0)=c=∫C` `>0` | **scheme/vertex-dependent = the order parameter** |

**§1 verdict.** The scheme-independent column is entirely **acid-vacuous** (identities / positivity / upper-bounds / no-pole). The only **sign-definite exclusion** of `m²=0` lives in the last row, `C(r²)`, which is **not** scheme-independent. This is not an accident; it is forced (§2).

---

## §2 · The core theorem — the **acid–scheme dichotomy** `[T]`

> **Theorem D (DSE acid–scheme dichotomy).** Partition the content of the truncated DSE (2)/(3) into **(U)** the STI/BQI-fixed (longitudinal, Ward, BQI, renormalization-constant) part and **(S)** the transverse completion (the full transverse vertex `Γ_⊥`, the pole residue `C(r²)`, the IR boundary condition). Then:
> **(i)** `U` is scheme- and gauge-independent (BQI, `[T/R]`), but every relation of `U` is an **identity** satisfied by the free massless field; by the **seagull identity**, `U` forces `m²(0)=0` — `U` is **acid-vacuous and certifies the gapless solution**. `[T]`
> **(ii)** The dynamical mass `m²(0)=c` is **defined entirely in `S`**: `c` is the `1/q²` pole residue of the transverse vertex, which `U` does not constrain (no anomalous slot). Hence the **exclusion of `m²=0` lives entirely in `S`**, and `S` is scheme/vertex/boundary-dependent — the **scaling** vs **decoupling** solutions of the *same* STI-satisfying system are two `S`-choices. `[T]`
> **(iii)** Therefore **no piece of the DSE is simultaneously scheme-independent and sign-definite-excluding-`m²=0`**: scheme-independent `⇒` in `U` `⇒` acid-vacuous; excludes `m²=0` `⇒` in `S` `⇒` scheme-dependent. The two requirements are in **exact complement**. `[NG]` for the "scheme-independent sign-definite DSE piece" class.
>
> *Proof.* (i) is `META-GAPFORCING-IDENTITY §3.1` (every exact identity is satisfied by the free/gapless field) specialized to the STI/BQI, plus the seagull identity (§1.3, `ABP 2016`). (ii) is definitional: `m²(0)≡` the part of `Π(0)` not fixed by `U`, and the seagull identity says the `U`-part is zero; the residue `C(r²)` is fixed by the inhomogeneous vertex SDE + homogeneous BSE whose scale is unprotected (`ABP 2008/2016/2022`, `W3-marginal-BPHZ`). (iii) is the immediate corollary: `U ∩ S = ∅`. ∎

**This is the precise form of the open slot.** The `PROFILE §3` slot asks for a relation that is *controller-free, intrinsic (⇒ built from `U`), and excludes the gapless solution*. Theorem D says **intrinsic/STI-fixed and excludes-gapless are complementary** at the DSE level: the intrinsic part *is* the gapless-certifying identity. The exclusion must come from `S` `=C(r²)`, which is the order parameter `=H1`. [T]/[NG]/[O].

---

## §3 · Gapless fixed-point test — **does `m²=0` solve the DSE? Yes, always; the STI does not exclude it — it forces it**

**3.1 `m²=0` is always a solution `[T]`.**
- **Free branch:** `Π_T=0 ⇒ D⁻¹=k² ⇒ m²(0)=0`. The free massless field is an exact solution of the truncated DSE (with free ghosts/vertices). Acid test clean but **vacuous** (it excludes nothing).
- **Scaling branch:** the IR power-law `D_T~(k²)^{2κ}`, `D_G~(k²)^{−κ}`, `κ∈[½,1)` solves the DSE *with vertices scaling as `(k²)^{(n−m)κ}`* (`Fischer–Pawlowski, PRD 75 (2007) 025012, hep-th/0609009; II: PRD 80 (2009) 025023, arXiv:0903.2193`). Here `m²(0)=0` too (`Δ⁻¹(0)=0`, no IR mass) — the would-be gap is in the *enhanced ghost*, not a gluon mass.
- **Decoupling branch:** `D_T(0)=1/m²<∞` finite. **Both branches solve the same STI-satisfying system** with different `S` (transverse vertex / boundary condition). `Fischer–Pawlowski` argue the *scaling* branch is unique **given** DSE+FRG consistency *and* the boundary conditions; that uniqueness is itself an `S`-choice, and the decoupling solution *exists* and matches the lattice propagator (`ABP 2008`). **The DSE does not by itself exclude the `m=0` branch.** `[T]` (existence) / `[NG]` (no selection).

**3.2 Does the nonabelian interaction or ST identity exclude `m²=0`? No — on both counts `[T]`.**
- **The STI does not exclude it; the seagull identity *forces* it** for regular vertices (§1.3). The STI is on the side of the gapless solution.
- **The nonabelian interaction is necessary but not sufficient.** The vertex pole is `∝g²C_A` (`C_A=N>0`, nonabelian — so `[IL-ii]` passes *in form*, U(1) has no such pole), but **no identity forces the pole**. Its residue `C(r²)` solves a **homogeneous** BSE (`m=C·K·…`, always admitting `C=0`), normalized by a separate inhomogeneous SDE — so `C=0` (gapless) is always a solution of the homogeneous BSE, and `C>0` is a *spontaneous* choice. `[T]`.

**3.3 Reduced-model witness `[V-here]`.** `w7/dse_fixpoint_witness.py` solves the 0-dim reduced gap equation `x=λ K(x) x` (`x=m²/scale²`):
- **regular vertex** `K=1/(1+x)`: `x=0` always solves; a positive branch `x*=λ−1` appears **only for `λ>1`** (a vertex-strength threshold), vanishing at `λ_c=1`.
- **pole vertex** `K=1/x+1/(1+x)`: `x*>0` for **every** `λ>0`, `x*→0` as `λ→0`.
Both have `x=0`; **which positive branch exists is decided by `K`, i.e. by the vertex/scheme.** Fixed-point iteration confirms `x=0` is the contraction attractor for `λ<λ_c` ("the provable branch collapses to the gapless solution," matching `MS-GAP-EQUATION`). *(Schematic; no YM number asserted.)*

**§3 verdict.** `m²=0` survives the DSE in both branches; the STI and the seagull identity *certify* it; the nonabelian interaction enables but does not force the pole. **[O]** the force is `C(r²)`.

---

## §4 · Uniform attainment across `β≈2.5` — a category statement plus a transmutation obstruction

**4.1 The DSE carries no `β` `[T].`** `β` is the lattice coupling; the continuum DSE has only the renormalized coupling `α_s(μ)`/the scale `Λ`. "Uniform across `β≈2.5` **without the truncation's ansatz**" is the statement that the positive branch is scheme-independent *and* survives the lattice→continuum matching. **That matching is the crossover certificate** (`PROOF-FRG-DSE-GAP-EQUATION §2`). The DSE cannot "see" `β`; the only bridge is a scheme-match to the lattice family `(β,L,a)`, i.e. `H1` itself. `[T]` (category) / `[O]` (the bridge).

**4.2 Even within the DSE, the mass is a transmutation, not a uniform constant `[T].`** Both branches have
```
    m²(0) ~ Λ² exp(−c/α_s)     (dimensional transmutation; c an O(1) number)
```
(`ATOM-FRG-mass`; `Cornwall 1982`). The exponent coefficient `c` is **scheme-dependent** (it is set by the exact two-loop `β`-function normalization and the vertex), and `m=0` at `α_s→0`. So the positive branch is **not bounded below uniformly in the coupling**: near the (vertex-dependent) critical strength it shrinks to `0` (`§3.3`: `x*→0` as `λ→λ_c` / `λ→0`). The "uniformity without the ansatz" the slot demands is exactly the statement that this transmutation scale is stable and scheme-independent across the crossover — which is the transport/crossover wall, not a DSE computation. `[T]`/`[O]`.

**4.3 The `0⁺⁺`-channel uniformity is worse `[T].`** Even granting `m_g>0`, the `0⁺⁺` two-gluon threshold sits at `s=0` until `m_g` lifts it, and the J=0 kernel is spin-momentum-**indefinite** with the Savvidy tachyon feeding it at second order (`BS-KERNEL-SIGN §3–5`; `BS-VARIATIONAL-LOWER §5`). So the `0⁺⁺` mass is not uniformly attained from the gluon DSE; the threshold obstruction is the same `m_g>0` (=`H1`) one rung up.

**§4 verdict.** Uniform attainment is not a DSE property — it *is* the crossover/transport certificate. `[O]`.

---

## §5 · Acid test + the three iron laws

**5.1 Free-field acid test `[T].`** On the free massless field `Π_T=0 ⇒ m²(0)=0`; the DSE (and every scheme-independent piece `U`) is satisfied and returns `0`. **Acid-clean but vacuous** — it excludes nothing. This is precisely `META-GAPFORCING §3.1`: an exact identity's solution set contains the gapless free field.

**5.2 The *discriminating* acid test (`:φ²:`, gapped but atom-free) — the DSE FAILS it `[T].`** The gapped free massive field `:φ²:` has a propagator mass `m>0` but **no physical `0⁺⁺` atom** (`f₁²=0`). A DSE "mass" `m²(0)=Π(0)>0` **cannot distinguish a dynamically generated mass from an input mass** — fed the massive free field it returns `m>0` and would falsely certify a "gap," exactly the `:φ²:` failure. **Worse, the gluon DSE mass is not even a positive spectral weight:** by **Oehme–Zimmermann** (`PRD 21 (1980) 1661`), the transverse gluon obeys `∫ρ_A(s)ds=0`, so `D(0)=∫ρ_A(s)/s ds>0` is compatible with `∫ρ_A=0` **only because `ρ_A` is sign-indefinite** — `m_app²=1/D(0)` is **not** a spectral atom. `[V-here]` (W3): `ρ=δ(s−a)−δ(s−b)` gives `∫ρ=0`, `D(0)=1/a−1/b>0`, apparent mass `²=ab/(b−a)`, but the weight at `b` is `−1`. **The DSE's massive solution is a gauge-fixed, sign-indefinite object, not the `0⁺⁺` edge atom.**

**5.3 The three iron laws `[T]/[NG].`**
- **(IL-i) non-positivity — FAIL.** The scheme-independent DSE content (`U`) is an identity/positivity; it is satisfied by the free field and, via seagull, *forces* `m²=0`. It cannot be strict (Theorem D(i)). The physically-lower-bounding `0⁺⁺` correlator is sign-**definite** (`O=¼trF²≥0` ⟹ `ρ_O≥0`) and gives **upper** bounds (`KEY-DIRECT-DERIVATION` Attack IV). `[NG]`.
- **(IL-ii) nonabelian — PASS in form, but not forcing.** The pole residue `∝g²C_A=N·g²` is nonabelian (U(1) has no pole), but its existence is unforced (`§3.2`). `[T]`.
- **(IL-iii) no small parameter — PASS (it is O(1)), but the bound is not uniform.** The DSE works at `α_s~O(1)`; but `m²~Λ²e^{−c/α_s}` is not uniformly bounded and `→0` at the scheme-dependent threshold (`§4.2`). `[T]`.

**5.4 The channel dichotomy `[T]`.** `1⁻⁻` gluon DSE: sign-indefinite `ρ_A` (OZ) → a "mass" that is **unphysical, no transfer**. `0⁺⁺` physical channel: sign-definite `ρ_O` → **upper bounds only**. The two give **opposite-direction** information; neither yields a scheme-independent sign-definite **lower** bound on `m_{0⁺⁺}`. `[T]/[NG]`.

---

## §6 · Verdict — constructed or the exact failing step

**Verdict [NG]/[O]/[T].**
- **Not constructed:** there is **no** scheme-independent, sign-definite DSE relation that excludes `m²=0` (Theorem D). The scheme-independent piece is the seagull identity, which *certifies* `m²=0`; the exclusion lives entirely in the transverse vertex pole `C(r²)`, which is the unprotected order parameter. **[NG]** for the "scheme-independent DSE exclusion" class.
- **Exact failing step [O]:**
  ```
  m²(0) = c = ∫ d⁴q  C(q²) · [pole-dressed kernel] ,
  C(r²) solves the HOMOGENEOUS BSE (C=0 always a solution),
        normalized by the INHOMOGENEOUS vertex SDE;
  C(r²) ∝ g²C_A > 0 is NOT fixed by the STI (no anomalous slot, W3-marginal-BPHZ).
  ```
  The one unproven input is **`C(r²)`** — a nonperturbative vertex residue, scheme/vertex-dependent, whose uniform survival across the crossover is the transport wall. `C(r²)` **is** the order parameter `=` `MS-GAP-EQUATION`'s `F′(0)>1` `=` the edge atom `f₁²` (`PROOF-INDUCED-MASS-GAUGE-INVARIANT`: the gauge-invariant PT mass *is* the edge atom) `=` `H1`.
- **Why this is not a re-discovery of a dead name.** The new content (`[V-here]` framing) is the **acid–scheme dichotomy** (Theorem D): it converts "the DSE is truncation-dependent" (a known negative) into the **exact complementarity** *scheme-independent `⟺` acid-vacuous-and-certifies-`m²=0`* versus *sign-definite-exclusion `⟺` transverse-pole `⟺` scheme-dependent*. It identifies **the seagull identity** (§1.3, verified `ABP 2016`) as the *mechanism* by which the scheme-independent piece certifies the gapless solution — the sharpest available statement that the slot's `P4` requirement (intrinsic/controller-free **and** excludes-gapless) is **exactly complementary at the DSE level.**
- **Honest boundary.** `[NG]` is for the *named class* (scheme-independent DSE pieces), not a global impossibility (`∀class ≠ ∀possible`, `KEY-EXISTS-OR-IMPOSSIBLE §5`). Theorem D is `[T]` at the level of the STI/seagull (verified literature) + `[V-here]` schematic witness. The residue `C(r²)` is a real object; bounding its uniform survivor is the crossed-out step. **Not a proof of Clay. Clay 0%. No QED.**

**Sharpest single line.** *The STI, the pinch technique, and the seagull identity all live on the side of `m²=0` — the scheme-independent DSE is an identity that the free field satisfies and that, for regular vertices, forces the gapless solution. The only sign-definite exclusion of `m²=0` is the `1/q²` vertex-pole residue `C(r²)`, which no STI fixes and which is the dynamical order parameter. The slot's demand — "intrinsic **and** excludes the gapless fixed point" — is, at the DSE level, a demand for a scheme-independent object that is not an identity; Theorem D says that object cannot exist: the intrinsic content is the identity, the exclusion is the scheme charge, and the scheme charge is `H1`.*

---

## §7 · References (ids verified this session; internal banked labels carried)

**External (verified).**
- C. D. Roberts, A. G. Williams, *Dyson–Schwinger equations and their application to hadronic physics*, Prog. Part. Nucl. Phys. **33** (1994) 477, hep-ph/9403224. `[R]`
- J. M. Cornwall, *Dynamical mass generation in continuum QCD*, Phys. Rev. D **26** (1982) 1453. `[R]`
- D. Binosi, J. Papavassiliou, *Pinch technique: theory and applications*, Phys. Rept. **479** (2009) 1, arXiv:0909.2536. `[verified]`
- A. C. Aguilar, D. Binosi, J. Papavassiliou, *Gluon and ghost propagators in the Landau gauge: deriving lattice results from Schwinger–Dyson equations*, Phys. Rev. D **78** (2008) 025010, arXiv:0802.1870. `[verified]`
- A. C. Aguilar, D. Binosi, J. Papavassiliou, *The gluon mass generation mechanism: a concise primer*, Front. Phys. **11** (2016) 111203, arXiv:1511.08361. `[verified]`
- A. C. Aguilar, D. Binosi, E. G. de Figueiredo, J. Papavassiliou, *Unified description of seagull cancellations and infrared finiteness of gluon propagators*, Phys. Rev. D **94** (2016) 045002, arXiv:1604.08456. `[verified — seagull identity §1.3]`
- A. C. Aguilar, M. N. Ferreira, J. Papavassiliou, *Exploring smoking-gun signals of the Schwinger mechanism in QCD*, Phys. Rev. D **105** (2022) 014030, arXiv:2109.12712. `[verified]`
- C. S. Fischer, J. M. Pawlowski, *Uniqueness of infrared asymptotics in Landau gauge Yang–Mills theory*, Phys. Rev. D **75** (2007) 025012, hep-th/0609009; *II*, Phys. Rev. D **80** (2009) 025023, arXiv:0903.2193. `[verified — scaling vs decoupling]`
- C. S. Fischer, *Infrared properties of QCD from Dyson–Schwinger equations*, J. Phys. G **32** (2006) R253, hep-ph/0605173. `[R]`
- R. Oehme, W. Zimmermann, *Relation between effective couplings for asymptotically free theories*, Phys. Rev. D **21** (1980) 1661 (`∫ρ_A=0`). `[R]`
- T. Kugo, I. Ojima, *Local covariant operator formalism … confinement*, Prog. Theor. Phys. Suppl. **66** (1979) 1. `[R]`
- M. Q. Huber, *Nonperturbative properties of Yang–Mills theories*, Phys. Rept. **879** (2020) 1, arXiv:1808.05227. `[R]`
- A. K. Cyrol, M. Mitter, J. M. Pawlowski, N. Strodthoff, PRD **94** (2016) 054005, arXiv:1605.01856. `[R]`
- V. N. Gribov, *Quantization of non-Abelian gauge theories*, Nucl. Phys. B **139** (1978) 1; I. M. Singer, *Some remarks on the Gribov ambiguity*, CMP **60** (1978) 7 (Landau-gauge `D_T` is gauge-fixing / Gribov-copy dependent). `[R]`

**Internal (absolute paths; banked labels carried).**
- `D:\Dev\math\research\unconditional\MS-GAP-EQUATION.md` (Schema C: `κ=F(κ)`, `κ=0` always, `F′(0)>1` pitchfork, tachyon high-momentum).
- `D:\Dev\math\research\unconditional\PROOF-FRG-DSE-GAP-EQUATION.md` (gauge/truncation dependence; OZ sign-indefinite; transfer fails; `β` category error).
- `D:\Dev\math\research\unconditional\PROOF-INDUCED-MASS-GAUGE-INVARIANT.md` (gauge-invariant PT mass = edge atom = gap; STI/anomaly cannot force it; forced pole ⇒ Goldstone `m=0`).
- `D:\Dev\math\research\unconditional\PROFILE.md`, `PROFILE-FACET-ANCHOR.md` (the slot; `P4`; Theorem SA; gapless-solution obstruction).
- `D:\Dev\math\research\unconditional\BS-KERNEL-SIGN.md`, `BS-VARIATIONAL-LOWER.md` (`0⁺⁺=1⁻⁻⊗1⁻⁻`; kernel indefinite; threshold at 0; `ρ=4m_g²` invalid).
- `D:\Dev\math\research\unconditional\META-GAPFORCING-IDENTITY.md` (§3.1 identity-forcing no-go; §3.5 READ≠FORCE).
- `D:\Dev\math\research\unconditional\KEY-DIRECT-DERIVATION.md` (§0 iron laws/acid test; Attack IV sign-definite `ρ_O` ⇒ upper bounds).
- `D:\Dev\math\research\unconditional\W3-marginal-BPHZ.md` (no anomalous slot for the transverse vertex / `1/q²` pole).
- `D:\Dev\math\research\unconditional\cpb-attack\10-dse\REPORT.md` (Makeenko–Migdal open hierarchy).
- `D:\Dev\math\research\unconditional\w7\dse_fixpoint_witness.py` (`[V-here]` W1–W3 numeric witness).

---

*Status: **DSE-route facet, verdict [NG] for the scheme-independent-DSE-exclusion class + [O] for the vertex-pole residue `C(r²)` + [T] Theorem D (acid–scheme dichotomy) + [V-here] schematic witnesses.** The `k→0` DSE relation is `D⁻¹_T(0)=Π_T(0)=0` for regular vertices (seagull identity) — it **certifies** the gapless solution; the sign-definite mass `m²(0)=∫C` requires the unprotected `1/q²` vertex pole, the order parameter. `m²=0` survives both scaling and decoupling branches. Uniform attainment is the crossover/transport certificate, not a DSE property. **No scheme-independent sign-definite DSE piece excludes `m²=0`; the exclusion is the scheme charge `C(r²)=H1`. Clay 0%. No QED.***
