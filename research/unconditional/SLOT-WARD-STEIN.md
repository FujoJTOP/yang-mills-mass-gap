# SLOT-WARD-STEIN — Ward + Stein + variance at `k→0`: the triple is ONE identity = "the coarse marginal is a normalised pdf"; it has a gapless solution with `G″` FREE, the variance term SUBTRACTS (wrong direction), and the only nonabelian sign source is the Haar Jacobian, whose promotion to `G″≥c` is the exact missing step (`(★)=H1`)

> 2026-11-24 · **ROUTE: Ward/Stein identity as a gap equation.** Attack the ONE open slot of `PROFILE.md` (a `k→0` sign-definite self-consistency, no gapless fixed point, uniform at the crossover, non-positivity, nonabelian, acid-clean). Use the **exact** identities — Ward (IBP) `E[af]=1`, Stein `E[HessV]=E[∇V⊗∇V]⪰0`, variance `G″=E[H|0]−Var(f|0)` — which hold for **every admissible measure, no controller**, and ask whether a **combination `R`** forces a `k→0` sign that only gapped measures satisfy.
>
> Labels (charter §2): **[T]** true theorem · **[C]** conditional theorem · **[O]** exact missing step · **[NG]** no-go for a named class · **[V-here]** new this round · **[R]** inherited. **No QED. Clay stays 0%.**

---

## VERDICT (判词)

> **The Ward + Stein + variance system, at `k→0`, is a SINGLE exact identity — "the coarse marginal `p(a)` is a normalised probability density."** Writing `m(a)=E[f|a]=−(log p)′(a)`, the three are:
> ```
>   Ward   :  E[af] = ∫ a m p da = 1            ⟺  ∫ p da = 1          (normalisation; since mp=−p′)
>   Stein  :  E[H] = E[f²]  ⟺  ∫ m′ p da = ∫ m² p da  ⟺  ∫(p′)²/p = ∫(p′)²/p  (Fisher=Gram; a TAUTOLOGY)
>   Var.   :  G″(a) := m′(a) = E[H|a] − Var(f|a);   G″(0) = −p″(0)/p(0)   (definition)
> ```
> Collapsing (`K0-WARD-LOWER §2`) gives the *one* relation `R`:
> ```
>   R :  E[H − f² | a] = −p″(a)/p(a),    hence    κ · Var(a) = 1 − η,
>        κ := G″(0) = m′(0),  η := 1 − κ·Var(a)  ("non-Gaussianity defect"),
> ```
> **`R` does not constrain `p` at all, hence does not constrain `G″(0)` at all.** `G″(0)=−p″(0)/p(0)` is **free**: the gapless free field (Gaussian marginal, `Var(a)→∞`) satisfies `R` with `G″→0`; a bimodal marginal satisfies `R` with `G″(0)<0`; a flat-top marginal with `G″≈0`. **No combination of these exact identities is `≥c>0` for gapped measures only** — they are **READ, not FORCE** (`META-GAPFORCING §3.1`, specialized here).
>
> **The nonabelian steelman.** The *one* sign-definite, acid-clean, genuinely **nonabelian** input is the Haar (Jacobian) curvature `Φ″⪰N/12` (`Φ=−log J`), which the **flat/abelian free field violates** (`Φ=0`). It gives `W_geo≥0` on the geodesic ball and the **`β=0` gap `≥N/12`** **[T]**. But it bounds `W_geo`, **not** `G″`: the variance identity has the **wrong sign** — `G″ = E[H|0] − Var(f|0) ≤ E[H|0]` — so flooring `G″` needs a **lower bound on the conditioned action Hessian `E[H|0]`**, which is **Savvidy-indefinite** at `β≈2.5`. Promoting the Haar floor to `G″≥c` is *exactly* `(★)`: `E[H|0] ≥ λ + Var(f|0)` — the local face of `H1` (`PROVE-Gpp §4`, `K0-WARD-LOWER §6`).
>
> **Verdict: NOT constructed — exact failing step.** `R` = "`p` is a pdf"; gapless Gaussian satisfies it with `G″→0`, bimodal with `G″<0`; the variance identity is a subtraction. Reduces to `H1` via `(★)`. The Ward/Stein route does **not** evade the wall; it states it in the sharpest identity form: the entire gap is the single-point comparison `E[H|0] > Var(f|0)`, i.e. "the conditional action Hessian at the mode exceeds the force variance by `≥λ`." **`H1`, not a 6th species.** Clay 0%.

---

## §1 · The exact system at `k→0` (setup and identifications)

**Setting** (`K0-WARD-LOWER §1`, `MS-GAP-EQUATION §1`, `INVENT-VARIANCE-INEQ §1`). Flat link-algebra chart `X=log U`, measure `μ ∝ e^{−V}dX`, `V=S_W+Φ`, `Φ=−log J` the Haar Jacobian, `J` the volume element of `SU(N)` in exponential coordinates. Coarse mode `a=⟨v,X⟩`, `|v|=1` (`k→0` uniform); conjugate force `f=∂_aV=⟨v,∇V⟩`; `H=Hess V`, `H_vv=vᵀHv`; `m(a)=E[f|a]`; marginal `p(a)`; `G=−log p`, `G″=m′`.

**Ward (IBP). [T]** Compact fiber ⇒ no boundary; for any `C¹ g`, `E[∂_a g]=E[g f]`. `g=a`: `E[af]=1`. Since `mp=−p′`, this is `∫a(−p′)da=∫p=1` — **normalisation only**.

**Stein. [T]** `E[HessV]=E[∇V⊗∇V]⪰0` (the RHS is a second moment). Coarse direction: `E[H]=E[f²]`. Conditioned: differentiating `m(a)=E[f|a]` under the conditional measure gives
```
E[H|a] = Var(f|a) + m′(a)  =:  Var(f|a) + G″(a).          (Stein-cond)
```
Integrate over `p`: `E[H]=E[Var(f|a)]+E[m′]` and `E[f²]=E[Var(f|a)]+Var(m)`, so `E[m′]=Var(m)` — the **Fisher identity** `∫(p′)²/p=∫(p′)²/p`, a tautology.

**Variance. [T]** `G″(a)=m′(a)=E[H|a]−Var(f|a)`; at `a=0`, `G″(0)=E[H|0]−Var(f|0)=−p″(0)/p(0)`.

**Collapse to one identity. [T]** (`K0-WARD-LOWER §2`) The whole Ward–SD tower `E[a^{2j}(H−f²)]=−2j(2j−1)E[a^{2j−2}]` (`j≥0`) forces the **pointwise** `E[H−f²|a]=−p″/p=G̃(a)`, i.e. `G″(0)=G̃(0)`. Every one of Ward / Stein / variance is a single integration by parts — they are **one** statement.

---

## §2 · Combining into `R`; is it sign-definite for gapped only?

**The relation `R`. [T]** With reflection symmetry (`a,f→−a,−f`, `p` even, `p′(0)=0`, `E[a]=0`), Ward `E[af]=1` plus the Gaussian comparison `κ·Var(a)=1 ⇔ p` Gaussian gives
```
  R :   κ · Var(a) = 1 − η,   κ = G″(0),   η = 1 − κ·Var(a)  (≈0 ⇔ Gaussian).      (R)
```
`R` is *exact* for every measure; `η` is the genuine non-Gaussianity defect. Equivalently `E[H−f²|a]=−p″/p` is the distributional content, and `κ=−p″(0)/p(0)` is read off.

**`R` ⟺ "`p` is a normalised pdf" — `G″` free. [T; V-here numerics]** The three identities impose *nothing* on `p` beyond `∫p=1`:
- Ward ⟺ `∫p=1` (normalisation, vacuous);
- Stein ⟺ `∫(p′)²/p=∫(p′)²/p` (Fisher = Gram, an identity);
- variance ⟺ the *definition* of `G″`.
Hence `p″(0)/p(0)`, i.e. `G″(0)`, is unconstrained. This is verified to quadrature precision (`w7/ward_stein_witness.py`, Part A):

| marginal `p(a)` | `G″(0)` | Ward | Stein residual |
|---|---|---|---|
| Gaussian `χ=0.25` | **+4.000** | 1.000000 | −4·10⁻¹⁶ |
| Gaussian `χ=4` | **+0.250** | 1.000000 | −3·10⁻¹⁷ |
| Gaussian `χ=10⁴` (gapless free limit) | **+1.0·10⁻⁴** | 1.000000 | −1·10⁻²⁰ |
| Bimodal `μ=1.5,σ=0.5` | **−32.0** | 1.000000 | −9·10⁻¹⁶ |
| Bimodal `μ=2.0,σ=0.35` | **−258.4** | 1.000000 | −2·10⁻¹⁵ |
| Flat-top mixture | **+2.78** | 0.999996 | 0 |

**Every** marginal satisfies Ward + Stein to `10⁻¹⁵`; `G″(0)` ranges from `−258` to `+4` to `→0`. (The bimodal formula is closed-form: `G″(0)=(σ²−μ²)/σ⁴`, negative for `μ>σ` — a **dip** at the mode.) **There is no sign-definite combination `≥c>0` only for gapped measures** — a sign-definite `R` would need a constraint the identities do not carry.

**Steelman of the nonabelian pin — and why it is not enough. [T/C/O]** The nonabelian structure supplies exactly one *positive, acid-clean* input the flat free field lacks: the compact-group Haar Jacobian satisfies `Φ″⪰c_Φ` (`c_Φ=N/12`: SU2 1/6, SU3 1/4), so `−log J` is convex and `W_geo:=E[Φ″|0]−Var(Φ′|0)≥0` on the geodesic ball `B_{π/2}` **[T]** (`PROVE-JACOBIAN-SAVES`, `PROVE-Gpp §2`). This is the **only** sign source in the corpus that the flat/abelian free field (`Φ=0`) violates — i.e. it is **acid-clean**. But it cannot be promoted to a `G″` floor, for three exact reasons:
1. **It bounds `W_geo`, not `G″`.** `G″ = brX + W_geo` with `brX:=bracket+X`, `bracket=E[H|0]−Var(G|0)` the action part, `X=−2Cov(G,Φ′|0)`. `W_geo≥0` gives `G″≥brX`; flooring `G″` needs `brX ≥ λ−W_geo` with `brX` **Savvidy-indefinite**.
2. **The variance term subtracts** (see `IL2`). `W_geo=E[Φ″|0]−Var(Φ′|0)`; the `−Var` erodes even the Jacobian floor from `N/12` to `≥0`.
3. **At `β=2.5` the data confirm the action part is the problem:** SU(3) `β=2.5`, `G″=+0.51±0.08` with `E[H]−Var=−0.03<0` (`W7` Part B): the positivity is `≈100%` **Haar-carried**, and the action bracket is at/below zero. So the Haar floor exists but the *action* does not cooperate.

---

## §3 · Gapless-fixed-point test (the free field satisfies `R`)

**Exhibit: the free massless field satisfies Ward + Stein + variance with `G″→0`. [T; V-here]** Take the coarse marginal Gaussian, `p_χ(a)=N(0,χ)`. Then `m(a)=a/χ`, `G″=1/χ`, `E[H|0]=1/χ`, `Var(f|0)=0` (no fibre); Ward `=χ/χ=1`; Stein `∫m′p=∫m²p=1/χ`. The gapless limit is `χ=Var(a)→∞` (the `k→0` zero mode / thermodynamic limit), in which `G″=1/χ→0`. Hence the **free gapless field is a solution of `R`** with `G″=0` — it is **not excluded**. `R` therefore has a gapless fixed point for every value of `G″` down to `0` (and, by the bimodal witness, below `0`).

**Does the nonabelian Ward exclude it? [T: no / partial]** No: the free **nonabelian** Gaussian field (the `g=0` limit of SU(N) YM on the algebra) satisfies every Ward/Stein/BRST identity — they are IBP, abelian/nonabelian-agnostic. The **only** thing that excludes the *flat* free field is the added **inequality** `Φ″⪰N/12` of the compact Haar measure (a property, not an identity); and that excludes the flat field, not the compact gapless class, and does not floor `G″` (§2). This is the exact content of `META-GAPFORCING §3.1` for this triple: an identity's solution set contains the gapless free field.

---

## §4 · Uniform attainment at `β≈2.5`

**The `k→0` target, in `R`-coordinates. [T]** `κ≥λ>0` uniformly `⇔ (1−η)/Var(a) ≥ λ ⇔ Var(a) ≤ (1−η)/λ`. The free-field dispersion gives `4 sin²(k/2) → 0` as `k→0` ⇒ `Var(a;k→0) → ∞` ⇒ `κ→0`; so a **uniform** `κ≥λ` is **exactly** a uniform **upper bound on the low-mode amplitude** `Var(a_{k→0}) ≤ 1/λ`, i.e. a transverse induced mass `m_eff²≥λ` (`PROVE-Gpp §5`, `PROOF-COARSE-POINCARE`). This *is* `H1`.

**At `β≈2.5` the floor is not attained rigorously. [V-here, real data]** From `h1_right` (SU(2), `L=8`, `k=π/2`; `w7/ward_stein_witness.py` Part C) and `W7`:

| `β` | `G″` measured | `±SE` | `1/Var(a)` | `LB_rig` | `G″/free` |
|---|---|---|---|---|---|
| 1.0 | +0.372 | 0.070 | 0.256 | **−0.021** | 0.63 |
| 1.5 | +0.106 | 0.066 | 0.261 | **−0.017** | 0.18 |
| 2.0 | +0.320 | 0.080 | 0.266 | **−0.038** | 0.55 |
| 2.2 | +0.223 | 0.078 | 0.281 | **−0.145** | 0.38 |
| 2.4 | +0.298 | 0.067 | 0.287 | **−0.213** | 0.51 |
| 2.5 (W7 SU2) | +0.274 | 0.004 | 0.27 | — | — |
| 2.5 (W7 SU3) | +0.510 | 0.080 | — | — | (`E[H]−Var=−0.03`) |

`G″>0` **measured** everywhere, but the **rigorous-style lower bound** `LB_rig` (keeping only the proven `W_geo≥0`, bracketing the rest statistically) is **negative at every `β`** (`−0.02…−0.21`). The free reference at the same `k` (`4 sin²(k/2)=0.586`) is **not** uniformly exceeded (`G″/free<1`). So uniform attainment at `β≈2.5` has **no rigorous `G″≥c>0`**; the measured positivity is near-critical and Haar-carried. Exact missing estimate: `(★)` `E[H|0] ≥ λ + Var(f|0)` uniformly in `(V,a,k→0)` — `H1`.

---

## §5 · Acid test + three iron laws

**Acid test. [T; V-here]** Run the triple on the **free massless field** (`χ→∞`): `G″=0`, Ward `=1`, Stein residual `→0` (Part A, `χ=10⁴` row). The identities **do not fail** for the gapless free field — they are **satisfied with `G″=0`**. So the *bare* Ward+Stein+variance is **acid-vacuous as a gap-forcer** (it certifies `0`). The **only** acid-clean piece is the added **`Φ″⪰N/12`** Haar inequality, which the flat free field violates — but it floors `W_geo≥0`, not `G″` (§2), and its promotion to `G″≥c` at `β=2.5` is `(★)=H1`.

**Iron laws (three).**
- **IL1 — Identity is READ, not FORCE. [T]** Any exact identity (Ward/Stein/variance or any combination) holds for **every admissible measure**, including the gapless free field; its solution set contains the gapless point `G″=0`; hence **no sign-definite `≥c>0` consequence** for gapped-only exists. (`META-GAPFORCING §3.1`; here verified: Part A, `G″` free across `±∞`.)
- **IL2 — Variance subtracts (wrong direction). [T]** `G″(a)=E[H|a]−Var(f|a) ≤ E[H|a]`; `Var(f|a)≥0` is a **subtraction**. Poincaré/Brascamp–Lieb/Stein **upper** bounds on `Var(f|0)` can only **cap** `G″`, never floor it; flooring `G″` requires a **lower** bound on the *conditioned* Hessian `E[H|0]`, precisely where Savvidy is indefinite (`λ_min=−g|B|<0`). Witness: 2D `B1` (`E[H|0]=3`, `Var=0.801`, `G″=2.199` — identity exact); `B2` (`E[H|0]=−1.6`, `Var=0.102`, `G″=−1.702<0`).
- **IL3 — The nonabelian sign is Haar-only and `β=0`. [T/C/O]** The unique positive, acid-clean, nonabelian input is `Φ″⪰N/12` (Haar Jacobian); it gives `W_geo≥0` on `B_{π/2}` and the `β=0` gap `≥N/12` **[T]**, but bounds `W_geo` not `G″`, is uncontrolled across the cut `|x|>π/2`, and at `β=2.5` is defeated by the Savvidy-indefinite action bracket (`brX`) — promotion `= (★) = H1` **[O]**.

---

## §6 · Verdict

**Exact failing step (not constructed).**
- The Ward + Stein + variance system at `k→0` **is one identity** `R: κ·Var(a)=1−η`, equivalent to "`p` is a normalised pdf." `G″(0)=−p″(0)/p(0)` is **free**.
- **Gapless fixed point exists:** the free massless field satisfies `R` with `G″→0`; a bimodal marginal with `G″<0`. **No sign-definite `≥c>0` only for gapped measures** is extractable.
- **Uniform attainment at `β≈2.5` fails rigorously:** `G″>0` measured, `LB_rig<0` at every `β`; the free reference is not uniformly beaten.
- **Acid:** the bare triple is satisfied (`G″=0`) by the gapless free field (acid-**vacuous**); the only acid-clean piece (`Φ″⪰N/12`) floors `W_geo`, not `G″`.
- **Reduction:** the single missing estimate is `(★)` `E[H|0] ≥ λ + Var(f|0)` — a lower bound on the *conditioned action Hessian* at the mode, Savvidy-indefinite at `β=2.5`. This is the **local face of `H1`** (`K0-WARD-LOWER §6`, `PROVE-Gpp §4`, `PROOF-COARSE-POINCARE`).

**Sharpest line.** The Ward/Stein/variance route does not build a 6th species; it is the cleanest possible statement that **the whole gap is one conditional comparison** — *the action Hessian at the mode must exceed the force variance by `λ`* — and the identities themselves are exactly the statement "the coarse marginal is a pdf," which the gapless Gaussian also satisfies. The nonabelian Haar curvature is the sole sign carrier the flat free field lacks, but it lives in `W_geo`, not `G″`, and the variance term subtracts: the transfer of its `N/12` positivity onto `G″` across the Savvidy crossover is the machine that does not exist. **`R` = READ, not FORCE; reduces to `H1`; Clay 0%; no QED.**

---

## §7 · Witness (reproducible)

`research/unconditional/w7/ward_stein_witness.py` (`py -3.12`, ~3 s, light: 1D/2D `numpy` quadrature + one JSON read; single job). Output `w7/ward_stein_witness.log`.
- **Part A** — 1D marginals (Gaussian `χ∈{0.25,1,4,10²,10⁴}`, bimodal, flat-top): Ward`=1`, Stein residual `<10⁻¹⁵` for **every** `p`; `G″(0)` spans `+4, +0.25, +10⁻⁴→0, −32, −258`.
- **Part B** — 2D (coarse `a` + fibre `y`): identity `G″(0)=E[H|0]−Var(f|0)` exact to `10⁻⁵`; `B1` (`E[H|0]=3>0`) gives `G″=2.199` with `Var=0.801` subtracted; `B2` (Savvidy-like `E[H|0]=−1.6`) gives `G″=−1.702<0`; both satisfy `G″≤E[H|0]` and Ward`=1`.
- **Part C** — real SU(2) `β∈[1.0,2.4]`, `L=8`: `G″>0` measured, `LB_rig<0` at every `β`; `G″/free<1`.

---

## §8 · References

**Internal (absolute paths; banked labels carried).**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall = Savvidy concavity / key = `G″(0)=H1`); `docs\PROOF-SKELETON.md`.
- `research\unconditional\PROFILE.md` (the one open slot = `H1`); `PROFILE-FACET-ANCHOR.md` (Theorem SA: every structural relation admits the gapless solution; Fredholm zero READ not FORCE); `META-GAPFORCING-IDENTITY.md` (§3.1 identity-forcing no-go; READ≠FORCE; §3.3 free-field compatibility).
- `research\unconditional\K0-WARD-LOWER.md` (Ward–SD tower collapses to `G̃=−p″/p`; `G″(0)=−p″(0)/p(0)`; sign uncontrolled; squeeze rename); `MS-GAP-EQUATION.md` (Ward/Stein gap equation `κ=F(κ)=(1−η)/Var(a)`; gapless attractor; `F′(0)>1` pitchfork absent); `INVENT-VARIANCE-INEQ.md` (`G″=S−D`; SCBL; `D=0⇔Gaussian`, free `G″=0`; `β=0` `≥N/12`); `PROVE-Gpp-positive.md` (Lemma A variance bound, Lemma B `W≥0`, `(★)` the one failing estimate; `G″≈1/Var`); `PROOF-JACOBIAN-SAVES.md` (`Φ″≥N/12`, `G″` decomposition); `PROOF-COARSE-POINCARE.md` (`H1`=Poincaré); `H1-IS-IT-THE-GAP.md`; `PROOF-IR-H1.md` (§4 fast-mode cannot gap slow-mode); `PROOF-SAVVIDY-EXACT-SIGN.md` (`λ_min=−g|B|<0`).
- Numerics: `research\unconditional\w7\ward_stein_witness.py` + `.log` (this work); `h1_right_L8_N200.json`, `h1_right_L8_N200.log` (SU(2) `G″` decomposition, `LB_rig<0`); `w7/k0_continuum_analyze.py` (`G″=bracket+W_geo+X`, `inv=1/Var`); `W7-Gpp0.md` (SU2 `β=2.5` +0.274, SU3 `β=2.5` +0.51 with `E[H]−Var=−0.03`); `FINITE-VOL-CERTIFY.md` (fixed-`(a,L)` free box gap `→0` as `L↑`; finite-vol ≠ Clay).

**External `[R]` (standard).** Stein (1972) / Schwinger–Dyson IBP; Brascamp–Lieb (1976) / Prékopa (1971) / Bakry–Émery (1985); Källén (1952) / Lehmann (1954) (KL spectral measure); Osterwalder–Schrader (1973) / Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Savvidy, *Phys. Lett. B* **71** (1977) 133; Shen–Zhu–Zhu, *CMP* **400** (2023) 805.

---

*Status: **exact failing step.** Ward + Stein + variance at `k→0` collapse to one identity `R: κ·Var(a)=1−η` = "the coarse marginal `p(a)` is a normalised pdf." `R` has a gapless solution (free field, `G″→0`) and admits `G″(0)<0` (bimodal); `G″(0)=−p″(0)/p(0)` is free (`[T]+[V-here]`, Part A). The variance identity subtracts (`G″≤E[H|0]`, IL2), so flooring `G″` needs a lower bound on the Savvidy-indefinite conditioned Hessian. The only nonabelian sign source (Haar `Φ″⪰N/12`) floors `W_geo≥0`, not `G″`; its promotion at `β≈2.5` is `(★)`, the local face of `H1`. Uniform attainment: `G″>0` measured but `LB_rig<0` at every `β`. **READ, not FORCE; reduces to `H1`. Clay 0%. No QED.***
