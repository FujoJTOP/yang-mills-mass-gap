# k0-TOPOLOGICAL-COUPLING — does global topology (`Q`, `θ`, `χ_t`) force a variance floor on the `k→0` mode `a_0`?

> 2026-11-23 · FORGE5 probe. Target: `Var_μ(a_0) ≤ 1/λ` (the gap in measure language; gapless ⟺ `Var(a_0)=∞`). Direction tested: the zero-momentum mode couples to the **global** instanton number `Q`, the `θ`-angle, and `χ_t=Var(Q)/V`, and this coupling transfers a bound onto `Var(a_0)`. **判词: 不过 — the arrow runs backward, and at `θ=0` parity zeroes the coupling outright.** Even a *proven floor* on `Var(a_0)` is compatible with gaplessness; the gap needs a *ceiling*. `χ_t>0` is a bulk **softness** (`Var(Q)=χ_t V→∞`), the opposite of the stiffness required. Of the 5: **#1, #4 pass; #2, #3, #5 fail — one object: positivity/moments supply the wrong side.**

## Why a "variance floor" is already the wrong object
Gap ⟺ `Var_μ(a_0) ≤ 1/λ < ∞`. Gapless ⟺ `Var(a_0)=∞`. A floor `Var(a_0) ≥ c>0` is satisfied by **both**; alone it is the signature of a soft mode, not a gap. The brief's own caution is right: `χ_t>0` is positivity of `Var(Q)`, and positivity is a lower/moment statement (`M_k≥0 ⇒ ratios are means ⇒ upper`, GAP-IN11 moment lemma), the wrong side for a ceiling.

## Exact failure I — the Cauchy–Schwarz arrow
`|Cov(a_0,Q)|² ≤ Var(a_0)·Var(Q)` gives **`Var(a_0) ≥ Cov²/Var(Q)`** — a floor. With `⟨Q⟩=0`, `χ_t>0`, `Var(Q)=χ_t V→∞`, so the floor `→0` as `V→∞`. The larger the susceptibility, the weaker the statement; maximal coupling still yields only `Var(a_0) ≥ 0`. A ceiling would need an **upper** bound on `Var(Q)`; `χ_t>0` supplies the opposite.

## Exact failure II — orthogonal decomposition merely diverges
Write `Q=c·a_0+R`, `Cov(a_0,R)=0`: `Var(Q)=c²Var(a_0)+Var(R) ⇒ Var(a_0) ≤ (χ_t V−Var R)/c² → ∞`. Right sign, but `∝V`: no uniform ceiling — and (III) `c=0` at `θ=0`.

## Exact failure III — parity kills the coupling at `θ=0`
`Q=∫tr F F̃` is a `0⁻⁺` pseudoscalar, `P`-/`CP`-odd; the residue `f₁²=⟨0|¼trF²|G⟩²` is `0⁺⁺`, `P`-even. The `θ=0` vacuum is `P`- and `CP`-invariant, so `Cov(a_0,Q)=⟨¼trF²·Q⟩_conn = 0` **exactly** at the Clay point — the bound reads `Var(a_0) ≥ 0`. At `θ≠0` a coupling survives, but that is not the Clay problem, and it is cross-channel: it touches the pseudoscalar (`0⁻⁺`) variance, never the scalar residue (`PROOF-1FORM-ANOMALY`).

## Exact failure IV — sign wall and circularity
`χ_t` is a positivity total (`ν≥0`), `M_{−1}` of the pseudoscalar spectral density: KL gives `χ_t=Σ f_n² m_n² ≥ f₁² m_{0⁻⁺}²`, i.e. `m_{0⁻⁺}² ≤ χ_t/f₁²` — an **upper** bound in its *own* channel (SCHUR-CONVEX; IN-11). Same `ν≥0` sign wall as B1. Non-circular? `χ_t>0` is **unproven on ℝ⁴** (instanton gas IR-divergent, non-dilute, B1/W2); assuming it presumes a dilute liquid with a scale `Λ` — smuggling the gap. Witten–Veneziano `χ_t∝(f_π²/2N_f)m_{η'}²` is the η'/`0⁻⁺` channel and reads `m²` *from* `χ_t` (upper relation), never a `0⁺⁺` floor.

## Five-point verdict
1. **Nonperturbative — PASS.** `χ_t→0` at `g=0` (`e^{−8π²/g²}`); genuinely evades EDGE's `g=0` no-go.
2. **Sign-controlling — FAIL.** Positivity → upper/moment; C-S → a floor that `→0`. No lower bound on `G″(0)`.
3. **Non-circular — FAIL.** `χ_t>0` unproven; dilute-liquid control or Witten–Veneziano imports the scale / the wrong channel.
4. **Non-SUSY + fails free field — PASS (narrow).** `χ_t=0`, `Var(a_0)=∞`; acid test `0≥0` clean — necessary, not sufficient.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Every transfer scales `∝χ_t V` (diverges) or `→0`; no uniform-in-`(V,a)` floor.

**Sharpest line.** Topology hands a **softness** — `χ_t>0`, `Var(Q)=χ_t V` unbounded, the hallmark of a delocalized (axion-like) global mode — exactly where the gap needs a **stiffness** `Var(a_0)≤1/λ`. Cauchy–Schwarz is monotone the wrong way, the orthogonal route diverges like `V`, and at `θ=0` parity zeroes the `a_0–Q` coupling outright. A variance floor is compatible with the very gaplessness we must exclude; it cannot certify the gap. **#2, #3, #5 fail (one object: positivity gives the wrong side); #1, #4 pass.**

*Refs.* B1-nonperturbative-residue §2 (`χ_t=M_{−1}`, upper, wrong channel/sign, `χ_t>0` unproven); GAP-IN11 (moment lemma, always upper); SCHUR-CONVEX (`χ_t≥f₁²m² ⇒ f₁²≤χ_t/m²`, wrong direction); IN-11-attack / BREAK-IN-11 (`m_{0⁻⁺}²≤χ_t/f₁²`); LEVER-nonconvex-edge (linear data never lower-bounds the edge); THETA-PI (`E(θ)` pins vacuum energy, not `m_{0⁺⁺}`; `θ=π≠θ=0`); PROOF-1FORM-ANOMALY (local `trF²` neutral, cross-channel); INVENT-VARIANCE-INEQ (`G″=S−D`, `S=1/Var`; gap ⟺ `G″` lower ⟺ `Var` upper); PROOF-IR-H1 (H1 = `Var(a)≤1/λ`); EDGE-f1sq-attack (free-field no-go). Ext.: Witten, Nucl. Phys. **B156** (1979) 269; Veneziano, Nucl. Phys. **B159** (1979) 213; Callan–Dashen–Gross, Phys. Rev. D **20** (1979) 3279. No fabricated lemma; this is a confirmed dead end (B1 §2) re-derived for the `k→0` target.
