# SLOT-SD-FIXPOINT — the variance identity `G″=S−D` as a self-consistency: **gapless is a fixed point, `𝒟` is not a function of `S`, and `W_geo` vanishes on the acid field** — the residual is exactly `H1`

> **VERDICT (no QED).** The exact identity `G″(0)=S−D` (`PROOF-JACOBIAN-SAVES`, `INVENT-VARIANCE-INEQ`) **cannot be turned into a self-consistency that excludes the gapless fixed point.** Three independent exact obstacles, any one sufficient: **(O1)** the DSE-closed fixpoint `κ=S(κ)` (defect set to 0 = Gaussian fiber) always has `κ=0` as a solution (free field: `S→0` as `Var(a_{k→0})→∞`); a positive branch exists **iff `F′(0)>1`**, a `k→0` pitchfork = Savvidy tachyon, **absent at β=2.5** `[NG at crossover]`. **(O2)** the exact defect `D` is **not a functional of `S`**: BL gives only the one-sided `K⪰0 ⇒ D≤0`, and even an indefinite fiber does *not* force `D>0` (toy: weight can sit where the quartic makes `K>0`) — `D` is a *measure-weighted overlap* of the cross-Hessian with `K^{−1}`, an independent object; pinning `D<S` is `H1` verbatim `[V-here]`. **(O3)** `W_geo≥0` (Prékopa; `Φ″≥N/12` pointwise) **is** sign-definite but is **zero on the free-Gaussian acid field** (a Gaussian on the algebra carries no Haar Jacobian `Φ`), so it does **not** repel the free-field zero; it is positive only on Haar-carrying measures (β=0, gapped already; β=2.5, where `D<S` is `H1`). **Self-consistent solution with `G″>0` at β=2.5 exists (measured `+0.28`) but is not forced — it coexists with the gapless solution; exclusivity/uniformity is `H1`.** **Failing step pinned, not constructed. Clay 0%.**

---

## §0 · Labels, the three iron laws, the acid test

Labels: **[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-here]** new this round · **[R]** inherited/banked.

**Three iron laws** (`KEY-DIRECT-DERIVATION §0`, carried): **(IL-i) non-positivity** — must fail on the free massless field (`f₁²=0`, `G″=0` there); any Schur-concave positivity only *upper*-bounds the edge. **(IL-ii) nonabelian** — must use the interaction (U(1) is massless). **(IL-iii) no small parameter** — must hold at `α_s~O(1)`, β≈2.5.
**Acid test** (`SRF-FREE-FIELD-ACID`, `KEY-DIRECT-DERIVATION §0`): `Φ(free massless)` must be `≤0` (here `=0`), and *discriminatingly* `≤0` on the gapped free massive `:φ²:` composite (`f₁²=0`).

---

## §1 · The self-consistency `D = 𝒟[S,μ]`, and what makes `𝒟` sign-definite

**The exact identity `[T]`** (`INVENT-VARIANCE-INEQ §1`, `PROOF-JACOBIAN-SAVES §1`). Flat chart `X=log U∈su(N)`, `μ∝e^{−V}dX`, `V=S_W+Φ`, `Φ=−log J`. Coarse mode `v` (`|v|=1`), `a=⟨v,X⟩`; condition on `a=0`; `h_n=vᵀHess V v`, `f=∂_nV=⟨v,∇V⟩`, `Π₀=I−vvᵀ`, `c=Π₀ Hess V v`, `K=Π₀ Hess V Π₀`. Then exactly
```
G″(0) = E_μ[h_n|0] − Var_μ(f|0) = S − D,
S = E_μ[h_n − ⟨K^{−1}c,c⟩|0],   D = Var_μ(f|0) − E_μ[⟨K^{−1}c,c⟩|0].
```
`S` is the Gaussian–Schur curvature (`S=1/Var(a)` exactly for a Gaussian); `D` is the Brascamp–Lieb defect of the conditioned fiber `μ₀`. **This is a determination, not a fixpoint:** at fixed `β`, `S` and `D` are two independent numbers.

**Making it a fixpoint `[T-structure]`.** A genuine self-consistency must let the measure depend on the coarse mass `κ:=G″(0)`. The only intrinsic closure is the DSE/FRG one (`MS-GAP-EQUATION`): replace the fiber by a Gaussian of mass `κ`, so `D_closed(κ)≡0` and `S_closed(κ)=1/Var(a;κ)`, `Var(a;κ)=Σ_k|w(k)|²/(k̂²+Π(k;κ))`. The fixpoint is
```
κ = S_closed(κ) = (1−η(κ)) / Var(a;κ),    η = non-Gaussianity defect  (Ward: κ·Var=1−η).
```
**Identifiability of `𝒟` `[V-here]`.** To keep the exact `D>0` one needs `D=𝒟[S,μ]`. But BL (`Brascamp–Lieb 1976`) supplies only the **implication** `K⪰0 (fiber log-concave) ⇒ D≤0`, with no functional form for `D=𝒟(S)`. Numerically even the implication's converse fails: the toy `V=½Aa²+½By²+(κ/4)y⁴+h a y²` (§4) gives `D=+0.153` at `B=−1` but `D=−0.572` at `B=−2` — **both indefinite fibers**, the sign set by whether the negative region *carries measure weight*, not by `sign(K)` or by `S`. So **`𝒟` is a measure-weighted overlap functional, independent of `S`; there is no intrinsic `D=𝒟[S]`.** `[V-here]`

**Is there a self-consistent `G″>0` at β=2.5?** Yes — measured `G″=+0.281` (§3). But (O1) says it is the *non-forced* branch: the gapless branch is always a solution, and the positive branch is selected only by the external measure (β), not by an intrinsic relation. `[T-cond]+[O]`

---

## §2 · Gapless-fixed-point test: does `G″=0` solve it? Is `W_geo` a sign-definite source at `G″=0`?

**Does `G″=0` solve the self-consistency? Yes, exactly `[T]`.** Free massless Gaussian field on the algebra: `D=0` (Gaussian fiber), and the `k→0` mode has `Var(a)→∞` so `S=1/Var→0`. Hence `G″=S−D=0−0=0`. Acid test: massless `G″=0`, massive `G″=m²` (§4B: `L=64`, `m²=0` → `S=0.0096→0`; `m²=0.5` → `S→0.51`). In the DSE fixpoint this is stronger: `F(0)=0`, so `κ=0` is **the contraction attractor** (`MS-GAP-EQUATION §2`, `BANACH-FIXED-POINT §2`, `[T]`). **The zero is a genuine fixed point, not a removable artifact.**

**Is `W_geo≥0` a sign-definite source at `G″=0` that repels the zero? No `[V-here]` + [T].** `W_geo=E[Φ″]−Var(Φ′)≥0` by Prékopa (`PROOF-JACOBIAN-SAVES §2`), with `Φ″≥N/12` **pointwise** (Lemma B). But this is a property of the *Haar-carrying* measure. On the **acid-test free field** the measure is the pure Gaussian `e^{−½⟨X,AX⟩}` on the algebra — it carries **no `Φ`** (the Haar Jacobian is absent in the linearization; `Φ→0`, and `W_geo` needs `Φ`). Hence
```
W_geo(free Gaussian) = 0   EXACTLY   ⇒   no source at the gapless fixed point.
```
So `W_geo` cannot repel the zero **at the zero**. It is positive only where the measure carries Haar weight:
- **β=0 (pure Haar, strong coupling):** `W_geo=Φ″≥N/12>0`, and `Var(a)` is **L-independent** (compact Haar, §4C: `Var(x₁)=3.72`, `G″→1/Var≈0.27` for any number of links). This is a genuine L-uniform source — **but β=0 is the strong-coupling end, already gapped; it is not the crossing** `[T]`.
- **β=2.5:** `W_geo=+0.22>0` numerically (§3), but `G″=S−D>0` needs `W_geo` (inside `S`) to **beat the positive `D`** — i.e. `D<S` = `H1` `[O]`.

**Conclusion of §2 `[V-here]`.** `G″=0` solves the self-consistency and is its attractor; `W_geo` is sign-definite but is **exactly zero on the field that must be repelled**, and positive only on measures where its dominance over `D` is already `H1`. The zero is **not** repelled by an intrinsic source.

---

## §3 · Uniform attainment at β≈2.5

**Measured (this round, `rerun_su2_L4.npz`, `slot_sd_fixpoint.py §D`) `[V-here measurement]`.** Per-mode conditioned identity `G″=(E[H]−VarG)+W`:

| mode | G″(KDE) | 1/Var | E[H]−VarG | E[Φ″] | W | bracket+W |
|---|---|---|---|---|---|---|
| 0 | +0.250 | 0.274 | −0.031 | 0.330 | +0.220 | +0.189 |
| 1 | +0.287 | 0.250 | −0.123 | 0.329 | +0.236 | +0.112 |
| 2 | +0.286 | 0.288 | +0.066 | 0.334 | +0.233 | +0.299 |
| 3 | +0.294 | 0.271 | −0.002 | 0.340 | +0.230 | +0.228 |
| 4 | +0.226 | 0.278 | −0.002 | 0.370 | +0.173 | +0.171 |
| 5 | +0.345 | 0.272 | −0.009 | 0.332 | +0.221 | +0.212 |
| **mean** | **+0.281** | **0.272** | **−0.017** | **0.339** | **+0.219** | **+0.202** |

`G″>0`; the **Haar-Jacobian term `W` carries it** (`E[H]−VarG≈0`, equipartition); volume-stable `L=4→L=6` (`1/Var=0.272` both, W7-RERUN). So **on the lattice the sign is attained at β=2.5.**

**Uniform `(L,a)` attainment — `[O]`.** What is proved is sign-attainment on a *fixed finite lattice*; what the gap needs is `G″(k→0)≥λ` **uniform in `L` and as `a→0`**. The exact missing estimate is (`PROOF-COARSE-POINCARE`, `KEY-DIRECT-DERIVATION §II`):
```
∃λ>0: Var_μ(∂_nV|a=0) < E_μ[h_n|a=0] − λ,  uniform in (L,a)   [⟺ D < S − λ].
```
This is `H1` verbatim. The S−D route **relabels** it (`D<S`); it does not supply the uniformity. The only L-uniform piece is the **pointwise** `Φ″≥N/12` (Lemma B), but it is a *lattice-compactness* curvature; bounding `Var(G|0)` from above (needed to convert `N/12` into `G″>0`) is the same `(◆)` = `H1` (`PROOF-JACOBIAN-SAVES §2`). **The S−D fixpoint is uniformly attained nowhere it is not already `H1`.** `[O]`

---

## §4 · Acid test + iron laws (exact, `slot_sd_fixpoint.py`)

**A. Toy `S−D` (exact 1-D quadrature).** `V=½Aa²+½By²+(κ/4)y⁴+h a y²`; `h_n=A`, `f=Aa+hy²`, `K=B+3κy²`, `c=2hy`. Identity `S−D=E[h_n|0]−Var(f|0)` holds to quadrature (two right columns equal).

| case | Kmin | S | D | G″ | G″ direct |
|---|---|---|---|---|---|
| free massless Gaussian | 1.0 | 0.000 | 0.000 | **0.000** | 0.000 |
| free massive Gaussian | 1.0 | 0.700 | 0.000 | **0.700** | 0.700 |
| gapped coarse + convex fiber | 1.0 | 0.450 | −0.022 | +0.472 | 0.472 |
| indefinite fiber w/ weight (B=−1) | −1.0 | 0.051 | **+0.153** | −0.102 | −0.102 |
| indefinite but weight at K>0 (B=−2) | −2.0 | −0.953 | **−0.572** | −0.381 | −0.381 |

**Acid test `[T]/[V-here]`.** Free Gaussian: `D=0` exactly, `G″=S`; massless→`0`, massive→`m²`. **Pass** (produces no gap on free; discriminates massive). The indefinite-fiber rows show (O2): `D`'s sign is **weight-dependent**, not a function of `sign(K)` or `S`.

**B. Free Gaussian lattice (exact dense).** `S=4sin²(π/L)+m²`: massless `L=4,8,16,32,64 → 2.0,0.586,0.152,0.038,0.0096→0` (no L-uniform free gap); massive `→m²`. `D=0`. `[T]`

**C. iid Haar SU(2), β=0 (pure-Jacobian source).** `Var(x₁)=3.72`; constant-mode curvature `→1/Var≈0.27`, **L-independent** (1,4,16,64,256 links all ≈0.27). `W_geo>0`, L-uniform — at the strong-coupling end. `[V-here]`

**D. β=2.5 real data** — §3 table. `[V-here]`

**Iron laws.**
- **(IL-i) non-positivity:** `G″=S−D` **passes the acid test** — free massless gives `G″=0` (not `>0`); `S=D=0` there, so the route is non-vacuous on free `[T]`. But the sign `S>D` at β=2.5 is a *Schur-concave variance inequality* that the free field **saturates** (`D=S=0`), so it cannot be made *strict* from positivity alone — proving it for YM is `H1` (`KEY-DIRECT-DERIVATION §II`) `[O]`.
- **(IL-ii) nonabelian:** the defect `D>0` is sourced by the **indefinite fiber `K=Π₀HessVΠ₀`** (Savvidy), which exists only for nonabelian `2 ad F_B` (`PROOF-SAVVIDY-EXACT-SIGN`); U(1) has `D=0`. **Pass in form** `[T]`.
- **(IL-iii) no small parameter:** the identity and `W_geo≥0` are nonperturbative; *but* the β=2.5 positive branch is only *measured*, and the proved sign holds only at β=0 and β<1/48 (`PROOF-P3-UNCONDITIONAL`). **Fail at the crossover** `[O]`.

---

## §5 · Fixpoint map and the pitchfork (`slot_sd_fixpoint.py §E`)

Minimal faithful family with `F(0)=0` (gapless always solves) and `F′(0)=A`:
```
F(κ) = A κ / (1 + κ/κ₀)   ⇒   κ=0 always a fixed point;  positive κ* = κ₀(A−1) iff A>1.
```
| A = F′(0) | positive κ* | regime |
|---|---|---|
| 0.5, 0.9, 1.0 | 0 | gapless only / marginal |
| 1.5 | 0.50 | positive branch |
| 3.0 | 2.00 | positive branch |

`F′(0)>1` **is** the pitchfork = a linear IR instability of the gapless vacuum = the `k→0` Savvidy tachyon (`MS-GAP-EQUATION §3`). At β=2.5 the coarse plane-wave block is **positive** (`+0.40β`, 0% negative; negatives 90–97% high-momentum, `W6-BE-at-2.5`) ⇒ `A<1` ⇒ **only the gapless branch is provable** `[T]`(criterion)`+[NG]`(at crossover). Changing the coordinate/transform (`G″`/`f₁²`/Fredholm/DSE) does not change `A`'s value — it is the same `H1` (`PROFILE §3`).

---

## §6 · Verdict — the exact failing step

**Not constructed. The self-consistency admits `G″=0`.** Precisely which line fails (any is sufficient; all are the same wall):

1. **(O1) Gapless solves it.** The DSE-closed fixpoint `κ=S_closed(κ)` has `κ=0` as its attractor for every `L` (free `F(0)=0`); the positive branch needs `F′(0)>1`, absent at β=2.5. `[NG]`
2. **(O2) No `D=𝒟[S]`.** `D` is an independent measure-weighted overlap (toy B=−1 vs B=−2); BL is one-sided. Bounding `D<S` is `H1` verbatim, `[O]`.
3. **(O3) `W_geo` does not repel.** It is `≥0` (Prékopa) but **`=0` on the free-Gaussian acid field** (no `Φ`); positive only on Haar measures where its dominance over `D` is `H1`. `[V-here]`+`[O]`.
4. **(Uniformity) `[O]`.** Even granting `G″>0` at β=2.5 (measured `+0.28`), the `(L,a)`-uniform `G″≥λ` is `H1`; the S−D route relabels, does not supply it.

**Honest status.** In the `S−D` coordinate the whole campaign reduces — again — to `D<S` = `Var(G|0)`-bound = `H1`, with the extra finding that **`W_geo` cannot rescue it because it vanishes exactly on the gapless field**. There is **no** self-consistent `S,D` relation that forces `S−D>0` and excludes `G″=0`. This is the `PROFILE`‑`ANCHOR` verdict reached through the variance coordinate: representational freedom (which decomposition `G″=S−D` / QV-sum) is illusory; the residual is `H1`. **Clay 0%. No QED.**

---

## §7 · Files, reproducibility, refs

**Witness:** `research/unconditional/w7/slot_sd_fixpoint.py` (single process, no heavy MCMC; Part A quadrature, B dense Gaussian, C iid-Haar, D loads `rerun_su2_L4.npz`, E fixpoint). Run: `py -3.12 slot_sd_fixpoint.py`.

**Internal (ground):** `PROOF-JACOBIAN-SAVES.md` (`G″` identity, Prékopa, `Φ″≥N/12`, (◆)); `INVENT-VARIANCE-INEQ.md` (SCBL `G″=S−D`, acid test); `PROFILE-FACET-ANCHOR.md` (§2.3 gap-equation, §2.6 Theorem SA, §3 residual DOF); `PROFILE.md` (single open slot = `H1`); `MS-GAP-EQUATION.md` (`κ=F(κ)`, `F′(0)>1` pitchfork, tachyon UV); `PROOF-COARSE-POINCARE.md` (H1 exact missing estimate); `KEY-DIRECT-DERIVATION.md` (§0 iron laws/acid; §II variance split = H1), `W7-RERUN.md` (`G″` decomposition table), `W6-BE-at-2.5.md` (negatives high-momentum), `PROOF-P3-UNCONDITIONAL.md` (β<1/48), `PROOF-SAVVIDY-EXACT-SIGN.md` (`λ_min=−g|B|`), `SRF-FREE-FIELD-ACID.md` (free-field acid), `BANACH-FIXED-POINT.md` (gapless attractor), `SCHUR-CONVEX.md` (free-field saturation).

**External `[R]`:** Brascamp–Lieb, *Adv. Math.* **14** (1976) 366; Prékopa, *Acta Sci. Math.* **32** (1971) 301; Savvidy, *Phys. Lett. B* **71** (1977) 133; Shen–Zhu–Zhu, *CMP* **400** (2023) 805.

*Status: **SLOT-SD-FIXPOINT attempt — FAILING STEP PINNED, not constructed.** The exact identity `G″=S−D` is `[T]` and acid-clean `[T]`; turned into a self-consistency it admits `G″=0` (O1), `𝒟` is not a function of `S` (O2), and `W_geo` is zero on the gapless field so cannot repel it (O3). Uniform attainment at β=2.5 is `[O]`=`H1`. Witness reproduces all numbers. Not a proof of Clay; not a global impossibility theorem (`∀class ≠ ∀possible`). Clay 0%; no QED.*
