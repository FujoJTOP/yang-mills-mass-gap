# BASE-FIELDMEASURE — the measure-inequality family is NOT a source: a classification, and the anomaly escape

> **VERDICT (no QED).** Of the natural **field-measure functionals** — FKG defect, log-Sobolev / spectral
> gap, entropy production, Bakry–Émery `Γ₂`, Wasserstein drift, Brascamp–Lieb defect, transfer entropy,
> Gaussian-domination defect — **none** is a κ-independent, sign-definite **source** `𝒮` in `κ=F(κ)+𝒮`.
> They fall into three provably-dead classes: **(A) exactly zero on equilibrium** (entropy production,
> transfer entropy); **(B) proportional to `κ`** — homogeneous, so absorbing them into `F` leaves the zero
> a solution (Poincaré/LSI gap, Bakry–Émery, Wasserstein gap, Gaussian-domination); **(C) sign-indefinite**
> (Brascamp–Lieb defect; FKG defect). **The only inhomogeneous field-measure objects are *anomalies*
> (`β(g)≠0` breaks scale invariance).** In the `0⁺⁺` channel the anomaly object is the **gluon condensate
> `⟨tr F²⟩`**, whose positivity-and-magnitude floor *is* `H1` and whose uniform `a→0` is the renormalon
> wall; the *other* anomaly object, **topological susceptibility `χ_top`**, is genuinely κ-independent and
> nonzero but **sources the `0⁻⁺` channel** (Witten–Veneziano / `η'`), not Clay's `0⁺⁺` — **channel
> mismatch**. **So: the base is classified — a source must be an anomaly, and the unique `0⁺⁺` anomaly is
> `H1`. No construct. Clay 0%.**

---

## §0 · Labels, iron laws, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-here]** new · **[R]** inherited.
**IL-i** non-positivity (fail on free massless: `f₁²=0`, `G″=0`). **IL-ii** nonabelian (U(1) massless). **IL-iii** no small parameter (`α_s~O(1)`, `β≈2.5`).
**Acid (`SRF-FREE-FIELD-ACID`)**: `𝒮(free massless)=0` *and* `𝒮(:φ²:)=0` (the gapped-atomless composite).

**The S-scheme** (from `SLOT-SYNTHESIS` §4). A source `𝒮(μ)` must be
**(S1)** acid-clean on *both* adversaries; **(S2)** `>0` on the Haar-carrying interacting measure at `β≈2.5`;
**(S3)** **inhomogeneous**: `κ=F(κ)+𝒮` with `𝒮` **κ-independent**, so `κ=0` is **not** a solution;
**(S4)** `κ*` bounded below uniformly in `(L,a)` across `[1/48,2.5]`.

---

## §1 · The scorecard — 8 candidates × 4 conditions

| Candidate | S1 acid (free / `:φ²:`) | S2 `>0` at β=2.5? | S3 κ-independent source? | S4 uniform? | dies at |
|---|---|---|---|---|---|
| **FKG / positive-association defect** | `0`/`0` (free Gaussian is FKG) | unknown, **sign-indefinite** | ✗ measure-global coefficient | ? | **S3** |
| **Log-Sobolev / Poincaré gap** | `→0` as `κ→0` (not clean) | yes (but `=κ`) | ✗ **homogeneous** (`=κ`, see §2A) | = `H1` | **S3** |
| **Entropy production** | `0` | **`0` (equilibrium)** | ✗ | — | **S2** |
| **Bakry–Émery `Γ₂`** | `0` (free: BE const `=0`) | yes (`=`inf Hess) | ✗ **homogeneous** (`=κ`) | = `H1` | **S3** |
| **Wasserstein drift** | `0` | yes (but `=`gap) | ✗ **homogeneous** | = `H1` | **S3** |
| **Brascamp–Lieb defect `D`** | `0`/`0` | **sign-indefinite** (§2B) | partial (independent of `S`) | = `H1` | **S2** |
| **Transfer entropy** | `0` | **`0` (equilibrium)** | ✗ | — | **S2** |
| **Gaussian-domination defect** | `0` | yes (but `=`sup Hess) | ✗ **homogeneous** | = `H1` | **S3** |

**0 / 8 meet S3 (the decisive condition).** The "source" is never κ-independent: it is either `0`, or it
*is* `κ`.

---

## §2 · The two decisive computations

### 2A · The Gaussian coarse source is *linear* in `κ` — homogeneous (so not a source) `[T]`

1D periodic free field, mass `m²=κ`. The `k=0` coarse mode `a` has `Var(a)=1/(Lκ)` (it projects onto the
zero mode), so
```
S := 1/Var(a) = L·κ .
```
Witness `w7/base_fieldmeasure.py` Part A (L=64):

| `κ` | `Var(a)` | `S=1/Var` | `S/κ` |
|---|---|---|---|
| 0.5 | 0.031 | 32 | **64** |
| 0.01 | 1.56 | 0.64 | **64** |
| 1e-4 | 156 | 0.0064 | **64** |

`S/κ = L` **exactly constant** ⇒ `S` is a **homogeneous** function of `κ`, `S→0` at `κ=0`. A source in
`κ=F(κ)+S` must survive at `κ=0`; `S` does not. **The Gaussian "source" is not a source — it is the
κ-dependent right-hand side.**

### 2B · The Brascamp–Lieb defect is *sign-indefinite* `[V-num]`

2-mode toy `V=½Aa²+½By²+(κ/4)y⁴+h a y²`, conditioned at `a=0`; `D=Var(f|0)−E[c²/K|0]`, `f=h y²`,
`c=2hy`, `K=B+3κy²` (witness Part B, `A=1,h=0.5,κ=1`, vary fiber curvature `B`):

| `B` | −3 | −2 | −1 | −0.5 | 0 | +1 | +3 |
|---|---|---|---|---|---|---|---|
| `D` | −0.094 | **+0.223** | −0.457 | **+0.236** | −0.198 | −0.062 | −0.025 |

`D` **changes sign** with `B` (the sign is set by whether the negative-curvature region carries measure
weight, not by `sign(K)` — matching `SLOT-SD-FIXPOINT` O2). **Not sign-definite ⇒ cannot be a
sign-definite source.** (For `B<0` the fiber is indefinite, `K` changes sign, and `E[c²/K]` can itself be
negative — `−0.055` at `B=−0.5`.)

### 2C · Equilibrium functionals vanish identically `[T]`

Entropy production `σ=∫J·∇log(p_ss/p) = 0` for **any** equilibrium measure (detailed balance ⇒ probability
current `J=0`); transfer entropy `T_{past→future}=0` at equilibrium (time reversal). Both are `0` on every
YM measure, hence cannot be `>0` at `β=2.5`.

---

## §3 · The classification (the actual "base")

> **Theorem BE (class no-go) `[NG]` for the measure-inequality family.** Any field-measure functional
> built from the Markov-semigroup / gradient-flow structure of `μ` — i.e. a functional of the generator
> `L=−∇V·∇+Δ` defined through a *homogeneous* inequality (Poincaré/LSI/Bakry–Émery/Wasserstein/FKG/
> Gaussian-domination), or an information-production rate (entropy production / transfer entropy) — is,
> respectively, **homogeneous of degree 1 in the Hessian** (class B: equals or is bounded by a multiple of
> `κ`, vanishes at `κ=0`) or **identically zero at equilibrium** (class A) or **sign-indefinite** (class C,
> Brascamp–Lieb). Therefore **it cannot appear as the κ-independent additive term `𝒮` in `κ=F(κ)+𝒮`.**

**Proof sketch.** Each defining inequality is invariant under `V→λV` (for class B) so the associated
constant scales with the Hessian and tends to `0` at the free massless point; class A vanish by detailed
balance; class C is sign-indefinite by §2B. `∎`

The classification says **where a source can and cannot live**: *not* in the measure-inequality family —
those are the outcome (`κ`) disguised, or zero, or indefinite.

---

## §4 · The escape: only an anomaly is inhomogeneous

The one structure in the field measure that is **not** homogeneous is **scale breaking**, `β(g)≠0`. The
trace anomaly
```
θ^μ_μ = (β(g)/2g) tr F² ,   ⟨θ^μ_μ⟩ = −4 ε_vac
```
introduces a scale `Λ` that is **independent of the coarse mass `κ`** — the only candidate that can be a
κ-independent source. Two anomaly objects exist, and both are blocked, in *different* ways:

- **Topological susceptibility `χ_top = ⟨Q²⟩/V`.** Genuinely κ-independent (topological, `Q∈Z`), nonzero
  for interacting YM (`≈(180 MeV)⁴`), `=0` for free field and `:φ²:`. It satisfies **S1, S2, S3** in
  principle — **but it sources the `0⁻⁺` (pseudoscalar, Witten–Veneziano / `η'`) channel via the anomaly,
  not the `0⁺⁺` scalar Clay needs.** *Channel mismatch* (the same no-go as `DEFORM-DECOUPLE`/defect
  channels): a nonzero `χ_top` does not floor the `tr F²` (=`0⁺⁺`) mass.
- **Gluon condensate `⟨tr F²⟩`** (`=⟨θ⟩`, the `0⁺⁺` anomaly). Right channel, κ-independent, acid-clean
  after normal ordering (`:trF²:` = 0 on the free measure). **But** proving `⟨tr F²⟩` is nonzero **with a
  magnitude floored uniformly in `(L,a)`** *is exactly* `H1`; and separating the physical condensate from
  the `a⁻⁴` UV divergence uniformly is the **renormalon wall** (`S4`).

**So the `0⁺⁺` source exists in principle — and it is the gluon condensate, which is `H1`.**

---

## §5 · Honest bottom line

- **No construct.** No measure-inequality functional is a κ-independent sign-definite source (§1–§3).
  The closest object to S1–S3 is the **topological susceptibility** — but the wrong channel — and the only
  right-channel object is the **gluon condensate**, whose S2 is `H1` and whose S4 is the renormalon wall.
- **The base is laid:** the classification `[NG]` rules out the entire measure-inequality family in one
  theorem, and localizes any future source to **anomalies in the `0⁺⁺` channel** — a sharp, reusable
  constraint (do not propose another Poincaré/LSI/FKG/Wasserstein functional as "the source"; it is
  homogeneous by §3).
- **What would fill it:** a proof that the `0⁺⁺` trace anomaly `⟨tr F²⟩` is floored uniformly in `(L,a)`
  — i.e. `H1` — or a dual/anomaly-inflow statement that transfers `χ_top`'s κ-independence into the
  `0⁺⁺` channel. Neither constructed.
- **Clay = 0%.** No QED.

*Labels:* §1 table and §3 theorem `[NG]` (measure-inequality class); §2A `[T]` (`S=Lκ` exactly on the free
field), §2B `[V-num]`, §2C `[T]`; §4 channel mismatch `[NG]`, condensate `=H1` `[O]`. **No gap-forcing
identity claimed.**

*Refs:* `SLOT-SYNTHESIS.md`, `SLOT-SD-FIXPOINT.md` (BL defect), `PROFILE.md` (P1 field-measure),
`DEFORM-DECOUPLE-4D.md` (channel mismatch), `AREALAW-GAP.md`/renormalon notes.
Witness: `w7/base_fieldmeasure.py` (runs clean, `py -3.12 -I`).
