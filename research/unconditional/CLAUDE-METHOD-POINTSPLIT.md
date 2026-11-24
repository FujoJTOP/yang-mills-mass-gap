# CLAUDE-METHOD-POINTSPLIT — the finite-time (point-split) susceptibility: it **cures the `a⁻⁴` UV degeneration**, but on the hard side it is **H1 in finite form** — the atom wall (`f₁²≥r₀`) is unmoved

> 2026-11-23 · Audit of a proposed fix to the **sufficient condition** `(f₁²≥r₀) ∧ (χ≤1/λ) ⇒ m₁²≥r₀λ`.
> The problem: for the dim-4 operator `tr F²`, `χ=∫d⁴x⟨OO⟩_c ~ c·a⁻⁴` (UV), so `χ≤1/λ` is unsatisfiable and `f₁²/χ→0`. The proposed fix: replace the full susceptibility by the **finite-time** (point-split-in-time) `H(t₀)=∫_{t₀}^∞C(t)dt`, giving `m₁²e^{m₁²t₀}≥f₁²/H(t₀)`.
> Iron rules (charter §2): every claim is labelled **真定理 / 条件定理 / no-go / 精确缺失步**. Nothing unproven is called proved; the free-field acid test is applied to every tool.

**判词 (verdict, up front).**

1. **(a) The UV degeneration is REAL — 真定理 (dimensional analysis / OPE), and stronger than stated.** For any gauge-invariant `O` of scaling dimension `Δ` the Källén–Lehmann spectral density is `ρ(s)~s^{Δ-2}` at large `s` (matched to `⟨O(x)O(0)⟩~|x|^{-2Δ}`), so
   `χ = ∫ s⁻¹dν(s) = ∫ s⁻¹ρ(s)ds ~ Λ_UV^{2Δ-4} = a^{-(2Δ-4)}.`
   For `tr F²` (`Δ=4`) this is `χ ~ a⁻⁴`; for **every** pure-YM gauge-invariant (`Δ≥4`, since `A` is dim 1, `F` dim 2, the lowest trace is `tr F²` dim 4) it diverges as `a^{-(2Δ-4)}≥a⁻⁴`. The divergence lives entirely in the **contact/coincidence region** (`x=0` and `|x|~a`); the IR (physical) part is `O(a⁴)` in lattice units. Corollary — **a 真 no-go for the literal H1 route on `tr F²`**: `χ≤1/λ` is **unsatisfiable**, so `(f₁²≥r₀)∧(χ≤1/λ)` can never hold. There is **no** contact-term subtraction that rescues the *criterion as written*; subtracting the contact gives a finite `χ_R`, which is a *different* (H1-species) premise, not the original.
   *Reconciliation with the measured `χ_lat≈0.04` (H1-NUMERIC-RIGHT §2):* the lattice computes the **dimensionless** `q≈a⁴·tr F²`, so `χ_lat = a⁴ χ_phys = c + O(a⁴)` — the measured `0.04`, `β`-insensitive, **is the UV contact coefficient `c`, not the IR susceptibility.** The IR piece is `O(a⁴)` and invisible. This *confirms* rather than contradicts "χ is gap-blind".

2. **(b) The point-split ESCAPES the degeneration but NOT the wall — it is H1 in finite form on the hard side.** The finite-time `H(t₀)` is **UV-finite for `t₀>0`** (真定理): `H(t₀)=∫dm ρ(m)e^{-m t₀}/(2m²) ~ t₀^{-4}` is cut off at `m~1/t₀`, with no `a⁻⁴`. So it converts an **unsatisfiable** sufficient condition into a **satisfiable** one — a genuine improvement in *well-posedness* (this is the real content of the fix, and it is correct). But the reduced criterion is **(I) atom `f₁²≥r₀`** (or `w₁≥r₀`) **∧ (II) finite-time tail `H(t₀)≤X`** — and:
   - **(II) is STRICTLY WEAKER than H1** (真定理): H1 `⇒` (II) since `∫_{t₀}^∞≤∫_0^∞=χ`; but (II) `⇏` H1 (the short-time part `∫_0^{t₀}C` — the contact — is unconstrained by (II)). So (II) is "H1 with the divergent short-time region excised."
   - **(II) is the SAME SPECIES as H1**: it is a one-sided **upper bound on a positive (completely-monotone, reflection-positive) correlator moment** — the corpus's named species "one-sided variance/upper bound." Positivity (RP, complete monotonicity) supplies **lower** bounds and monotonicity, never the needed tail **upper** bound; that requires Gaussian-domination / spectral-moment control, i.e. the same Brascamp–Lieb / coarse-Poincaré machinery as H1 (and Lemma SC proves exactly this at `β<1/48`).
   - **(I) is unchanged** — it is verbatim `H*(r₀)`, the unprotected residue = the gap.
   **Net: ESCAPES the `a⁻⁴` degeneration; SAME WALL (`I`∧ upper-bound species) on the content.** The fix relocates, it does not remove, the unprotected atom.

3. **(c) The numeric kill-switch is well-posed and normalization-invariant, but it is decisive about the ROUTE's viability, not about the gap, and the current `L≤8` ensembles cannot resolve the atom needed for it.** It cleanly separates "atom share survives as `a→0`" (route alive) from "atom share →0" (route dead); but a surviving `m*>0` is *numerical evidence*, not proof, and the criterion is *sufficient-not-necessary* (failing it does not disprove the gap). Reaching the signal needs **new** `L≥16–24` scaling-window ensembles.

---

## §0 · Definitions — one convention, fixed (the overloaded `a`, `ν`, `m` the whole confusion; the 1/(2m) the audit caught)

A gauge-invariant `0⁺⁺` interpolator `O` (the A1g-projected `tr F²` / plaquette density). Euclidean. Two representations of the same physical spectral function `ρ` — **do not mix them** (this is the audit's genuine catch):

| repr. | correlator | measure | susceptibility | single atom ⇒ |
|---|---|---|---|---|
| **proper-time / momentum** (the paper's, body2) | `C(p)=∫_0^∞ dν(s)/(p²+s)`, `s=m²` | `dν` on `s=m²` | `χ=∫ s⁻¹dν(s)` | `f₁²/χ=m₁²` **exactly** |
| **Euclidean time** (what the lattice measures) | `C_E(t)=∑_n w_n e^{-m_n t}`, `w_n=|⟨0|O|n⟩|²/(2m_n)` | weights `w_n` | `χ=2∫_0^∞C_E(t)dt = ∫s⁻¹dν` | `w₁e^{-m₁t}`, `χ=w₁/m₁` |

- `f₁²:=ν({m₁²})=|⟨0|O|G⟩|²` (the atom weight; dim 8 for `Δ=4`); in Euclidean weights `w₁=f₁²/(2m₁)` (dim of `C_E`).
- `m₁=inf supp ν` — the lightest `0⁺⁺` mass = the **gap**.
- `M₀=∫dν` — total spectral weight (diverges for `Δ=4`; GAP-G2-edge already records this).
- The bridge between the columns is exactly the `1/(2m)` the audit flagged: `C_E(t)=∫dm ρ(m)e^{-mt}/(2m)` while the proper-time object is `∫dν(s)e^{-sτ}`. **`χ` is the same number in both columns**, but the *finite-time integrals differ* (`e^{-mτ}` vs `e^{-m²τ}`). Using `χ=∫s⁻¹dν(s)` to integrate the Euclidean `C_E(t)=∑w e^{-mt}` is the mixed-representation error.

**The squeeze** (真定理, either column, algebra done not quoted):
```
f₁²/χ  ≤  m₁²  ≤  M₀/χ          (proper-time;  f₁²=ν({m₁²}), M₀=∫dν)      [paper's form]
w₁/χ   ≤  m₁   ≤  M_w/χ          (Euclidean;     w₁=|⟨0|O|G⟩|²/(2m₁), M_w=∑w_n)
```
Lower side = the edge-atom term (`χ≥f₁²/m₁²`); upper side = residue-weighted harmonic mean (`m₁²≤M₀/χ`; harmonic mean `≥` min). They coincide iff `ν` is a single atom.

---

## §1 · The UV degeneration — verified, and forced by dimension (真定理 + 真 no-go)

**Claim (audit).** `⟨O(x)O(0)⟩_c ~ c/|x|⁸` ⇒ `χ~a⁻⁴` ⇒ `f₁²/χ ∝ a⁴→0`, normalization-invariant.

**Proof of the scaling (真定理, OPE/dimensional analysis).** `O` has scaling dimension `Δ`. Matching `⟨O(x)O(0)⟩~|x|^{-2Δ}` (dim `2Δ`) to the KL integral `⟨O(x)O(0)⟩=∫_0^∞ds ρ(s)Δ(x;√s)` with `Δ(x;m)` the free propagator gives, at `s|x|²≪1`, `∫^{1/|x|²}ds ρ(s)·(1/|x|²) ~ |x|^{-2Δ}`, i.e. `ρ(s)~s^{Δ-2}`. Then
```
χ = ∫_0^{Λ²} s⁻¹ρ(s)ds ~ ∫^{Λ²} s^{Δ-3}ds ~ Λ^{2Δ-4} = a^{-(2Δ-4)}.          (真定理)
```
- `Δ=4` (`tr F²`, `tr F F̃`): `χ~a⁻⁴`. ✓ (audit correct)
- **Pure YM has no gauge-invariant of dimension `<4`** (`A` dim 1 is not invariant; `F` dim 2 is Lie-algebra-valued; the lowest trace is `tr F²` dim 4). Hence **every** pure-YM gauge-invariant has `χ ≥ c a⁻⁴ →∞`. ✓ Stronger than the audit's single-operator statement.
- The divergence sits in the **contact/coincidence region** (`x=0`, `|x|~a`); the IR piece is `O(a⁴)` in lattice units (see below). 

**Consequence — 真 no-go (the literal H1 route on any pure-YM gauge-invariant).** `χ≤1/λ` with `λ>0` bounded below is **unsatisfiable**. Together with `χ≥f₁²/m₁²` this gives only `m₁² ≥ f₁²/χ → 0` — vacuous. So the sufficient condition `(f₁²≥r₀)∧(χ≤1/λ)` is **never satisfied** for `tr F²`; the literal H1 route is dead. This is a no-go for the **criterion**, not for the physics (the gap `m₁²` is unaffected — the *bound* degenerates).

**Normalization invariance (真).** Under `O→cO`: `f₁²→c²f₁²`, `χ→c²χ`, so `f₁²/χ` is invariant. ✓ The degeneration is not a normalization artifact.

**Contact-term subtraction — what it does and does not fix.** On the lattice, `⟨O(x)O(y)⟩⊃c(a)δ⁴(x−y)+` power-law contact `~|x|^{-8}`. Subtracting the contact (`χ_R:=χ−contact`) gives a **finite** `χ_R` (IR-dominated, `~f₁²/m₁²+`tail). The squeeze `f₁²/χ_R≤m₁²` is then well-posed — **but** `χ_R≤1/λ` is a *one-sided upper bound on the IR susceptibility = H1-species*, not the original H1 (and not proven). So subtraction **moves** the problem from "divergent, false" to "finite, same species, open" — it does **not** escape the wall. *(This is precisely what the point-split implements in the time domain, §2.)*

**Reconciliation with the measured `χ_lat≈0.04`. 真 + consistency check.** The lattice measures the dimensionless `q≈a⁴trF²/4`; then `χ_lat=a⁴χ_phys=c+O(g²)+O(a²)`, where `c` is the OPE contact coefficient. The measured `χ_lat=0.043–0.055` (H1-NUMERIC-RIGHT §2), `β`-insensitive over `β=1.0–2.4`, is **`c`** — the UV contact — not the IR susceptibility. The IR susceptibility is `a⁴·(f₁²/m₁²+…)=O(a⁴)` in lattice units, swamped by `c`. This *explains* why "χ is finite and gap-blind": it is the contact term. It corroborates both the audit and the earlier "susceptibility is the wrong object" verdict.

---

## §2 · The point-split / finite-time susceptibility — derivation, representation fix, acid test

Define the **finite-time susceptibility** (a time cut on the Euclidean-time integral — a "time split", not a spatial point-split; see §4):
```
H_E(t₀) := ∫_{t₀}^∞ C_E(t)dt = ∑_n w_n e^{-m_n t₀}/m_n,       w_n=|⟨0|O|n⟩|²/(2m_n)≥0.    (真定理, KL/OS)
```
*(Proper-time analogue: `H_P(τ₀)=∫_{τ₀}^∞[∫e^{-sτ}dν(s)]dτ=∫s⁻¹e^{-sτ₀}dν(s)`; `H_P(0)=χ`. The two columns differ by the `1/(2m)` and the `e^{-mτ}` vs `e^{-m²τ}` — §0.)*

**Key identity — UV-finiteness (真定理).** In spectral form `H_E(t₀)=∫_0^∞dm ρ(m)e^{-m t₀}/(2m²)`. With `ρ(m)~m^{2Δ-3}` (`Δ=4`), the integrand is `~m^{2Δ-5}e^{-mt₀}`; the `e^{-mt₀}` cuts off `m≳1/t₀`, so
```
H_E(t₀) ~ t₀^{-(2Δ-4)}   (Δ=4: ~ t₀^{-4}),   finite for every t₀>0,  NO a⁻⁴.
```
The `t₀→0` limit re-creates the `a⁻⁴` divergence (`∫_0^{t₀}` is the contact region); **the whole point is the cut at `t₀>0`.**

**Atom lower bound and the finite-time squeeze (真定理).** Since all `w_n≥0`,
```
H_E(t₀) ≥ w₁e^{-m₁t₀}/m₁   ⇒   m₁ e^{m₁ t₀} ≥ w₁/H_E(t₀).                       (Euclidean, 真定理)
```
Proper-time twin (the prompt's formula, in the paper's `s=m²` column):
```
H_P(τ₀) ≥ f₁²e^{-m₁²τ₀}/m₁² ⇒   m₁² e^{m₁²τ₀} ≥ f₁²/H_P(τ₀).                    (proper-time, 真定理)
```
**Single-atom equality** (both columns): if `ν` is one atom, `H_E(t₀)=w₁e^{-m₁t₀}/m₁` and `m₁e^{m₁t₀}=w₁/H_E(t₀)` **exactly**; likewise `m₁²e^{m₁²τ₀}=f₁²/H_P(τ₀)`. ✓

**Acid test (charter §1) — fails on the free massless field (non-vacuous), as required.**
- **Free massless `φ` (`m₁=0`, acid):** the atom sits at `m₁=0`, so `H_E(t₀)=∑ w_n e^{-m_nt₀}/m_n` has a `1/m₁=∞` term ⇒ `H_E(t₀)=∞` ⇒ the criterion's premise `H_E≤X` **fails**. Correctly refuses to certify a zero gap. ✓
- **Free massive `φ` (`m₁=m>0`):** single atom ⇒ exact equality; criterion certifies the (real) free gap. ✓
- **Free massive composite `:φ²:` (gapped, `f₁²=0`):** `w₁=0`, so `(I) w₁≥r₀` **fails** ⇒ no false certification (and the theory *is* gapped — the criterion is sufficient-not-necessary, as the corpus's `H*` already is). ✓

So the point-split is a **valid, non-vacuous tool**: YM-positive (uses the positive atom), free-massless-failing. It passes the acid test where H1's literal form fails (H1 needs the `χ≤1/λ` that is false here).

**Note on `t₀`-dependence (the subtlety that decides §3).** The bound `m₁² e^{m₁ t₀} ≥ w₁/H_E(t₀)` holds for **every** `t₀>0`. But:
- at `t₀→0⁺`, `H_E~t₀^{-4}` (perturbative, large) and the bound is **weak** (small `m*`);
- at `t₀≫1/m₁`, `H_E(t₀)` is dominated by the atom itself (`~w₁e^{-m₁t₀}/m₁`) and the bound becomes an **identity** `m₁e^{m₁t₀}≈m₁e^{m₁t₀}` — tautological;
- the **only useful window** is `t₀ ~ 1/m₁`, which is **`t₀` set by the unknown answer `m₁`**.

This `t₀~1/m₁` requirement is the second face of the stuck step (§3c, §5).

---

## §3 · THE CRUX — escape or H1-in-finite-form?

### §3.1 (II) is *strictly weaker* than H1 — 真定理

Let (II) be `∃t₀>0 fixed: H_E(t₀)≤X` and H1 be `χ=2∫_0^∞C_E≤1/λ`.
- **H1 ⇒ (II)** (with `X=1/λ`): `H_E(t₀)=∫_{t₀}^∞C_E≤∫_0^∞C_E=χ/2≤1/2λ`. ✓
- **(II) ⇏ H1**: the short-time integral `∫_0^{t₀}C_E` — exactly the contact region that makes `χ` diverge — is unconstrained by `(II)`. Add a contact bump on `[0,t₀]`: `χ→∞`, `H_E(t₀)` fixed. ✓
Hence **(II) ⟸ H1, strictly**: the point-split criterion is a **strictly weaker premise** than H1. For a *sufficient* condition this is **good** (we want the weakest satisfiable premise): H1 is *false* here (`χ~a⁻⁴`), while (II) is *finite and possibly true*. **This is the genuine contribution of the fix.**

### §3.2 Same species as H1 — the wall is the direction, not the divergence

`H_E(t₀)=∫_{t₀}^∞C_E` is an integral of a **positive, completely-monotone, reflection-positive** correlator (`C_E(t)=∫e^{-mt}dμ≥0`, decreasing, convex). **Positivity gives `C_E≥0`, monotonicity, and lower bounds (`ν≥0`); it never gives an upper bound on `∫_{t₀}^∞C_E`.** The needed tail **upper** bound is the **same one-sided variance/moment upper bound** that is H1's named species:
- At **small** `t₀` (`a≪t₀≪1/m₁`): `H_E(t₀)≈χ_{IR}` — an upper bound on the **physical** susceptibility = **H1 on `χ_R`** (the contact-subtracted H1), verbatim.
- At **`t₀~1/m₁`**: `H_E(t₀)` is a spectral-weight upper bound on the light states = the same "control the light-mode weight" species.
- It is **not** derivable from positivity alone: any naive bound (`H_E≤C_E(t₀)/m₁`, or `≤Ae^{-μt₀}`) is **circular** (needs `m₁`/`μ` = the gap). The only unconditional version is at **strong coupling**, where **Lemma SC** (`Var_ν(φ)≤(C/λ)∫|∇φ|²`, proven `β<1/48`) supplies exactly H1 ⇒ (II) — i.e. (II) inherits H1's *proof* only in the same strong-coupling regime, and is open at the scaling window `β≈2.5` for the **same reason H1 is**. So **(II) is H1 in finite form**: same species, same regime wall, strictly weaker statement.

### §3.3 Is there a `t₀` window where (II) is *provably easier*? (short-time OPE vs long-time need) — NO

The tempting idea: compute `H_E(t₀)` from a **short-time OPE**. But `H_E(t₀)=∫_{t₀}^∞` needs the **long-time** (nonperturbative) tail; the OPE controls the **short-time** `∫_0^{t₀}` — the **complement**. So:
- OPE bounds the **wrong integral** (it gives the contact/subtraction, i.e. `χ_R`, not the bound on the IR tail).
- At small `t₀`, `H_E~t₀^{-4}` is *perturbatively computable* but **large** ⇒ `X` large ⇒ `m*` tiny ⇒ **weak bound**.
- At the useful `t₀~1/m₁`, `H_E` is **nonperturbative** and bounding it **is** the gap.
Hence there is **no window where (II) is provably easier**: the computable window gives a useless (weak) bound; the useful window needs the answer. **精确缺失步** (the same one-sidedness as H1, now finite).

### §3.4 Verdict

| question | answer |
|---|---|
| Does point-split cure the `a⁻⁴` degeneration? | **YES** (真定理; `H_E(t₀)` UV-finite for `t₀>0`) — genuine, correct. |
| Is (II) weaker than H1? | **YES, strictly** (真定理) — so the *premise* is better-posed (satisfiable vs unsatisfiable). |
| Is (II) a *different species* from H1? | **NO** — same one-sided upper-bound species; at small `t₀` literally H1-on-`χ_R`. |
| Does point-split touch (I), the atom `f₁²≥r₀`? | **NO** — unchanged = `H*(r₀)` = the gap. |
| Net | **ESCAPES the degeneration; SAME WALL on the content.** It is **H1 in finite form**, with the important upgrade that the finite form is *satisfiable*. |

---

## §4 · Gauge invariance and OS — no obstruction (真)

`O` is the gauge-invariant `0⁺⁺` (`A1g`, `C=+`) plaquette trace; `C_E(t)=∑_x⟨O(x,t)O(0)⟩_c` is a reflection-positive two-point function at `(t,n)`-separation, so (Osterwalder–Schrader 1973/75; Osterwalder–Seiler 1978) it satisfies `C_E(t)≥0` and has the KL representation with `ν≥0`. 

- **The "point split" here is a TIME cut** (`t₀>0` as the lower limit of `∫dt`), **not** a spatial displacement of the operator. A time cut is a positive linear functional of an RP correlator (`H_E(t₀)=∫_{t₀}^∞C_E≥0`), so it is **OS-compatible with no obstruction**; it does not change the operator, the gauge content, or reflection positivity.
- A genuine **spatial** point-split of `tr F²` (separated operators) would require a **Wilson line** between the split points to stay gauge-invariant; this is standard and OS-compatible too, but is **not** what the criterion uses (the susceptibility is already zero-momentum/`∑_x`, so no spatial separation arises).
- The only "obstruction" is the one already noted: `H_E(t₀)` is a *derived functional*, not itself an OS-correlator — harmless, since the inequality is a positive combination of RP quantities.

---

## §5 · Precise statements, labels, exact stuck step

**Definitions (fixed).** `C_E(t)=∑_x⟨O(x,t)O(0)⟩_c=∑_n w_n e^{-m_n t}`, `w_n=|⟨0|O|n⟩|²/(2m_n)`, `m₁=inf supp ν`, `w₁=|⟨0|O|G⟩|²/(2m₁)`. Measurable, normalization-invariant ratio: `R(t₀,a):=w₁(a)/H_E(t₀,a)` (both scale `c²` under `O→cO`).

**(I) Atom / residue floor.** `∃r₀>0`, uniform in `(a,L)`: `w₁≥r₀` (equivalently `f₁²≥2m₁r₀`, i.e. `H*(r₀)`).
**Status: 条件定理 (open).** This is the unprotected residue = the gap content. Unchanged by the point-split. *(Numerically seen at `β=2.2`: `w₁≈7.1` for bare `tr F²`, COMPUTE-F1SQ; not a continuum statement.)*

**(II) Finite-time tail bound.** `∃t₀>0` (fixed physical) and `X<∞`, uniform in `(a,L)`: `H_E(t₀)=∫_{t₀}^∞C_E(t)dt≤X`.
**Status: 条件定理 (open); 真定理 that it is UV-finite.** A one-sided upper bound on a positive RP correlator tail — **H1's species in finite form**. Proven only at strong coupling (via Lemma SC ⇒ H1 ⇒ (II)); open at `β≈2.5` for H1's reason.

**Conditional theorem (真, given (I),(II)).** `(I)∧(II)` at a common `t₀` ⇒
```
m₁ e^{m₁ t₀} ≥ w₁/H_E(t₀) ≥ r₀/X  ⇒  m₁ ≥ m* > 0,   where m* e^{m* t₀}=r₀/X.
```
(Proper-time twin: `m₁²e^{m₁²τ₀}≥f₁²/H_P≥r₀'/X'` ⇒ `m₁≥m*`.) Since `m↦m e^{mt₀}` is strictly increasing, `m*` is well-defined and `>0`. ✓ **So a gap follows from (I)∧(II)** — a clean **条件定理**.

**The exact missing step (精确缺失步).** Two halves, only one of which the fix helps:
```
       (I)  w₁≥r₀  : the unprotected atom  = H* = the gap.            [UNCHANGED by point-split]
       (II) H_E(t₀)≤X  at a FIXED physical t₀~1/m₁ :                  [H1-in-finite-form]
            — one-sided UPPER bound on a positive RP correlator tail;
            — NOT from positivity; circular without the gap;
            — the useful window t₀~1/m₁ depends on the unknown m₁.
```
**What the fix DID:** replaced the **unsatisfiable** `χ≤1/λ` (真 no-go on `tr F²`, §1) with the **finite, satisfiable** `H_E(t₀)≤X` — the premise is now well-posed. **What it did NOT do:** it left **(I)** untouched (= the gap), and **(II)** remains the same one-sided upper-bound species as H1 (§3.2). **It is H1 in finite form: better-posed, same species, atom wall unmoved.**

---

## §6 · The numeric kill-switch (precise protocol)

**Object (gauge-invariant).** `O(x,t)=` the `A1g,C=+`-projected spatial-plaquette trace, `q(x)=1−(1/6)∑_{μ<ν}½Re tr U_p` (or the standard `∑` over the three spatial planes `{12,13,23}` and both time layers), **APE-smeared spatial links** (as in COMPUTE-F1SQ: levels `s=0,2,4,8`, `α=0.5`) to kill the contact and expose the atom. Correlator `C_E(t)=∑_x⟨q(x,t)q(0,0)⟩_c`, connected, zero spatial momentum, periodic time `T`.

**τ range.** Fixed **physical** `t₀=τ/a` with `τ` in the glueball-correlation window `τ∈[τ_min,τ_max]`, `τ_min≳2a` (past contact), `τ_max≲1/m_G` (before noise); for `m_G≈1.5 GeV`, `τ~0.1–0.25 fm`. On the lattice: `t∈[2,4]` at `a≈0.08 fm`. **Scan `t₀` over a few values** to expose the `t₀`-dependence of §3.4.

**Measured per ensemble (fixed `β` ⇒ fixed `a`, single volume `L`):**
1. Fit `C_E(t)=w₁(e^{-m₁t}+e^{-m₁(T-t)})+(\text{excited})` over the plateau ⇒ `w₁(a)`, `m₁(a)`.
2. `H_E(t₀,a)=∫_{t₀}^{T/2}C_E(t)dt` (trapezoid from data, or from the fit).
3. Invariant ratio `R(t₀,a)=w₁(a)/H_E(t₀,a)`.
4. Solve `m*(t₀,a)e^{m*(t₀,a)t₀}=R(t₀,a)` for the derived lower bound `m*`.

**Across `a` (a β-ladder into the scaling window) and `L` (finite-size stability).**

**Pass / fail.**
- **VERIFY-UV (do on existing `L≤8`, β-scan 1.0–2.4, NOW):** confirm `χ_lat=V·Var(q̄)≈const=c` while `χ_phys=χ_lat/a⁴` grows `~a⁻⁴`, and `H_E(t₀)<∞` for `t₀>0`. **Pass = audit's `a⁻⁴` confirmed; contact-dominated χ confirmed.**
- **KILL-SWITCH (needs NEW `L≥16–24` scaling-window ensembles):**
  - **PASS (route survives):** `m*(t₀,a)≥m*>0` bounded below uniformly as `a→0`, for **some** fixed physical `t₀` in the window — i.e. `w₁(a)/H_E(t₀,a)` does **not** decay to 0 (and `w₁(a)` does not vanish).
  - **FAIL (route FALSIFIED):** `m*(t₀,a)→0` as `a→0` for **every** `t₀` in the window — equivalently `w₁(a)/H_E(t₀,a)→0` (the atom's share of the finite-time correlator vanishes) **or** `w₁(a)→0` (no surviving atom). Then **(I)∧(II) cannot yield a uniform gap, and the route is dead.**
- **Scope of the verdict.** FAIL falsifies the **route/criterion** (its certifying power); it does **not** disprove the gap (the criterion is sufficient-not-necessary — a real gap with `w₁→0` in the continuum is possible, precisely RESIDUE-CONDITIONAL's open "uniform `r₀` surviving `a→0`"). PASS is **numerical evidence**, not proof.

**Feasibility (honest).** Existing `L≤8`, `β≤2.4` data are **contact-like** for the `0⁺⁺` (`m₀≳π`, `ξ≲a`; H1-NUMERIC-RIGHT §4) and **cannot resolve `w₁`** in the window — COMPUTE-F1SQ needed `β=2.2, L=6, N=4000`, smearing, to see `w₁≈7.1`; the continuum statement needs `a→0`. So the **VERIFY-UV** leg runs on existing data; the **KILL-SWITCH** leg needs **new** ensembles (`L≥16–24`, `β≈2.35–2.5`, several `a`). Statistics: delete-block jackknife; **one job at a time** (charter §4 memory rule).

**Is it decisive?** **Well-posed and decisive about the ROUTE:** it cleanly separates "atom survives / tail bounded" (route viable) from "atom share →0" (route dead), via a normalization-invariant, free-field-clean ratio — and the FALSIFY direction (`m*→0`) is unambiguous. **Not decisive about Clay:** PASS is evidence, and FALSIFY kills the criterion, not the gap. So the numeric is the *right* experiment for the *claimed* improvement, and it can genuinely kill the route — which is exactly the honest use (charter §2: over-claiming forbidden, falsification valuable).

---

## §7 · Scorecard and sharpest line

| test | verdict |
|---|---|
| Is `χ~a⁻⁴` real? | **YES, 真定理** (dimension `Δ≥4` ⇒ `χ~a^{-(2Δ-4)}`); contact-dominated; measured `χ_lat≈c`. |
| Is literal H1 (`χ≤1/λ`) achievable for `tr F²`? | **NO — 真 no-go** (unsatisfiable). |
| Does point-split cure the degeneration? | **YES** (真定理: `H_E(t₀)` UV-finite for `t₀>0`). |
| Is (II) a *new* species, or H1-in-finite-form? | **H1 in finite form** — same one-sided upper-bound species; strictly weaker (真). |
| Does point-split touch the atom `f₁²≥r₀`? | **NO** — unchanged = the gap (`H*(r₀)`). |
| Acid test (free massless fails, YM-positive)? | **PASS** (criterion rejects the massless field; single-atom equality; `:φ²:` correctly uncertified). |
| OS / gauge invariance obstruction? | **None** (time cut on an RP correlator). |
| Numeric kill-switch well-posed / decisive? | **Well-posed & route-decisive**; needs `L≥16–24`; not gap-decisive; FALSIFY = route dead. |
| Clay | **0%** — unchanged; the atom `f₁²≥r₀` (with `a→0`) is still the wall. |

**Sharpest line.** The point split is a **correct fix to a real bug** — `χ~a⁻⁴` is genuine (a 真 no-go for literal H1 on the dim-4 gauge-invariant `tr F²`, and on *every* pure-YM gauge-invariant since `Δ≥4`), and `H_E(t₀)=∫_{t₀}^∞C_E` is genuinely UV-finite, so the sufficient condition becomes **satisfiable** where it was not. But `H_E(t₀)≤X` is **H1 with the divergent short-time contact excised** — the **same one-sided upper-bound species** (at small `t₀` literally H1-on-`χ_R`; provable only where H1 is, at `β<1/48` via Lemma SC). The hard half, the atom `f₁²≥r₀`, is **verbatim the gap and is untouched**. So the method **escapes the degeneration, not the wall**: it converts an *unsatisfiable* sufficient condition into a *satisfiable* one whose remaining content **is** `H*(r₀)`. Exact stuck step: **(I) `w₁≥r₀`** [unprotected atom, `a→0`-uniform] **∧ (II) `H_E(t₀)≤X` at a fixed physical `t₀~1/m₁`** [H1-species tail bound, circular without the gap].

---

*Refs.* **KL / OS:** Källén 1952; Lehmann 1954; Osterwalder–Schrader, CMP **31** (1973) 83 & **42** (1975) 281; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Seiler, LNP **159** (1982). **Strong-coupling gap:** Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737. **Internal — the degeneration:** `GAP-G2-edge.md` (dim-4 total moment diverges; only the anomaly log is scheme-independent), `H1-NUMERIC-RIGHT.md` §2–§3 (`Var(q̄)=χ_lat/V`, `χ_lat≈0.04` finite, `G″_dens=V/χ→∞`), `H1-IS-IT-THE-GAP.md` (squeeze `f₁²/χ≤m₁²≤M₀/χ`; H1 = inverse susceptibility; thin-atom witness), `REDUCTION-AUDIT.md` (H1 necessary-not-sufficient; `χ=∫s⁻¹dν` is a susceptibility). **Internal — the atom/species:** `RESIDUE-CONDITIONAL.md` (`H*(r₀)`; uniform `r₀` across `a→0` is the open content), `PROOF-ISOLATION.md` (L1 = I1∧I2; (P3) volume-uniform residue), `PROOF-F1SQ.md` (gapped `:φ²:` has `f₁²=0`; isolation `Λ>R̄`), `GAP-G2-edge.md` (residue–level-gap lemma; moments give only upper bounds), `G2-ASSEMBLE-100.md` §1–§2 (one species, curvature face H1 / edge face H*). **Internal — the proof at strong coupling:** `PROOF-COARSE-POINCARE.md` (Lemma SC ⇒ H1 ⇐ β<1/48). **Internal — measured residue:** `COMPUTE-F1SQ-SU2.md` (`w₁≈7.1` at β=2.2, L=6; contact-like at β=1.0). **Internal — tool registry / master:** `docs/YANG-MILLS-MASTER.md` §3,§5; `H1-DIRECT-SPLIT.md` (`G″=E[H|0]−Var(f|0)`, `G″·Var=1` iff Gaussian).
