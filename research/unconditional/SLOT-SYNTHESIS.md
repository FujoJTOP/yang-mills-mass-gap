# SLOT-SYNTHESIS — the complete face of the `k→0` sign-definite self-consistency (8 routes, one profile)

> **VERDICT (no QED).** Eight independent attacks — Birman–Schwinger, Dyson–Schwinger, Riccati/ERG,
> `G″=S−D`, Ward/Stein, Perron–Frobenius transfer, numerical bootstrap, quadratic-minorant fixpoint —
> were each told to *construct* a `k→0` sign-definite relation `R` that (a) excludes the gapless fixed
> point and (b) attains uniformly at `β≈2.5`. **None constructed it.** All eight converge, to the exact
> failing step, on **one and the same object**: a **sign-definite, `Q`-independent *source at the zero*** —
> equivalently spontaneous scale breaking = dimensional transmutation = `H1 = I2 = H*`. The "equation"
> is infinite in every coordinate (vertex tower / Makeenko–Migdal / A3–Fredholm / field-mode data), but
> its **irreducible core collapses to a single sign** — and that single sign is the Clay gap. **Clay 0%.**

---

## §0 · Labels, iron laws, acid test (carried unchanged)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-here]** new this round · **[R]** inherited.
**IL-i** non-positivity (must fail on free massless: `f₁²=0`, `G″=0`). **IL-ii** nonabelian (U(1) massless). **IL-iii** no small parameter (must hold at `α_s~O(1)`, `β≈2.5`).
**Acid**: `Φ(free massless)≤0` *and discriminatingly* `≤0` on the gapped-atomless `:φ²:` composite (`f₁²=0`).

---

## §1 · The scorecard (8 routes × 4 tests)

| # | Route (file) | the relation `R` | `F(0)=0`? (gapless solves) | sign-definite source at `κ=0`? | uniform at `β≈2.5`? | closest verdict |
|---|---|---|---|---|---|---|
| 1 | **Birman–Schwinger** `SLOT-BS-FIXPOINT` | `κ=F(κ)`, root of `det(1−K_BS)`, edge self-ref `θ=4κ` | **YES** (`F(0)=0`) | needs `F(0)<0 ⇔ K_BS(θ)>1` = `0⁺⁺` tachyon — **absent** | NO (`=I2`, non-contractive) | **[NG]+[O]** |
| 2 | **Dyson–Schwinger** `SLOT-DSE-FIXPOINT` | `m²=Π_T(0)`, 2PI tower | **YES**; seagull **forces** `m²=0` | scheme-indep. piece *certifies* zero; escape = vertex-pole residue `C(r²)` (scheme-dep.) | NO (lattice↔continuum = H1) | **[NG]+[O]** (Thm D: acid–scheme dichotomy) |
| 3 | **Riccati / ERG** `SLOT-RICCATI-FLOW` | `∂_kκ=2σ₀κ²+(3σ₂−2β₁)κ+C` | **YES** (`κ≡0`, free traj.) | `κ²` marginal; source `C` sign-indefinite at crossover | NO (`κ*` `k`-dependent) | **[NG]+[O]** |
| 4 | **`G″=S−D`** `SLOT-SD-FIXPOINT` | `κ=S_closed(κ)`, `D=0` (Gaussian fiber) | **YES** (`κ=0` contraction attractor) | `W_geo=0` **on the acid field**; needs `D<S` | NO (`H1`) | **[NG]+[O]** |
| 5 | **Ward / Stein** `SLOT-WARD-STEIN` | `κ·Var(a)=1−η` (`R` ≡ "p is a normalized pdf") | **YES** (zero-content; `G″` spans ±) | none extractable (READ not FORCE) | NO (`E[H\|0]≥λ+Var`) | **[O]** (IL-i read-not-force) |
| 6 | **Perron–Frobenius** `SLOT-PF-FIXPOINT` | `ρ=λ₁/λ₀`, no relation `ρ=F(ρ)` exists | **YES** uniformly (`ρ_L→1`) | single cone **sign-blind to `λ₁`**; needs 2nd cone or uniform `Φ≥δ` | NO (coarse Poincaré) | **[NG]+[O]** |
| 7 | **Bootstrap** `SLOT-BOOTSTRAP` | `Δ_min(K)` over Stieltjes moments | **YES** (`Δ_min(K)=0` ∀K) | moments carry **upper** bounds only; need edge datum `f₁²≥r₀`+`Λ` | NO (needs field-measure data) | **[NG]** finite-moment class |
| 8 | **Quadratic minorant** `SLOT-QUADRATIC-FIXPOINT` | `Q=max{q⪰0:q≤E[Hess]}` | **YES** (`Q=0 ⇔ ∂_vV≡0`) | pure YM source-free at `Q=0`; anomaly **homogeneous ⇒ cancels** | NO (`Q*=b−g²Var<0`) | **[O]** |

**0 / 8 constructed. 8 / 8 pinned the failing step.**

---

## §2 · The face: six convergences

The eight reports do not merely fail alike; they fail *identically*. Strip the coordinates and the same
six features appear on every one:

1. **The zero is universally a fixed point — and usually the *attractor*.**
   `F(0)=0` in every route: BS (`F(0)=0`), DSE (`m²=0` forced), Riccati (`κ≡0` free trajectory),
   S−D (`κ=0` contraction attractor), Ward/Stein (`G″` unconstrained, zero admissible), PF (`ρ_L→1`),
   bootstrap (`Δ_min=0` ∀K), quadratic (`Q=0` admissible). **No coordinate removes the gapless solution.**

2. **Every scheme-independent / positivity-clean / exact identity is READ, not FORCE.**
   Ward/Stein collapses to "`p` is a normalized pdf"; DSE's one scheme-independent piece (seagull)
   *certifies* `m=0`; PF's single cone is sign-blind to `λ₁`; bootstrap's exact moments bound the floor
   from *above* only; the variance identity `G″=E[H|0]−Var(f|0)` *subtracts*. **These are the `𝔉_*` of
   `PROFILE` P2: lsc, `≤0`, carry zero content.** They hold *at* the zero.

3. **The discriminating piece is always on the scheme/measure-dependent side = `H1`.**
   BS → vertex-pole residue / `K_BS(θ)>1`; DSE → `C(r²)`, unprotected by STI; Riccati → `F′(0)>1`;
   S−D → `D<S`; Ward/Stein → `E[H|0]≥λ+Var`; PF → uniform Cheeger `Φ≥δ`; bootstrap → edge datum
   `f₁²≥r₀`+isolation `Λ`; quadratic → `Q`-independent additive source. **All eight are the same
   inequality in eight coordinates.** (This is `PROFILE` §3's "representational freedom is illusory",
   now measured, not asserted.)

4. **What is missing is a *source at the zero*, not another identity.**
   The needed object is a term that is **nonzero exactly where the field is free/gapless** and has a
   **definite sign** there, so the relation is *inhomogeneous* (`κ=F(κ)` becomes `κ=F(κ)+c`, `c≠0`).
   Every route states it: quadratic → "pure YM source-free at `Q=0`"; S−D → "`W_geo` is exactly zero
   on the field it must repel"; Riccati → "source `C` sign-indefinite at `κ=0`"; Ward/Stein → "`R` is
   READ not FORCE"; DSE → "the acid-clean piece *is* the gapless solution". **A homogeneous relation
   always keeps its zero; only an inhomogeneous sign-definite source can exclude it.**

5. **That source, if it exists, *is* spontaneous scale breaking = dimensional transmutation.**
   Pure YM has no bare mass (quadratic); the trace anomaly — the sole scale-breaker — enters
   **homogeneously and cancels** (quadratic, Ward/Stein); the only sign-definite nonabelian object that
   is acid-clean (`Haar Φ″⪰N/12`) floors `W_geo`, *not* `G″`, and its `β=2.5` promotion is the same
   inequality (S−D, Ward/Stein). So the missing source is **the mass scale the nonabelian measure
   manufactures at `α_s~O(1)`** — which is, verbatim, the statement that the theory gaps.

6. **The required escape mechanism is a `k→0` `0⁺⁺` instability — and at `β=2.5` it is absent.**
   Exclusion needs a tachyon/negative mode in the `0⁺⁺` channel (BS: `K_BS(θ)>1`; Riccati/S−D:
   `F′(0)>1` = Savvidy tachyon). But the Savvidy negative mode is **`1⁻⁻` transverse magnetic**, projected
   out of `0⁺⁺`; the coarse block is measured **positive** (`G″≈+0.25…+0.51`), the rigorous one-sided
   bound still touches zero, and 90–97% of the negative weight is high-momentum (W6-BE). **Necessary-but-
   not-sufficient: correct channel handling removes the wrong-sector sign without creating a `k→0` sign.**

---

## §3 · "How big is the equation?" — resolved by the eight

Each route reports how its `R` fails to close, i.e. how big the exact relation is:

- **Riccati** does not close: `(R)` couples to `φ₃, λ, ∂²σ, ∂²β` = the **vertex tower**.
- **DSE** does not close: full gluon+ghost loops need the **infinite Makeenko–Migdal hierarchy**; the
  2PI kernel trades a vertex truncation for a kernel truncation — same infinity.
- **Birman–Schwinger** does not close: `det(1−K_BS)` = the **Fredholm `D_O(z)` = A3**, which is *read*,
  making the route **circular**.
- **Bootstrap** does not close: finite moments are provably floor-blind (`Δ_min(K)=0` ∀K); a closing
  input must be **edge-sensitive field-measure data**, which is not a moment.

⇒ **The exact `R` is coextensive with the whole QFT (infinite-dimensional in every coordinate).** But the
eight collapse its *irreducible core* to **one scalar sign condition at `κ=0`**: *"is there a
sign-definite, `Q`-independent source at the zero?"* **The equation is infinite; its core is one sign;
and that one sign is Clay.**

---

## §4 · What would actually fill the slot (the sharpest honest target)

To fill it, construct a **field-measure functional** `𝒮(μ)` (not a spectral/scheme statistic) such that:

- **(S1)** `𝒮(μ_free massless)=0` *and* `𝒮(:φ²:)=0` (acid-clean on both adversary classes), yet
- **(S2)** `𝒮(μ)>0` for the Haar-carrying interacting measure at `β≈2.5`, with
- **(S3)** `𝒮` **inhomogeneous** in the coarse mass: `κ=F(κ)+𝒮` with `𝒮` independent of `κ` (a genuine
  source, not a coefficient), so `κ=0` is **not** a solution, and
- **(S4)** the resulting `κ*` bounded below uniformly in `(L,a)` across `[1/48, 2.5]` (= `I2`, the
  `a→0` zero-finding uniformity = the transport certificate).

`(S1)` excludes all positivity/scheme tools (they are homogeneous or vanish on the adversary). `(S3)` is
the missing structural ingredient — *every* route's failure is "no `κ`-independent source". `(S4)` is the
transport wall. **No candidate `𝒮` was constructed; this is the exact missing step.**

---

## §5 · Honest bottom line

- **The slot is not filled.** Eight routes, eight honest `[NG]`/`[O]`. No QED.
- **The slot is now fully *drawn*.** The "五官" is complete: *(i)* zero is always a fixed point;
  *(ii)* all exact/scheme-clean identities are READ-not-FORCE and hold at the zero; *(iii)* the
  discriminating piece is always scheme/measure-dependent `=H1`; *(iv)* the missing object is a
  **sign-definite source at the zero**; *(v)* that source is spontaneous scale breaking; *(vi)* the
  would-be escape (a `k→0` `0⁺⁺` tachyon) is absent at `β=2.5`.
- **The equation's size:** infinite in every coordinate, **core = one sign**. This is the sharp resolution
  of "how big is it": as a *system* it is the whole theory; as a *decision* it is a single inhomogeneity.
- **`PROFILE` is confirmed and sharpened:** P1 (field measure, not spectral data) — bootstrap's
  `Δ_min=0`∀K is the proof; P4 (excludes the zero via a `k→0` sign) — all eight fail exactly there.
- **Clay = 0%.** The wall did not move; it was measured from eight directions and found to be one point.

*Labels:* every route `[NG]` (its named class) + `[O]` (exact failing step). Structural claims (`F(0)=0`,
read-not-force, acid-scheme dichotomy, single-cone sign-blindness, `Δ_min(K)=0`) are `[T]/[V-num]` as
detailed in each file. **No claim of gap-forcing identity; no QED.**

*Refs:* `SLOT-BS-FIXPOINT.md`, `SLOT-DSE-FIXPOINT.md`, `SLOT-RICCATI-FLOW.md`, `SLOT-SD-FIXPOINT.md`,
`SLOT-WARD-STEIN.md`, `SLOT-PF-FIXPOINT.md`, `SLOT-BOOTSTRAP.md`, `SLOT-QUADRATIC-FIXPOINT.md`;
`PROFILE.md` (+ facets), `KEY-DIRECT-DERIVATION.md`, `META-GAPFORCING-IDENTITY.md`, `A3`/`PROOF-F1SQ.md`.
Witnesses: `w7/slot_*fixpoint.py`, `w7/dse_fixpoint_witness.py`, `w7/riccati_flow_witness.py`,
`w7/ward_stein_witness.py`, `w7/slot_bootstrap.py`.
