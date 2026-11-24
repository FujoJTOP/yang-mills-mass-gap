# INVENT-HOUSEKEEPING — a driven (non-equilibrium) entropy-production functional: outside Theorem BE, but killed by the **reversibility dichotomy**

> **VERDICT (no QED).** The housekeeping-heat / entropy-production species is genuinely **outside**
> Theorem BE (it is not a gradient-flow functional), but it **fails for an orthogonal reason**. The
> YM measure has a **real** action ⇒ **detailed balance** ⇒ the *intrinsic* entropy production is
> **identically zero** (Theorem HK, `[T]`); this is exactly the Theorem-BE fate continued, now in the
> non-equilibrium coordinate. The only way to get a nonzero functional is an **external non-conservative
> drive** `F`, which is (a) **controller-dependent** (violates `PROFILE` P4: the species must be
> controller-free) and (b) **maximally acid-dirty** — the driven functional grows as `m→0` (witness:
> `σ = 2γ²/m`, *largest* on the gapless field). The equilibrium (drive-free) "excess" EP reduces to the
> **Fisher metric = Theorem BE** (homogeneous). **Clay 0%.**

---

## §0 · Labels, iron laws, acid test (carried)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric · **[R]** inherited.
**IL-i** non-positivity (fail on free massless) · **IL-ii** nonabelian · **IL-iii** no small parameter. **Acid**: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.
**P4 requirement** (`PROFILE`): the source must be **intrinsic / controller-free**.

---

## §1 · The species, defined

Overdamped driven dynamics on a coarse gauge mode `X`:
```
dX = b(X) dt + √(2D) dW,     b = −D∇S_W + F,     F non-conservative (curl F ≠ 0).
```
Stationary density `π` (Fokker–Planck). **Total entropy-production rate** (Speck–Seifert / Hatano–Sasa):
```
σ_total = ⟨ (b − D∇log π)ᵀ D⁻¹ (b − D∇log π) ⟩_π = Tr[ M Σ Mᵀ D⁻¹ ],   M = A + D Σ⁻¹,
```
for linear drift `b = A X` and Gaussian `π = N(0,Σ)` (`AΣ + ΣAᵀ + 2D = 0`). **Detailed balance**
`⇔ b = D∇log π ⇔ M = 0 ⇔ σ_total = 0`.

---

## §2 · Theorem HK (reversibility rigidity) `[T]`

**The Yang–Mills equilibrium measure is reversible.** The Wilson action `S_W = −(β/N)Re tr U_□` is **real**,
so `μ ∝ e^{−S_W}dU` satisfies detailed balance w.r.t. the gradient drift `b = −D∇S_W`: with
`D∇log μ = −D∇S_W = b`, one has `M = b − D∇log μ ≡ 0`, hence
```
σ_total(intrinsic) ≡ 0   for ALL β, ALL coarse mass κ.
```
**The intrinsic entropy-production functional is zero at the very measure it must describe** — vacuous,
and acid-vacuous (zero on free massless *and* `:φ²:`). The non-equilibrium coordinate does **not** escape
Theorem BE; it re-derives it (`σ = 0 ⇔ reversible ⇔ gradient-flow`).

**Witness `[V-num]` (A).** Reversible OU (`γ=0`): `σ = 0.000e+00` at `m = 0.03, 0.3, 1, 3` (exactly).

---

## §3 · The only nonzero route — an external drive — is controller-dependent AND acid-dirty `[NG]/[O]`

To get `σ > 0` one must add a **non-conservative drive** `F = γR` (`R` = rotation, `curl F = −2γ ≠ 0`,
invented). Then `A = −m I + γR` and:

- **(O-i) Controller-necessity.** `σ` is a function of the **invented** `F`. Witness (C): `σ/γ² = 2.000`
  exactly for all `γ` tested — **`σ` is *dialed* by `γ`, not determined by YM**. This is precisely the
  `PROFILE` P4 violation the whole campaign forbids (an imported anchor, not an intrinsic relation).
- **(O-ii) Acid-dirty, in the worst direction.** Witness (B): `σ = 2γ²/m` — at fixed `γ=0.2`,
  `σ(m=0.03)=2.67`, `σ(m=0.3)=0.267`, `σ(m=1)=0.08`, `σ(m=3)=0.027`. **`σ` is *largest* on the light /
  gapless field** (`m→0`), *exactly opposite* to the acid requirement (the functional must **vanish** on
  the free massless field, not blow up there). The functional measures the **drive's** dissipation, not
  the field's gap.

So the driven species carries a nonzero number — but it is a number about `F`, sign-positive on the
free field, and diverging as the field becomes massless. **It cannot source the gap; it sources the
drive.**

---

## §4 · The equilibrium (drive-free) "excess" EP reduces to Fisher = Theorem BE `[T]/[NG]`

Without a drive, the only non-equilibrium object left is the **non-adiabatic (excess) entropy production**
between two equilibrium states, which for `π_κ, π_{κ′}` (1D Gaussian `N(0,1/m)`) is
```
σ_na = (m − m′)² / m       (witness D: 0, 0.04, 1.0 at m′=1, 1.2, 2)
```
i.e. the **Fisher information metric of the mass**. It is `0` at `m′=m`, and (relative-change)² —
**homogeneous**, with no `κ`-independent additive structure. This is **exactly Theorem BE's class**
(information geometry = LSI/Poincaré), so it adds nothing new: the drive-free non-equilibrium functional
*is* the equilibrium functional relabelled.

---

## §5 · The one loophole — intrinsic non-reversibility — is topological and `θ=0`-killed `[NG]`

Is there a **non-reversible dynamics intrinsic to the YM measure** (no external drive)? A real action
forces detailed balance, so the only intrinsic non-reversibility is a **complex measure**, i.e. a `θ`-term
`e^{−S_W + iθ∫tr F∧F}` with `θ ≠ 0`. This is:
- **topological** (a total-derivative / winding term) — `PROFILE` P4 requires **IR-support** sensitivity,
  which a topological term does not supply;
- **`θ = 0`** in the physical theory ⇒ real measure ⇒ reversible ⇒ `σ = 0` (back to Theorem HK);
- and the associated complex-Langevin dynamics has the known **wrong-convergence** failure for nonabelian
  gauge theory.

So the loophole does not yield a source either.

---

## §6 · Gap-response variant `[NG]`

One might hope the **response** `dσ/dκ` (how EP reacts to a coarse-mass perturbation) is a sign-definite
`κ`-independent source. It is not: in the equilibrium limit the response is the Fisher metric (homogeneous,
Theorem BE); in the driven case it inherits the drive's `κ`-scaling (`∝ γ²/m²`), i.e. it tracks the
**invented drive** and diverges on the light field — the same two failures (controller-dependence, acid).

---

## §7 · Honest bottom line

- The **non-equilibrium entropy-production species is genuinely outside Theorem BE** — it is not a
  gradient-flow functional. **But it does not help.**
- **Exact failing step (the "reversibility dichotomy"):**
  *intrinsic ⇒ `σ ≡ 0` (Theorem HK, Theorem BE continued); external drive ⇒ controller-dependent (`P4`) +
  acid-dirty (`σ = 2γ²/m`, maximal at `m→0`).*
- **What would be needed:** a `κ`-independent source that is **intrinsic** yet **nonzero on a reversible
  measure** — impossible unless the measure itself breaks detailed balance intrinsically, which only a
  topological `θ`-term does, and that is not IR-support-sensitive and vanishes at the physical `θ=0`.
- **The wall persists**, now in the non-equilibrium coordinate. **Clay 0%; no QED.**

*Labels:* Theorem HK `[T]`; driven route `[NG]`+`[O]`; equilibrium-excess `[T]` (Fisher = Theorem BE);
θ-loophole `[NG]`; witness `[V-num]`.
*Refs:* `BASE-FIELDMEASURE.md` (Theorem BE), `ROOF-SYNTHESIS.md`, `PROFILE.md` (P4), `IR-RESPONSE-TEST.md`;
witness `w7/invent_housekeeping.py`.
