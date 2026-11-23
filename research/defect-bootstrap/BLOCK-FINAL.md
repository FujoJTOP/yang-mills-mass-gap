# Block — adjudication: defect-channel ₂F₁ is correct; the "impossibility" was the wrong operator

> Two independent derivations disagreed (`BLOCK-EXPLICIT.md` vs `BLOCK-CROSSCHECK.md`). Settled by
> direct numerics. **The defect-channel block is a single ₂F₁** (E's result); A's "no block meets
> both ends" is true only for the **harmonic free-KG** operator — the wrong one.

## The disagreement

- **A (`BLOCK-EXPLICIT.md`):** `P_{Δ_j−1}` and `Q_{Δ_j−1}` are the two solutions of the *harmonic* H²
  Laplacian `f''+coth(d)f'=Δ_j(Δ_j−1)f`; no `αP+βQ` holds both OPE ends (proven); two channels
  related by crossing. Correctly flags `U-ODE.md`'s `f = P`.
- **E (`BLOCK-CROSSCHECK.md`):** the defect-channel block is a **single** ₂F₁
  `f(u) = (1−u)^{Δ_j}·₂F₁(Δ_j−Δ_φ, Δ_j−Δ_φ; 2Δ_j; 1−u)/N`, `N=Γ(2Δ_j)Γ(2Δ_φ)/Γ(Δ_j+Δ_φ)²`,
  governed by the **Δ_φ-deformed defect-channel Casimir**
  `u(u−1)²f''+(1−2Δ_φ)(u−1)²f'+Δ_φ²(u−1)f = Δ_j(Δ_j−1)f` — both OPE limits built in.

## Adjudication (independent numerics, Δ_φ=2)

| `Δ_j` | `f(0)` | `f/(1−u)^{Δ_j}` @ u→1 | ODE residual |
|---|---|---|---|
| 1.0 | 0.99999999 | → 1/N = 0.667 | ≤5e−8 |
| 1.5 | 0.99999999 | → 1/N = 0.920 | ≤5e−7 |
| 2.0 | 0.99999999 | → 1/N = 1.000 | ≤2e−7 |
| 3.0 | 0.99999999 | → 1/N = 0.800 | ≤1e−7 |

- `f(0)→1` (bulk OPE, normalized) ✓; `f/(1−u)^{Δ_j}` → finite constant (defect OPE, pure decay) ✓;
  ODE residual ~1e−7 ✓. **Hand-check** `Δ_φ=Δ_j=2`: `f=(1−u)²` solves the ODE exactly.
- ODE exponents: `{Δ_j, 1−Δ_j}` at `u=1` (defect primary + shadow), `{0, 2Δ_φ}` at `u=0`
  (bulk identity + double-trace `[φφ]`) — the physics is built into the operator.
- **A used the harmonic free-KG Laplacian** — the very operator `RADIAL-SOLVE.md` used and the
  `RESOLUTION` probe flagged as ≠ the conformal Casimir. Its "impossibility" is true for that
  equation and irrelevant to the physical defect channel. E's single ₂F₁ exists ⇒ refuted.

## What survives from each

- **From E:** the defect-channel closed form + normalization `N` (verified 3 ways: E's scipy,
  my independent run, a hand check). **This is the defect-channel block.**
- **From A:** the two-channel/crossing structure is right — bulk-OPE = bulk sum, defect-OPE = defect
  sum, related by crossing. But A's `Q` is **not** the defect-channel block (it is a harmonic-equation
  solution); the defect-channel block is E's ₂F₁.

## Channel bookkeeping (correct)

| channel | block | status |
|---|---|---|
| **defect** (defect ops `Ô`) | E's ₂F₁ `f(u)` above | **SOLVED** (closed form + norm, verified) |
| **bulk** (bulk ops `O`) | Dolan–Osborn bulk block restricted to defect kinematics | standard import (`SDPB-INPUT.md`); defect-kinematics reduction unverified |
| **crossing** | bulk-OPE sum = defect-OPE sum | the SDP constraint |

## Honest gaps

1. E **imports** the 1D defect-block ODE from the defect-OPE→1D-CFT identification — not re-derived
   from the 4D group Casimir (E flags this; validated by the exponents + numerics).
2. The **bulk-channel** block is still the Dolan–Osborn import; its defect-kinematics reduction is
   unverified.
3. The dynamical coefficients `a_j²` / `C^{φφ}_O` are the bootstrap unknowns (not kinematics).
