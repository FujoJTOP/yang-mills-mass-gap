# KAPPA-OVERDETERMINED — the single-κ rigid string is a ONE-datum fit, not a predictive spectrum: it fails the excited 0⁺⁺ (25% under) or the 2⁺⁺ (21% over, predicts a degeneracy not seen) — κ is under-determined, hence still = H1

> **VERDICT (no QED).** κ is fixed by ONE number (`m_{0⁺⁺}/√σ≈3.93` → κ≈0.196). Every other lattice observable either is in tension with the one-parameter rigid string or is κ-blind at leading order. **The model is not over-determined — it is under-determined**: its two natural "next 0⁺⁺" predictions bracket the lattice value but neither matches, and the closed-string tower predicts a 2⁺⁺/0⁺⁺\* degeneracy the lattice does not show. So κ cannot be pinned from the spectrum; it remains a fitted effective parameter, i.e. `H1` in string coordinates.

---

## §0 · Labels, and the honest status of the numbers

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric · **[V-model]** model-level. Lattice values are **representativeliterature** (Morningstar–Peardon 1999, hep-lat/9901004; SU(3), √σ≈440 MeV): 0⁺⁺≈1730, 2⁺⁺≈2400, 0⁺⁺\*≈2670 MeV. The excited-scalar 0⁺⁺\* carries a large error (±~180 MeV → ratio ±~0.15); the qualitative tensions below exceed it.

---

## §1 · The model's predictions (from CONSTRUCT-MASTERFIELD / FLUXTUBE-BREATHING / LUSCHER / WORLDSHEET-CFT)

Rigid closed flux loop `E(R)=2πσR+πκ/R`, interior minimum `R*=√(κ/2σ)`:
- **lightest 0⁺⁺** (loop ground state): `m_{0⁺⁺}/√σ = 2π√(2κ)`;
- **radial / breathing 0⁺⁺**: `m_breath/√σ = 2/√κ`;
- **closed-string oscillator tower** (`L=2πR*`, `c=D−2=2`): `ΔE=(N+N̄−1/6)/R*`, first level `(5/6)/R* = (5/6)√(2/κ)·√σ`, giving a **second 0⁺⁺ and a 2⁺⁺ degenerate** at that level.

These are three predictions from **one** κ. Over-determination = test whether they are simultaneously consistent with lattice.

---

## §2 · Over-determination table (witness `w7/kappa_overdetermined.py`)

| observable | lattice (√σ units) | model prediction at κ=0.196 | status |
|---|---|---|---|
| 0⁺⁺ lightest | 3.93 | 3.93 | **fit** (defines κ) |
| 0⁺⁺\* (excited scalar) | 6.07 (ratio 1.55) | radial: 4.52 (1.15); oscillator: 6.60 (1.68) | **tension**: radial −25.5%, oscillator +8.6% — neither matches |
| 2⁺⁺ | 5.45 (ratio 1.39) | oscillator tower puts 2⁺⁺ = 0⁺⁺\* at 6.60 | **tension**: +21%, and predicts a 2⁺⁺/0⁺⁺\* degeneracy the lattice splits (5.45 vs 6.07) |
| torelon (closed string) | `M=σL−π/(3L)` | same at leading order | **κ-blind** (κ enters only at O(1/L³)) |
| flux-tube width `w²(R)` | log slope `1/π` universal | log slope κ-independent; κ shifts only intercept | **κ-blind** (degenerate with UV core `r₀`; CONSTRUCT-WIDTH-KAPPA-LATTICE) |

**Only one row fixes κ.** Three rows are in tension or blind. There is no second, independent, κ-sensitive, model-consistent observable.

---

## §3 · Verdict

- **Not over-determined.** κ is fixed by `m_{0⁺⁺}` alone. The "spectrum" cannot constrain it further: the two candidate excited 0⁺⁺ (radial 1.15, oscillator 1.68) **bracket** the lattice 1.55 but neither matches, so the model gives **no unique excited-state prediction** to pin κ.
- **Soft falsification as a spectrum model.** If the first excited 0⁺⁺ is the breathing mode → **25% under**; if it is the closed-string oscillator level → the model then predicts **2⁺⁺ = 0⁺⁺\*** (degenerate) at 6.60√σ while the lattice has 2⁺⁺ = 5.45 (a 0.62√σ split) → **21% over**. Either way one excited state is badly wrong. The one-parameter rigid string describes the lightest 0⁺⁺ (by construction) but **does not predict the glueball spectrum**.
- **κ therefore remains under-determined** — it is a fitted effective parameter for the lightest mode, not a quantity the data pin. It must come from the theory (the microscopic flux-tube core), i.e. **κ = H1** in string coordinates.

---

## §4 · What the "ferocious attack" establishes

The construction-round hope was: "one model family, one missing input κ, over-determined by the glueball spectrum." The over-determination **fails in the honest direction**: the data **under-determine** κ (one datum), and the model's spectrum is in tension. This is a real, useful negative:
- it **closes** the route "pin κ from many glueball states";
- it **confirms**, one more time and now empirically, that κ is not a phenomenological handle but the microscopic nonperturbative core = `H1`;
- it leaves the rigid string as a **one-parameter phenomenological description of the lightest 0⁺⁺**, not a predictive spectrum.

**Clay 0%; no QED.**

*Refs:* `CONSTRUCT-MASTERFIELD.md`, `CONSTRUCT-FLUXTUBE-BREATHING.md`, `CONSTRUCT-LUSCHER.md`,
`CONSTRUCT-WORLDSHEET-CFT.md`, `CONSTRUCT-WIDTH-KAPPA-LATTICE.md`, `CONSTRUCT-KAPPA-DERIVE.md`;
lattice: Morningstar–Peardon, hep-lat/9901004. Witness: `w7/kappa_overdetermined.py`.
