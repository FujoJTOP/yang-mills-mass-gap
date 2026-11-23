# BS-VARIATIONAL-LOWER — does the J=0 two-gluon Bethe–Salpeter structure turn Temple's `E₂` into a *kinematic* two-body threshold?

> 2026-11-23 · 4D pure SU(N) YM · angle: variational LOWER bound on the `0⁺⁺` two-gluon BS eigenvalue via Temple/Weinstein, "next level" steelmanned as the kinematic two-gluon threshold `(2m_g)²`.
> **判词 (Verdict): NOT 5/5.** The only way Temple escapes `E₂`-circularity (TEMPLE-RIGOROUS-H §2) is a **kinematic** edge `ρ=4m_g²`. Unavailable. In the physical `0⁺⁺` channel the level above the lightest glueball is the **first radial excitation `0⁺⁺*`** — a dynamical BS bound state, not a continuum threshold — so `E₂=m²_{0⁺⁺*}` and `ρ` is again `I2`. The constituent "two-gluon threshold" `4m_g²` is either non-spectral or sits **above** the true `E₂` (invalid `ρ`). The two-body structure **relabels** the wall as `m_g>0`; it does not cross it.

## 1. Setting and Temple's requirement
Physical `0⁺⁺` is `1⁻⁻⊗1⁻⁻⊃0⁺⁺` (K0 §2). Spectrum `m₁²<m₂²<…`, target `m₁²=m²_{0⁺⁺}`. Temple: if `ρ≤m₂²` and `m₁²≤R<ρ`, then `m₁² ≥ R−Var/(ρ−R)`, `R=⟨ψ,Hψ⟩`, `Var=‖(H−R)ψ‖²`. The RHS **increases** in `ρ`, so sharpness `=` the best provable **lower bound on the second level** `m₂²`. The whole question: is `m₂²` kinematic?

## 2. The steelman worth taking seriously
For a genuine two-body problem with **asymptotic** constituents of mass `m_g`, a two-constituent continuum opens at `s_th=(2m_g)²`; if the ground state is the **only** bound state, `m₂²=s_th` is **kinematic** — exactly the Bazley–Kato situation that gives rigorous He/helium lower bounds. Temple then gives a real binding-energy bound, and the free-field test passes **with a nonzero result** (`m_g>0`⇒`ρ>0`; `m_g=0`⇒vacuous). This is the entire appeal of the two-body route, and it is genuinely how atomic bounds work.

## 3. Attack 1 — "next level = continuum" is FALSE for `0⁺⁺`
The channel has a **tower of bound states**: `0⁺⁺`, `0⁺⁺*` (radial, ≈2 GeV), `f₀(2200)`, … Each is itself a `1⁻¹⊗1⁻¹` composite; the BS equation `G₀⁻¹φ=Kφ` has radial solutions `n=0,1,…` with masses **below** the multigluon continuum. So `m₂²=m²_{0⁺⁺*}`, a **dynamical** eigenvalue. Verified on a two-body radial Schrödinger `V=−g e^{−r/a}` (Dirichlet box, `N=3000`): a one-bound-state well (`g=1.5`) has `m₂²`=continuum onset and kinematic `ρ=0` valid; a three-bound-state well (`g=6.0`, YM-like) has `E[1]=−0.572` far below the continuum `0`, and `ρ=0` **overclaims** — trial `ψ=0.6φ₀+0.7φ₁` gives `rhs=−1.74 > E₀=−2.16` (**false lower bound**), while `ρ=E[1]` gives exactly `E₀`. The excited bound state, not the threshold, is the next level.

## 4. Attack 2 — the "constituent" threshold is not physical
`m_g` is the gauge-fixed dressed-gluon mass (`D(p²)~1/(p²+m_g²)`), **not an observable**, and single gluons are **not asymptotic** (confinement). The physical continuum opens at the **two-glueball / string-breaking** scale `(2m₁)²` — kinematic only *given* the gap `m₁`, i.e. circular (`=I2`). In the confining two-body picture the flux tube **breaks** before any free two-gluon cut: the threshold is the **dynamical string-breaking scale**, not `4m_g²`.

## 5. Attack 3 — even granting `m_g`, `ρ=4m_g²` has the wrong sign
Temple needs `ρ ≤ m₂²`; both fail: (i) `m₂²=m²_{0⁺⁺*}<(2m₁)²`; if `m_{0⁺⁺*}<2m_g` then `4m_g²>m₂²` ⇒ **invalid `ρ`** (numerics: overclaim grows with admixture, `−1.74, −1.04, −0.63` vs `E₀=−2.16`); (ii) lattice/DSE `2m_g≈1.0–1.3 < m_{0⁺⁺}≈1.7` GeV ⇒ window `[m₁²,4m_g²)` is **empty**, the clean sub-threshold bound state Temple needs is absent.

## 6. Acid test — passed only by borrowing the gap
Form fails free field (`m_g=0`⇒vacuous), so not sign-blind in form; but all positivity comes from `m_g>0`, which is either (a) a gauge-fixed artifact with no rigorous gauge-invariant lower bound, or (b) the vector-channel gap `m_{1⁻⁻}>0` — a **different, unconquered** gap. The tool manufactures `m₁>0` only by importing a mass gap through `m_g`: the same circularity in constituent clothing.

## Five-point scorecard
1 Nonperturbative **PASS** · 2 sign-correct `m₁²≥λ` **PASS** (direction: the only lower-bound species) · 3 non-circular **FAIL** (`ρ≤m₂²=m²_{0⁺⁺*}` dynamical `=I2`, or `m_g>0` gauge/vector gap) · 4 fails free field **PASS (vacuous)** · 5 uniform `(V,a)` **FAIL** (`m_g` scheme-dependent; `Var_⊥` grows against the continuum).

## Exact stuck step
```
Temple:  m₁² ≥ R − Var/(ρ−R),  valid iff ρ ≤ m₂²,  m₁²≤R<ρ
need:    ρ ≤ m₂² = min( m²_{0⁺⁺*}, (2m₁)² )   [excited glueball — DYNAMICAL = I2]
OR       gauge-invariant m_g>0 with ρ=4m_g² ≤ m₂² (vector-channel gap — separate wall)
naive:   ρ = 4m_g²  ⇒ INVALID when m_{0⁺⁺*}<2m_g (numerically overclaims)
```
**Sharpest line.** The reformulation does real work — it names the target as a `J=1⊗J=1` composite and shows the *single*-gluon `E₂` was the wrong picture — but the `0⁺⁺` channel's second level is a **radial excitation, not a continuum threshold**. Bazley–Kato work because the next level above the ground state is the *kinematic* ionization continuum; here confinement replaces that cut with a **string-breaking** scale and fills the gap with excited glueballs. The threshold is **dynamical, ≈ the gap**; `E₂`-circularity is intact; the BS route relabels `H*` as "`m₂²` kinematic" only after assuming the atom it is meant to produce. **No 5/5; Clay stays 0%.**

*Refs.* Temple, *Proc. LMS* (1928); Weinstein, *Phys. Rev.* **40** (1932) 797; Kato, *J. Phys. Soc. Japan* **4** (1949) 334; Bazley–Kato, *J. Math. Mech.* **10** (1961) 291; Morningstar–Peardon, *PRD* **60** (1999) 034509; Cornwall, *PRD* **26** (1982) 1453. Internal: K0-TRANSVERSE-DECOUPLING §2; TEMPLE-WEINSTEIN §1–3; TEMPLE-RIGOROUS-H §2; GAP-G2-edge. Numerics: 2-body radial Schrödinger `N=3000` — 1-bound-state well: kinematic `ρ` valid; 3-bound-state well: kinematic `ρ=0` overclaims by `0.42–1.53`, `ρ=E[1]` exact. No fabricated bound.
