# STEKLOV — is the coarse stiffness a Dirichlet-to-Neumann / Steklov eigenvalue, and does its bottom give `G″(k→0)≥λ` uniformly?

> 2026-11-23 · NEW SPECIES (response / caloric frame). **判词: the frame names the target exactly and passes the free-field acid test, then lands on the wall in new coordinates.** The coarse stiffness *is* a Steklov eigenvalue: the Gaussian–Schur curvature `S=E_μ[h_n−⟨K⁻¹c,c⟩|0]` of the exact split `G″=S−D` (A2, INVENT-VARIANCE-INEQ) is literally the Dirichlet-to-Neumann eigenvalue of the block for the fiber Hessian `K` (`S=1/Var_μ(a)` on the Gaussian), so `G″ = σ_Steklov − D`. Two independent deaths: **(a)** the DtN operator is **Savvidy-indefinite** (`K` indefinite ⇒ no positive `λ` bottom); **(b)** the only positive Steklov gap `σ_1` controls the **stiff** (boundary) modes, while `G″(k→0)` lives at the **constant** mode `σ_0=0`. Weighting DtN by the fiber measure reproduces `S−D=G″`=H1. **#1 PASS; #2 FAIL (wrong mode + indefinite); #3 FAIL; #4 PASS (clean); #5 FAIL.**

**§1 The object is real — the coarse stiffness IS a Steklov eigenvalue.** Split `T=coarse⊕fiber`; `K=H_ff` (fiber Hessian at fixed boundary), `c=Π₀ Hess V v`, coarse bare `h_n`. The block DtN/Schur operator is `Λ=H_cc−c*K⁻¹c` on coarse boundary data; its Rayleigh eigenvalue `σ(v)=min_int[form]/‖v‖²_∂` is the Schur curvature — exactly `S(v)` in the exact identity `G″=E[h_n]−Var(f)=S−D`. The measure average `E_μ[·|0]` makes it the *weighted* Steklov eigenvalue the proposal asks for. **Measured:** β=2.5 SU3, `S` carries `G″=+0.51`, signed defect `D/S≈2%`.

**§2 Two independent deaths, in DtN coordinates.** *(a) Indefinite operator.* `Λ` is a bona fide **positive** DtN iff `K⪰0`; at Savvidy `K` is indefinite (W6-BE: 90–97% of negative BE directions in the fiber), so `Λ` is a **Krein-indefinite** DtN — not bottomed by any positive `λ`, and `D>0` is exactly the fiber-Savvidy injection. *(b) Wrong mode.* `σ=S=1/Var_μ(a)`. Block Steklov spectrum is `σ_0=0` (massless/Lebesgue: constant boundary data is harmonic, `∂_ν=0`, `Var(a)=∞`) `≤σ_1≤…`, with `σ_1(bΩ)=σ_1(Ω)/b→0`. The coarse target `G″(k→0)` *is* the constant-mode lift `σ_0=1/Var_ν(a_{k→0})`; a Steklov bound `σ≥σ_1>0` controls the modes **above** the constant — the stiff, block-positive / Lemma-SC / β<1/48 / Jacobian-half (`S^Φ≥N/12`) gap — **not** the constant mode. Lifting `σ_0` by the *signed* YM measure (`D>0`) is mass generation = H1.

**§3 Monotonicity — genuine asset, circular closure.** DtN is monotone (`A_1≥A_2⇒Λ_1≥Λ_2`): a reference `A_ref≤A_YM` with uniform positive constant-mode DtN would lift by order. But `A_YM=−D²+2adF_B` has `λ_min=−g|B|` (transverse-magnetic tachyon), so an `A_ref` that is *both* `≤A_YM` *and* `≥λ>0` on all modes needs `A_ref−A_YM≥λ+g|B|>0` on that mode — contradiction. `A_ref≤A_YM` thus forces a **tachyonic** reference (`σ_0≤−g|B|+…`, no uniform lift); a positive reference mass equals `g|B|`, the tachyon scale — circular. This is PROOF-SAVVIDY-WGEO-DOMINATION in DtN language.

**Five-point scorecard.**
1. **Nonperturbative — PASS.** Exact DtN/Schur + A2 identity; no loop expansion.
2. **Sign-controlling LOWER bound — FAIL.** Indefinite `K` ⇒ no `λ`; `σ_1` controls the wrong (stiff) mode.
3. **Non-circular — FAIL.** Measure weight = `S−D=G″`; monotone reference = tachyon.
4. **Non-SUSY + fails free field — PASS (cleanest).** Massless: `K=−Δ⪰0`, `D=0`, `σ_0=0`, `σ_1∝1/b→0`, `G″=0−0=0` — spectrum accumulates at 0, no spurious gap; β=0 YM `≥N/12`; the discriminator is precisely the signed measure.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Only `σ_1∝1/b` (stiff, β<1/48) is provable; it shrinks in the IR and is not the constant-mode lift.

**Exact missing step.**
```
∃λ: σ_0^{YM}(k→0) = S − D ≥ λ  uniform in (L,a),
i.e. signed fiber defect D bounded so D ≤ S − λ at the constant mode (D/S≈2%).
```

**Sharpest line.** Steklov is a legitimate new species: it *names* the coarse stiffness as a Dirichlet-to-Neumann eigenvalue (`G″=σ_Steklov−D`) and *passes* the free-field test (massless Steklov accumulates at 0, `σ_0=1/∞=0`). But the DtN operator is Savvidy-indefinite — no positive bottom — and its only positive gap `σ_1` controls the stiff boundary modes, not the constant/coarse mode `σ_0=0` where `G″(k→0)` lives; lifting `σ_0` by the signed measure is `S−D=G″`=H1. New frame, same gap; the genuinely new content is the dictionary `G″=σ_Steklov−D` and the **mode-mismatch no-go**: the provable Steklov gap is the *stiff* gap, the target is the *constant*-mode lift.

*Refs.* Girouard–Polterovich, J. Spectral Theory **7** (2017) 321 (`σ_k(bΩ)=σ_k(Ω)/b`; accumulation at 0); Weinstock, J. Rational Mech. Anal. **3** (1954) 745 (isoperimetric Steklov = *upper* bound — wrong direction); Steklov 1902; DtN monotonicity (standard). Internal: INVENT-VARIANCE-INEQ (A2 `G″=S−D`, `S`=Gaussian–Schur/DtN, `S=1/Var_μ(a)`, `D/S≈2%`); PROOF-COARSE-POINCARE (`G″` exact; `H_ff` indefinite); OTHER-H1-POINCARE; B3-ISOPERIMETRIC; PROOF-DIRICHLET-COMPARISON (fiber-form kernel = coarse functions); PROOF-SAVVIDY-EXACT-SIGN (`λ_min=−g|B|`); W6-BE-at-2.5 (fiber 90–97% negatives); ANDERSON-LOCALIZATION (sibling frame); C8.
