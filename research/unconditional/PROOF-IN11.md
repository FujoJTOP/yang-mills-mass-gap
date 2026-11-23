# PROOF-IN11 — strong-coupling mass ordering: a membership theorem, not the ordering

> **Verdict.** Strong coupling **proves the weak form (b′) — but not IN-11 proper.** For `0<β<β₁`, the lightest bulk state provably lies in `{0⁺⁺,2⁺⁺}` and every other `J^{PC}` is strictly heavier with a gap that diverges as `β→0`: a genuine theorem (convergent OS/SZZ expansion + a finite cubic-group lemma). But **the two lightest channels are exactly degenerate at leading order** (both are single-plaquette operators, `a_c=1`), so *which* of `0⁺⁺,2⁺⁺` is lightest — and its isolation — is a subleading `O(1)` dynamical splitting no symmetry fixes. Thus even at strong coupling "the lightest is `0⁺⁺` and isolated" is **not** a theorem; the lightest object is the degenerate *pair*. And the proven membership does not cross the analyticity wall to `β≈2.5`; persistence stays empirical.

## The theorem that does hold (strong coupling, weak form)

Leading strong-coupling asymptotics of the convergent character/cluster expansion give the zero-momentum mass
```
m_c(β) = a_c · ln(1/u(β)) + O(1),   u(β) → 0⁺ as β→0,
```
with `a_c` the minimal plaquette content per unit *time* of the lightest tube interpolating channel `c`. The single-plaquette operator, under the cubic group `O_h`, spans `A1g ⊕ Eg ⊕ T2g` — i.e. exactly `0⁺⁺ ⊕ 2⁺⁺`; **every** other `J^{PC}` needs `a_c ≥ 2` (`0⁻⁺` requires the `ε`-contraction `tr FF̃` of two field strengths; C-odd needs ≥3 gluons; `J≥3` needs transverse width). Hence
```
m_c − m_{0⁺⁺} ≥ (a_c − 1)·ln(1/u) + O(1) → +∞   (c ∉ {0⁺⁺,2⁺⁺}),
```
and, by the Osterwalder–Seiler / SZZ convergence (`β<1/48`, Hooft-normalised, λ>96), there is `β₁>0` such that for all `0<β<β₁`
```
inf_c Δ_c = min(Δ_{0⁺⁺}, Δ_{2⁺⁺}).
```
This **proves (b′)** at strong coupling: no `0⁻⁺`, odd-C, or spin-≥3 state is the lightest. It is the first non-empirical statement on IN-11 in the dossier — and it is sharp: the bound is an exact leading-order coefficient gap, not a loose constant.

## Why the core ordering fails even here

`a_{0⁺⁺}=a_{2⁺⁺}=1`. The two anomaly-protected channels share the **same** minimal flux tube (one plaquette per time slice), so `m_{2⁺⁺}/m_{0⁺⁺} → 1` and continuity/analyticity deliver nothing. The ordering and isolation of `0⁺⁺` live entirely in the `O(1)` subleading splitting — a finite strong-coupling computation in principle (finite order + OS remainder), but (a) not carried out to the needed uniformity, and (b) even if done, confined to `β<β₁`. At strong coupling the lightest object is the **degenerate pair** `{0⁺⁺,2⁺⁺}`, not an isolated `0⁺⁺`; isolation *across* these two channels is precisely the missing subleading datum.

## Persistence wall

`β₁` sits at the analyticity radius (`λ>96`), ~13× in `λ` below the crossover `λ≈7.3` (`β_W≈2.47`, SU(3); ATOM-strong-coupling-extend). `m_c/m_{c'}` is analytic only on `(0,β₁)`; nothing crosses Savvidy/crossover. The `β≈2.5` ordering (`0⁺⁺` 1710 < `2⁺⁺` 2390 < `0⁻⁺` 2560 MeV; Morningstar–Peardon) — even in its weak-form, matching the strong-coupling membership — is a controlled *extrapolation*, not a theorem. No rigorous continuation exists.

## (ii) no new symmetry; (iii) OZI orders nothing

The 1-form center `Z_N^{(1)}` charges Wilson **lines**; every glueball is N-ality-0 → blind (projected adjoint ≡ 1). Positivity/Källén–Lehmann is intra-sector: at `θ=0`, `P,C` make every mixed two-point kernel vanish identically, so no cross-channel inequality exists (T1; C1; GAP-IN11). Planarity/large-N: both channels single-trace, `m∼N⁰`, identical genus structure, and no controlled 4D `1/N`. **None orders the degenerate pair.**

**Sharpest line.** Strong coupling yields *membership* `{0⁺⁺,2⁺⁺}`, never the *order within it*: the two channels carry the same minimal flux tube and are leading-order degenerate, so "lightest = `0⁺⁺`, isolated" is a subleading dynamical splitting invisible to every symmetry, anomaly, and leading-order count — and the proven membership stops at the analyticity wall, short of `β≈2.5`.

*Refs.* Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Shen–Zhu–Zhu, arXiv:2204.12737 / CMP **400** (2023) 805; Münster–Weisz (strong-coupling convergence); Morningstar–Peardon, hep-lat/9901004; `ATOM-strong-coupling-extend`; `GAP-IN11`; `C1-derive-INb`; `T1-INb-symmetry`; `IN-3-attack`.
