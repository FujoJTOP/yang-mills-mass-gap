# MONOTONE-COUPLING — can a monotone stochastic coupling to a proven-gapped reference lower-bound the β=2.5 gap?

> 2026-11-23 · Attack of the **comparison-to-a-gapped-reference** route: `μ_{2.5} ≼ ν_ref` (or `≽`), `ν_ref` gapped with a **proven** gap, `gap` monotone under the order ⇒ `gap(2.5) ≥ gap(ref) > 0`. **判词 (verdict): NO-GO, 2/5.** The route fails on **(b)**, direction, independently of **(a)**, existence, and **(c)**, survival. Every gapped reference (strong-coupling lattice; mass-deformed) is gapped **exactly where its gap is *larger* than at 2.5**: the gap in lattice units is **decreasing in β** (AF scaling `√σ ~ √(−log(β/2N²))`, `gap(strong) > gap(2.5)`) and **increasing in the deformation `m_0`** (`gap(YM)=inf_{m_0}`). A monotone comparison therefore yields `gap(2.5) ≤ gap(ref)` — an **upper** bound, the reverse of the target. **No `5/5`.**

## 1 · The exact requirement
A lower bound needs `gap` monotone **downward** toward the physical point: `β↑` (strong→2.5) must **not increase** the gap, and `m_0↓→0` must **not increase** it, with a `β`/`m_0`-uniform floor `>0`. The bounded quantity must be the gap of the **less-gapped** (physical) end; the proven reference must sit on the **more-gapped** end and dominate it *from below*. No named reference does.

## 2 · (b) The direction wall — the decisive failure
The physical theory at 2.5 is the **spread-out end** of the order: crossing from strong coupling lengthens the correlation length, and in lattice units the mass shrinks (`m·a→0` toward the continuum). Stochastic order measures **concentration**; the gap is large where the measure is concentrated (strong coupling) and small where it is spread (2.5). So order and gap move **oppositely**: any correct `μ_{2.5} ≼ ν_strong` reads "physical is the less-localized, smaller-gap side," which is precisely the side you cannot lower-bound from a larger-gap reference. Mass deformation is the same arrow: `(½)m_0²|A|²` **lifts** the spectrum, `gap(m_0)` ↑ in `m_0`, so `gap(YM)=lim_{m_0→0}gap(m_0)` is the **infimum** and every proven `gap(m_0>0)` is an **upper** bound. This is the gap-valued twin of PROOF-TRANSPORT-MONOTONE: the monotone objects are *totals*, the atom-coercive ones are Schur-convex — opposite sides of the order.

## 3 · (a)+(d) Order exists or not; FKG ≠ Ginibre
**FKG/positive association fails on SU(N)** (compact ⇒ torsion ⇒ no group-compatible lattice order; PROOF-COMPARISON-FKG §1), so the *easy* Strassen order is unavailable — but this is **not** the binding obstruction. A rigorous monotone inequality **does** exist without FKG: **Ginibre** (1970) gives `⟨W⟩` and strip correlators monotone in `β` via character positivity (the same structure OS'78 use). It orders the **wrong object** (means/connected correlators — anti-decay, FKG-sign) and, by §2, the **wrong way** (σ decreasing in `β`). Order-free couplings escape the order entirely and collapse to the wall: the **Wasserstein/reflection coupling** has drift rate `≡G″=H1` (D2, von Renesse–Sturm sharp); the **Dirichlet-form comparison** (Deuschel–Stroock) finds every gapped reference form has the *opposite* scale in its kernel, gap `=0` on `L²(μ_Λ)` (PROOF-DIRICHLET-COMPARISON). So (d): no FKG blocks the easy order but nothing else — the surviving orders are `H1` in costume.

## 4 · (c) References and survival
- **Strong-coupling lattice** (OS'78; SZZ `β<1/48`; Borgs–Seiler `β<12/17≈0.706`): *proven*, but the cluster radius `β₁=O(1/N)≈1/48` is `~10²` below 2.5. Analytic β-continuation (PROOF-NO-MELT-PHASE: mass is β-analytic) carries the **value**, not a **sign-controlled uniform bound** — the ATOM/extend + PIROGOV-SINAI §2(B) wall.
- **Mass deformation** (steelman): a massive reference `gap≥m_0` is trivially proven. But "surviving `m_0→0`" is **circular** — a uniform-in-`m_0` floor is the YM gap.
- **Center-vortex / abelian `Z_N`**: `Z_N` gauge has FKG and a strong-coupling gap, but abelian projection is a **gauge-fixed map, not a measure coupling** (no domination of full `SU(N)`; Singer/Gribov), and its payout is `σ>0` in the **1-form channel** (`√σ`, spin-1), not the `0⁺⁺` residue (PROOF-CENTER-VORTEX; ABELIAN-PROJECTION).

## 5 · Acid test (letter PASS, spirit = the wrong sign)
Free massless field is the `β→∞` gapless end. If the route "worked" it would falsely certify a gap; instead it returns `gap(free)=0` (the strong-coupling `≤` is trivially true) — it is an **upper-bound machine**, exactly GIRSANOV's "letter-PASS, sign-wrong" profile. Non-vacuous, but aimed at the opposite side.

## Five-point scorecard
1. **Nonperturbative — PASS.** Coupling comparisons are measure-level; no loops.
2. **Sign-controlling `m_eff²(0)≥λ` — FAIL.** Direction wrong: gap ↓ in `β`, ↑ in `m_0`; reference gap is an *upper* bound.
3. **Non-circular — FAIL.** Reference must be gapped at 2.5 (=gap); mass `m_0→0` survival (=gap); order-free couplings = H1 (`G″`/Dirichlet).
4. **Non-SUSY + fails free field — PASS** (returns 0; no spurious gap).
5. **Uniform `≥λ` over `[1/48,2.5]` — FAIL.** No reference survives to 2.5; `β₁≈1/48`; deformation = infimum.

**2/5 — NOT `5/5`.**

## Exact failing estimate
```
target:  gap(2.5) ≥ gap(ref)>0  via  μ_{2.5} ≼ ν_ref  and  gap monotone toward physical
strong :  gap(β)  ~  √σ(β),  σ a²=−log(β/2N²) ⇒ gap DECREASING in β ;  gap(2.5) ≤ gap(1/48)   [UPPER]
mass   :  gap(m₀) ↑ in m₀           ⇒ gap(YM)=inf_{m₀};  gap(m₀>0) ≥ gap(0)                     [UPPER]
radius :  cluster β₁=O(1/N)≈1/48 ≪ 2.5 ; analytic continuation ≠ sign-controlled bound
order  :  FKG absent (torsion); Ginibre orders ⟨W⟩ (means, wrong object);
          D2 drift ≡ G″=H1 ; Dirichlet-form gap on L²(μ_Λ)=0 ; free field → 0 (wrong sign)
```

**Sharpest line.** The route needs a monotone order in which the **physical** theory is on the *dominated* side and the gap points **down** to it — but every proven-gapped reference sits on the **more-concentrated, larger-gap** side, and the gap has the **wrong monotone sign** (it *shrinks* toward 2.5 and *grows* with the deformation). Absence of FKG is a decoy: the rigorous order that survives (Ginibre) still orders *means*, and the order-free couplings are `H1`/`H*` re-spelled. Granting a perfect monotone coupling converts a strong-coupling gap into an **upper** bound on the physical one — the reverse of Clay. The missing object is not an order but a **uniform floor at the spread-out end**, which *is* the gap.

*Refs.* PROOF-COMPARISON-FKG (no FKG on SU(N); wrong inequality type; `C(β)≈15`); PROOF-DIRICHLET-COMPARISON (Deuschel–Stroock; gapped reference gap=0 on `L²(μ_Λ)`); D2-COUPLING-DOEBLIN (drift rate `≡G″=H1`, von Renesse–Sturm sharp); PROOF-TRANSPORT-MONOTONE (monotone=total vs atom-coercive=Schur-convex, opposite sides); PROOF-CENTER-VORTEX / ABELIAN-PROJECTION (`σ>0` is 1-form/`√σ`, cross-channel; gauge-fixed); PROOF-P3-UNCONDITIONAL / SZZ `β<1/48`; Borgs–Seiler, CMP 91 (1983) 329; PROOF-NO-MELT-PHASE (analytic value ≠ sign-controlled bound); ATOM-strong-coupling-extend; PIROGOV-SINAI §2(B); GIRSANOV (letter-PASS/sign-wrong precedent); GOLDS-POS-MASS (no positive-mass atom from positivity). External: Ginibre, *CMP* **16** (1970) 310; Harris 1960; Fortuin–Kasteleyn–Ginibre 1971; Strassen 1965; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Deuschel–Stroock, *Large Deviations* (1989) Ch. 5; Holley–Stroock, *J. Stat. Phys.* **46** (1987) 1159; Tomboulis–Yaffe, *CMP* **100** (1985) 313; Shen–Zhu–Zhu, *CMP* **400** (2023) 805.
