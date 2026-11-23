# DEFECT-PEIERLS — the center-vortex / monopole Peierls bound is rigorous only in the **dilute** direction, yields the **1-form `√σ`** where it bites, and its "same-channel `0⁺⁺`" avatar is an identity (compressibility `=` susceptibility `=` H1), not an estimate

> 2026-11-23 · NEW SPECIES: topological-defect (center-vortex / monopole) condensation via a rigorous **Peierls/contour bound**, targeting `G″(k→0)≥λ` in the **same** `0⁺⁺` channel as the gap. **判词: 1-form only.** The Peierls estimate bounds the *suppression* of large defect surfaces — the **dilute / non-condensed** direction — via `P(surface)≤e^{−σ_v A}` (`σ_v>0`). Condensation (the confinement-producing statement) is its **opposite**, `σ_v≤0`, and needs a **positive attractive** contour measure; non-abelian SU(N) has none (signed Weingarten, no FKG). Where the sign *is* controlled (`β<β₁` / `Z_N` strong coupling) area law **and** gap arrive **together** (=P3), never as a vortex-condensation theorem, never in `0⁺⁺`. Five-point: **#1 PASS(regime), #2 FAIL, #3 FAIL, #4 PASS(genuine), #5 FAIL ⇒ 2/5.**

## §1 What a Peierls bound delivers
Peierls pits entropy (`≈e^{cA}` surfaces) against energy (`e^{−σ_v A}`). It is **one-sided**: the sign-controlled statement is `P(large defect)≤e^{−σ_v A}`, `σ_v>0` — defects are **dilute** (non-condensed). Condensation (`⟨W⟩` area law) needs the opposite, `σ_v≤0`, i.e. an **attractive** positive fugacity. `Tomboulis–Yaffe` turn "the vortex is cheap (`F_v` perimeter-law)" into a *sufficient* `σ_W≥cρ̂>0` — but "the vortex is cheap" **is** confinement restated. The controllable regime is strong coupling `β<β₁` (SZZ `β<1/48`; Borgs–Seiler `0.706`), where Osterwalder–Seiler already give area law ∧ clustering from one KP expansion — P3 renamed.

## §2 Death 1 — channel: **1-form, not `0⁺⁺`**
The defect order parameter is the magnetic `Z_N^{(1)}` 't Hooft loop; the area law is electric-`Z_N^{(1)}` confinement — the **flux-tube** sector, gap `√σ`. `G″(k→0)=1/Var(a₀)` is the neutral **local** `tr F²` (`0⁺⁺`, 0-form): a point carries no 1-form charge, so every local conserved current is 1-form-neutral and no Ward/Peierls identity transports `σ` into the `0⁺⁺` residue (`PROOF-1FORM-ANOMALY`). The strongest same-channel candidate — the **vortex-density** correlator — dies by *gauge invariance*: it is defined only after abelian projection (Singer; Elitzur), so not the physical `tr F²`; gauge-invariant, it is the non-local 't Hooft loop (1-form) again. Its susceptibility `∫⟨ρ_vρ_v⟩_c=Var(N_v)` is an **identity** (=H1), not a lower bound. **1-form only.**

## §3 Death 2 — sign/phase dichotomy
Rigor ⟺ dilute (small signed weight) **or** strong-coupling P3. The **condensing** direction (`σ_v≤0`) needs a positive measure; SU(N) has none (torsion ⇒ no FKG; signed Weingarten — `C8`, `PIROGOV-SINAI` §2). At `β=2.5` the fugacity is `O(1)`, the Peierls expansion does not converge, and the Savvidy reference `−D²+2adF` is indefinite. "Vortices condense at 2.5" ⟺ "confinement at 2.5" — input = answer, circular (`PROOF-CENTER-VORTEX` §3).

## §4 Free-field acid test — genuine PASS
The free massless 4D field has **no finite-action defect to condense**: zero vortex/monopole fugacity, no Peierls gap, no area law, correctly **no mass** (`f₁²=0`). The tool is neutral on the free field — it fabricates no gap, so it is **non-vacuous**. Non-SUSY ✓. Necessary, not sufficient: one nonzero fugacity (interacting) yields a *positive* Peierls cost — exactly the wanted discriminator.

## §5 Five-point
1. **Nonperturbative — PASS, regime-limited** (`β<β₁` / Borgs–Seiler / `Z_N`; crossover open).
2. **Sign-controlling LOWER bound `G″≥λ` — FAIL.** Bounds `P(contour)` (1-form, an *upper* bound on `⟨W⟩`); the target is a *lower* bound on the neutral correlator — wrong sign, wrong channel.
3. **Non-circular — FAIL.** Condensation ⟺ area law ⟺ confinement (=P3's `β₁`); the bridge to `G″` is `m_{0⁺⁺}=c√σ` (measured, not theorem).
4. **Non-SUSY + fails free field — PASS (genuine, §4).**
5. **Uniform `λ` — FAIL.** Radius `β₁`/`0.706`; `β=2.5` non-convergent; no `0⁺⁺` `λ`.

## Exact failing estimate
```
Peierls:   P(surface)≤e^{−σ_v A}, σ_v>0  ⇒ DILUTE = non-condensed  [sign-controlled]
Condense:  σ_v≤0 (F_v≤perim) ⇒ ⟨W⟩ area law — needs POSITIVE attractive measure [SU(N): SIGNED]
radius:    β<β₁ (SZZ 1/48; Borgs–Seiler 12/17) ⇒ area law∧gap = P3 renamed
β=2.5:     fugacity O(1), no convergence; Savvidy Hess −D²+2adF indefinite ⇒ no reference
channel:   't Hooft loop = 1-form; a point carries no 1-form charge ⇒ tr F² blind ⇒ no 0⁺⁺ transport
same-ch:   vortex density gauge-dependent (Singer/Elitzur); susceptibility = Var(N_v) = H1 (identity)
free:      zero defects, σ_v undefined, no gap ⇒ acid PASS
```

**Sharpest line.** The Peierls/contour bound is rigorous exactly in the **dilute** direction — large defects are *suppressed*, the **non-condensed** phase — whereas **condensation**, which would force the area law, is the *opposite* inequality and demands a **positive attractive** fugacity SU(N) lacks (signed Weingarten, no FKG). Where the sign *is* controlled (`β<β₁`, Borgs–Seiler, `Z_N`) area law arrives **with** the gap (Osterwalder–Seiler) — P3 in disguise — and everywhere it is the **1-form `√σ`**, blind to the neutral `tr F²`. The same-channel hope dies twice: the vortex density is not gauge-invariant (not `tr F²`), and its gauge-invariant avatar ('t Hooft loop) is 1-form. **确切断于 #2（错道/错号，1-form 上界）与 #3（凝聚=禁闭，循环），#5 半径止于 β₁；#1、#4 通过而不触正界。`G″(k→0)≥λ` 未达，墙未动。**

*Refs.* Peierls, PCPS **32** (1936) 477; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Borgs–Seiler, CMP **91** (1983) 329; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; 't Hooft, NPB **138** (1978) 1 / **153** (1979) 141; Tomboulis–Yaffe, CMP **100** (1985) 313; Elitzur, PRD **12** (1975) 3978; Singer, CMP **60** (1978) 7; Greensite, Prog. Part. Nucl. Phys. **51** (2003) 1; Gaiotto–Kapustin–Seiberg–Willett, arXiv:1412.5148; Morningstar–Peardon, PRD **60** (1999) 034509. Internal: PROOF-CENTER-VORTEX; PROOF-CONFINEMENT-RANDOM-SURFACE; AREALAW-GAP; PIROGOV-SINAI; C8-SPECIES-NOGO; FORGE5-AREALAW-GAP-LINK; FORGE5-CONFINING-MONOPOLE; PROOF-1FORM-ANOMALY; PROOF-STRING-TENSION-STABILIZER; PROOF-P3-UNCONDITIONAL.
