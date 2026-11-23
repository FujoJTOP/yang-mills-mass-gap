# FORGE5-SL-SCALE-SUM — summing the backward cost over RG scales: does UV suppression converge the SL charge to O(1)?

> 2026-11-23 · FORGE5 direction: attack constraint #5 of HARDFORGE-SL-FULL by summing the per-scale backward cost `osc=Σ_b osc_b` explicitly with the UV asymptotics. **Verdict (判词): 4/5 stands — the sum does NOT converge to O(1) at β=2.5.** The proposed mechanism (running `g_b→0` + trap emptying the Savvidy cone) is **real but is a deep-UV artifact**. The sum is UV-convergent — `O(1)` in `(L,a)`, hence uniform #5 — *iff* the cone edge `κ_b=√(2g_b|F_b|)` falls below the shell lower edge at the cutoff, which needs `g_b<1/(√2e)≈0.26`, i.e. `1/g²>2e²≈14.8`, i.e. `β>59` (SU(2)). At β=2.5, `1/g²=β/2N=0.625`, the cutoff coupling is `g_J≈1.27`, `κ_J≈1.8>Λ_{J−1}=π/e≈1.16` — the UV shell stays **partially concave**. Concave modes persist to the Brillouin edge (W6: 22–49% negative, 90–97% high-momentum), so `N^{conc}_b=Θ(N_□)` and `Σ_b osc_b=Θ(N_□)`; cost `e^{Θ(N_□)}`. **Secondary honest finding: the burden is bulk/mid-to-UV, not the k→0 IR mode** — the k→0 scale carries a single O(1) charge, never the Θ(N_□) bulk; and H1's k→0 sign (`G″_∞<0`, W6 −0.165) is a separate defect no osc bound touches.

## §0 Setting and the explicit sum
Shell `b` (`b=0` IR, `b=J` UV), `Λ_b=Λe^b`, `Λ_J=1/a`, `J=log(L/a)`; mode count `N_b≈VΛ_b⁴/(2π)⁴∝e^{4b}` (fixed physical `V`, rising toward UV), `Σ_b N_b=Θ(N_□)`. Exact Savvidy data (PROOF-SAVVIDY-EXACT-SIGN): depth `δ_b:=−λ_min|_b=2g_b|F_b|`; cone `{k²<δ_b}`, edge `κ_b=√δ_b`. Each *concave* mode carries an O(1) compactness-regulated charge (HARDFORGE §4: `∫_0^{t_b}½dτ/(λ_b+τ)`, pole at `τ=σ_b`, `|x|≤π`). Hence, with `N(<κ)=κ⁴V/(32π²)`:
```
osc_b ≃ c·N^{conc}_b ,   N^{conc}_b = N(<min(Λ_b,κ_b)) − N(<Λ_{b−1}) .        (S)
```
(The task's tighter form `osc_b∝t_b·N^{conc}_b·δ_b` multiplies by `t_bδ_b=O(1)` at β=2.5; it changes no verdict below.)

## §1 Three regimes; the convergence criterion
Since `Λ_b=eΛ_{b−1}`:
- `κ_b>Λ_b`: whole shell concave, `N^{conc}_b=N_b∝e^{4b}`;
- `Λ_{b−1}<κ_b<Λ_b`: partial, `N^{conc}_b=N(<κ_b)−N(<Λ_{b−1})`;
- `κ_b<Λ_{b−1}`: shell **empty**, `N^{conc}_b=0`.
Convergence to `O(1)` in `(L,a)` requires regime (iii) to set in at a `b*` **fixed** as `L` (hence `J=log(L/a)`) grows — i.e. `κ_b<Λ_{b−1}` for `b>b*`, `b*` independent of `J`.

## §2 The conditional 5/5 branch (deep UV)
Grant the AF scaling `g_b²=1/(2β₀b)` and the *physical* field scaling `|F_b|~g_bΛ_b²` (field at a fixed fraction of `Λ_b²`). Then `κ_b=√2 g_bΛ_b`, and `κ_b<Λ_{b−1}=Λ_b/e ⟺ √2g_b<1/e ⟺ b>b*:=e²/β₀`, a **fixed** shell count (SU(2): `b*≈159`), independent of `(L,a)`. Concave modes confined to `b≤b*`; `osc=Σ_{b≤b*}cN(<κ_b)∝Σ g_b⁴e^{4b}~e^{4b*}=e^{Θ(1/β₀)}`. This is **`O(1)` in `(L,a)`** — volume-independent, so #5's uniformity *would* hold (a genuine, if microscopic, gap `m²≳λe^{−e^{Θ(1/β₀)}}`). **This is the real 5/5 branch — and it lives at β≳59.**

## §3 Why β=2.5 is outside that branch (the 4/5 reality)
§2 needs (a) asymptotic running *at the cutoff* and (b) fields saturating `~g_bΛ_b²`. Both fail at β=2.5. (a) `1/g²=0.625`; over the available `J≲7` shells the running moves `1/g²` by only `2β₀J≈0.09J`, leaving `g_J≈1.27` — far from the emptying value `0.26`. Reaching it needs `Δ(1/g²)≈14.2`, i.e. **~152 more shells** (`L/a≳e^{150}`) or **β>59**. (b) The measured β=2.5 field is equipartition `|F|~O(g)`, not `~g_bΛ_b²`; so `κ_J=√(2·1.27·|F|)≈1.6–2.0 > Λ_{J−1}=1.16`. **The trap's stated bound `κ_b<Λ_b` (1.8<3.14) is satisfied but is the wrong inequality** — emptying a shell needs the strictly stronger `κ_b<Λ_{b−1}=Λ_b/e`. This factor-`e` gap is the precise defect: the trap lives on the *retained band* `|k|≥Λ_b`, never on the shell `[Λ_{b−1},Λ_b]` itself. Consequently concave modes survive to the cutoff at O(1) fraction (W6), `N^{conc}=Θ(N_□)`, and
```
Σ_b osc_b = c·Θ(N_□) = Θ(N_□)  ⟹  cost e^{Θ(N_□)}.        (DIVERGES)
```

## §4 Does the k→0 IR scale carry the cost? No — the bulk does
The k→0 mode is a *single* mode: one O(1) charge. The Θ(N_□) of §3 is carried by mid/UV shells (`N_b∝e^{4b}` maximized near the cutoff, matching W6's 90–97% high-momentum). So the cost is **volume-extensive and bulk**, not IR-localized at "the few IR scales." Equally important, this decouples `osc` from H1: H1 is the sign `G″_∞(k→0)≥λ`, measured **negative** at β=2.5 (W6 −0.165) — a defect in the coarse Hessian-vs-variance that **no** oscillation bound repairs. The IR scale fails H1 on its own; the bulk fails `osc`.

## Five-point verdict (honest)
1. **Nonperturbative — PASS.** Exact `Hess+tI`, BL, Holley–Stroock; no `g`/`1/N` expansion.
2. **Sign control on the tilted cone — PASS.** `Hess S+T⪰ε` value-free.
3. **Non-circular — PASS (structure).** `δ_b` is bare curvature.
4. **Non-SUSY + fails free field — PASS.** Free field `δ_b≡0⇒osc=0`.
5. **Uniform `m²(k→0)≥λ` — FAIL at β=2.5.** The scale sum is **not** UV-convergent: `Σ_b osc_b=Θ(N_□)` because the cutoff coupling is O(1) (β=2.5 ≪ the β≳59 emptying threshold) and the concave bulk survives to the Brillouin edge. **Exact failure:** `κ_J≳Λ_{J−1}` at the cutoff — the trap bound `κ_b<Λ_b` is the wrong inequality (need `<Λ_{b−1}=Λ_b/e`), so the proposed UV suppression never engages in the physical regime. Convergence to `O(1)` is real but confined to β≳59, not β=2.5.

*Refs.* HARDFORGE-SL-FULL (§4 `e^{osc}`, `osc≲½π²Σt_bN_b`, per-mode compactness charge); PROOF-MULTISCALE-4D (§2 trap `√(2g_j|F_B|)<Λ_j` ⇒ positive only on `|k|≥Λ_j`); PROOF-SAVVIDY-EXACT-SIGN (`−g|B|`, `N_neg∝(gB)²`); W6-BE-at-2.5 (−1.12β, 22–49% neg, 90–97% high-momentum, `G″(k→0)≈−0.165`); W4-polchinski-target; PROOF-COARSE-POINCARE (H1); Gross–Wilczek/Politzer (AF, `β₀=11N/48π²`).
