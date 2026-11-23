# ABC2-NEW-INEQ — shell-conditional Brascamp–Lieb transport (STC): the UV shells close, the transport to k→0 is H1

> 2026-11-23 · SHAPE C round 2. Target: bound `Var_μ(⟨v,∇S⟩|0)` on the **typical-field, high-momentum, g=O(1)** shell — the tail WSTC (high-field, measure-starved) missed. Composite: **(a)** shell-conditional BL + AF closes the top shell; **(b)** transport/IR-isolation to `k→0`; **(c)** Ward–Stein second-moment non-Gaussianity. **VERDICT: a genuinely new composite (STC) and a true #1/#4 theorem — NOT 5/5, and the open cell after it is exactly H1.** STC proves `G″(k_b)≥Λ_{b−1}²−δ_b>0` on every **cone-free** shell (`δ_b<Λ_{b−1}²`), a correct lower bound — but at β=2.5 the sign-carrying shells (W6: 90–97% high-momentum) are exactly the concave ones (`δ_b>Λ_{b−1}²`) where STC is silent, and the transport to `k→0` is a **zero-sum**: block-diagonal ⇒ IR-isolated ⇒ H1 untouched; background-coupled ⇒ the indefinite inverse ⇒ H1. **#2 UV-only, #3 fails at closing, #5 fails.**

## 1. The inequality (boxed)
Shell `Π_b=[Λ_{b−1},Λ_b]`, `Λ_b=Λe^b`, `Λ_J=π/a`; probe `v_b` in `Π_b`, `f_b=⟨v_b,∇S⟩`. Shell Hessian `K_b=k²+2a\,d\,F_B(x)`, floor `K_b⪰Λ_{b−1}²−δ_b`, `δ_b=2g_b|F_b|` with `|F_b|` the **typical** (measure-dominant) field — not `sup_x` (fixes FORGE5-SL-UV-CONE-EMPTY K1) and not the high-field tail (fixes WSTC). Where `δ_b<Λ_{b−1}²`:
```
   STC:  Var_μ(f_b) ≤ (Λ_{b−1}²−δ_b)^{-1}  ⇒  G″(k_b) ≥ λ_b := Λ_{b−1}²−δ_b > 0.
```

## 2. (a) AF closes the top shell — only where the negatives are not
`δ_b<Λ_{b−1}² ⟺ 2g_b|F_b|/Λ_{b−1}²<1`. **Saturating field** `|F_b|~g_bΛ_b²` ⇒ ratio `2g_b²e²<1 ⟺ g_b<1/(√2e)≈0.26 ⟺ β≳59` (SU2, verified `1/g²>2e²=14.78`) — identical to FORGE5-SL-SCALE-SUM. At β=2.5, `g_J=1.27`, `δ_J=2.5–3.8>Λ_{J−1}²=1.34`, `κ_J=√δ_J≈1.6–1.9>Λ_{J−1}=1.16`: the top shell is **concave**. **Equipartition field** `|F_b|~O(1)` does force `δ_b/Λ_{b−1}²→0` on the *very highest* shells — but those are lattice-cutoff modes that **decouple as a→0**, carrying no physical mass (CONTINUUM-VARIANCE-DEFECT: AF kills `D`'s UV numerator, never the IR `D<S`). The physical negatives (W6, `k~1–2/a_phys`) are *low* momentum relative to `π/a`, where `g` is the IR coupling, not `g(π/a)→0`. **AF closes exactly the shells the continuum removes.**

## 3. (b) Transport to k→0 — a zero-sum
`G″(k→0)` needs `Var(f|0)` for a k→0 probe, not a per-shell sum. Two cases, no third:
- **Block-diagonal (Gaussian / translation-invariant):** `G″` diagonal in momentum; `G″(k→0)` is the k→0 entry alone, the high-shell `λ_b` are *other* entries. **Fast modes cannot gap the slow mode** (PROOF-IR-H1 §4, exact). IR-isolated ⇒ STC's UV bounds say nothing about H1.
- **Background-coupled:** the Savvidy `F_B` mixes k→0 with high k; `G″(k→0)` then needs the *inverse* of the full indefinite `K` (SCBL `G″=S−D`, `D>0⇔K` indefinite). That inverse **is** the gap. Coupling ⇒ circular.

The transport is `0` or `H1`: no regime certifies the k→0 curvature without the indefinite block.

## 4. (c) Ward–Stein second-moment — wrong direction
`E[af]=1 ⟺ G″Var(a)=1 ⟺` marginal Gaussian (ABC §1). Schwinger–Dyson `E[H]=E[f²]` is **global**; the slice `Var(f|0)` is not pinned by it (conditioning *is* the gap). Cramér–Rao `I(p)Var(a)≥1` (p the coarse marginal), equality iff Gaussian, measures non-Gaussianity by the slack `I−1/Var`. But H1 is `Var(a)≤1/λ` (upper bound); Cramér–Rao gives `Var(a)≥1/I` (**lower** bound). Bounding the slack below (forbidding the Gaussian edge) is a Fisher lower bound = log-Sobolev-hard = H1. Wrong-direction inequality.

## 5. Acid test / non-Schur / non-log-concave
Free massless: `F≡0⇒δ_b≡0`, every shell BL-closed (`Var≤1/Λ²`), yet the zero mode gives `Var(a)=∞⇒λ=0`. **No spurious λ** — fails free field as required (`f₁²=0`). Non-Schur (shell BL, never `H_ff^{−1}`, which is indefinite); non-global-log-concave (typical-field per shell, silent on concave shells — honest).

## 6. Five-point verdict
1. **Nonperturbative — PASS.** Per-shell BL + exact AF, no loop expansion.
2. **Sign-controlling LOWER bound — PARTIAL.** `G″(k_b)≥λ_b>0` genuine on **cone-free** shells (UV); **silent** on the β=2.5 sign-carriers.
3. **Non-circular — FAIL at closing.** UV closure clean; transport is IR-isolated (inert) or coupled (H1).
4. **Non-SUSY + fails free field — PASS.** `F=0⇒λ=0`, no positivity.
5. **Uniform m_eff²(k→0)≥λ — FAIL.** No λ at k→0; UV bounds are `G″(k_b)`, the wrong momentum.

**Sharpest line.** STC is a real sixth-species composite — typical-field, shell-conditional BL with AF-decoupled top — and it proves a true theorem: **every cone-free shell carries `G″(k_b)≥Λ_{b−1}²−δ_b>0`**, uniformly in L. But the gap's sign is decided on the shells where `δ_b>Λ_{b−1}²` (concave high-momentum at β=2.5), where STC is mute, and the transport to `k→0` collapses to the exact dichotomy of PROOF-IR-H1 §4: fast modes either decouple (cannot gap the slow mode ⇒ H1 untouched) or couple through the Savvidy background (⇒ the indefinite inverse = H1). **The open cell is exactly H1 — not smaller.**

*Refs.* ABC-NEW-INEQ (Ward, `G″Var=1`); FORGE5-CONE-VARIANCE (`E[H]=E[f²]`); FORGE5-SL-SCALE-SUM (`β>59` emptying); FORGE5-SL-UV-CONE-EMPTY (K1 sup-blind); FORGE5-SL-IR-SINGLE (k→0 = H1); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`, `D>0⇔K` indefinite); PROOF-IR-H1 (§4 fast≠slow); PROOF-COARSE-POINCARE (H1); PROOF-SAVVIDY-EXACT-SIGN (`k²+2adF`, `λ_min=−g|B|`); CONTINUUM-VARIANCE-DEFECT (UV numerator ≠ decision); W6-BE-at-2.5 (90–97% high-momentum); Brascamp–Lieb 1976; Cramér–Rao; Gross–Wilczek/Politzer.
