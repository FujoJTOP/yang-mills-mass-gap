# B4-VARIANCE-MAG — the magnitude route is NOT small: `Var(G|0)≈2` is equipartition-pinned and L-stable (2.01→2.05), the CS threshold is 0.8–1.3, so Cauchy–Schwarz cannot certify even the *present* gap; and "make it <1" under RG *is* H1

> **判词 (verdict).** Round 4 attacks the **magnitude** `Var_μ(⟨v,∇S⟩|0)`, not its sign. To certify `G″>0` via Cauchy–Schwarz needs `Var(G|0)<W_geo²/(4VarΦ′)≈0.8–1.3`; measured is **≈2.0**, and all three routes fail informatively. **(a)** No RG mechanism pushes it below 1 — the force variance is pinned at `≈E[H]≈2` by **equipartition** (a UV/marginal, not IR, property); an RG argument that *made* it <1 would, by the exact identity, *be* the emergence of IR stiffness, i.e. H1. **(b)** The exact `E[f²]=E[H]` route is gap-blind (free field saturates it with `Var(a)=∞`), its conditioning correction being the whole gap (SUB1). **(c)** The L-scan shows `2.01` (L=4) → `2.05` (L=6): **volume-stable**, not a finite-L artifact shrinking below threshold. Deeper: CS bounds `|X|` by the *magnitude* `√(VarG·VarΦ′)`, discarding exactly the near-zero correlation (`corr(a,∇S)≈0`) that makes `X≈0.04` and carries the sign — it replaces `corr=0` by worst-case `±1`. So the magnitude route is not merely H1-renamed; via CS it is **sign-blind and numerically short (2 vs 1) even where the gap is present** (`G″≈0.4>0`).

## §1 Threshold is below the physical value
Exact identity: `G″=Δ_dyn+W_geo+X`, `Δ_dyn=E[H|0]−Var(G|0)`, `X=−2Cov(G,Φ′|0)`, `W_geo≥N/12`. CS: `|X|≤2√(VarG·VarΦ′)`. With `Δ_dyn≈0` (equipartition), `VarG≈2`, `VarΦ′≈0.02`, `W_geo≈0.32`: `|X|≤0.40>W_geo ⇒ G″≳−0.08<0`. Closing needs `VarG<W_geo²/(4VarΦ′)≈0.8–1.3`. **The measured 2 sits 1.5–2.5× above threshold** — CS is not "one inequality short," it is short *where `G″` is already positive*.

## §2 Route (a): RG/flow decrease — fails, circular
`Var(G|0)≈E[H|0]` (equipartition, measured `E[H]−Var∈[−0.06,+0.02]`). `E[H|0]` is the per-mode action stiffness, `O(2)`, *UV/marginal* — AF suppresses `g`, not `E[H]`. For `Var(G|0)` to fall below 1 under coarse-graining the effective stiffness must *grow* as `k→0`, i.e. an IR curvature must appear — that growth **is** `G″(k→0)≥λ` **is** H1. Equipartition blocks the Gaussian alternative (small stiffness ⇒ small `VarG` but small `G″`, no sign).

## §3 Route (b): `E[f²]=E[H]` + conditioning — fails (SUB1)
The identity is exact for every measure, any β — and **gap-blind**: massless `k=0` has `H=0`, so `Var(G|0)=0≤C` trivially while `Var(a)=∞`. Any bound `Var(G|0)≤C` is therefore *satisfied by the free field* (vacuous) and does **not** fail the free-field acid test. The conditioning correction `D=Var(f|0)−E[K⁻¹c,c]` is `Var(G|0)`; bounding it above *is* H1.

## §4 Route (c): L-scan — the 2 is real
Per-mode `VarG`: L=4 `{2.12,1.87,2.05,2.09,1.96,1.98}` (mean **2.01**); L=6 `{1.91,1.86,2.57,2.33,1.54,2.09}` (**2.05**). **No L-decrease**; `E[H]` flat (`2.03`/`1.99`). Finite-L artifact hypothesis disfavoured. Caveats that do *not* help: no run reaches true `k→0` (`π/2,π/3`), single β=5.7 — but if `VarG` dropped at `k→0`, that drop would *be* the gap, not a free lunch.

## §5 Does it separate from H1? No — two ways
**(i) Defect reading:** `Var(G|0)−E[H|0]<W_full` is (◆), verbatim H1. **(ii) Magnitude/CS reading:** needs `VarG<1`, but `VarG≈2` is equipartition/UV-pinned and L-stable, so no honest uniform bound is <1, and CS cannot certify a gap that *exists* (`G″≈0.4`). CS is the wrong functional — it uses `√(VarG·VarΦ′)` (worst-case `corr=±1`) while the physics is `corr(a,∇S)≈0`; the sign lives in near-orthogonality, a magnitude-invisible fact. **The magnitude route is strictly worse than H1: circular read as the defect, numerically dead read as the magnitude.**

## Five-point verdict
1. **Nonperturbative — PASS** (exact identity/RG, any β).
2. **Sign-controlling = LOWER bound — FAIL.** `VarG≤C` is an *upper* bound; CS converts it only at `C<1`, not met (`2>1`), so `G″≥−0.08` indefinite.
3. **Non-circular — FAIL.** `VarG<1` uniform `⟹` IR stiffness `⟹` H1; L-scan no rescue.
4. **Non-SUSY + FAILS free field — FAIL (vacuous).** Free field `Var(G|0)=0≤C` satisfies the bound; no YM positivity extracted.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Magnitude O(2), threshold O(1); no uniform `λ`.

## Exact failing estimate
```
need  Var_μ(G|0) < W_geo²/(4VarΦ′) ≈ 0.8–1.3 ;   measured VarG ≈ 2.0  (L4 2.01, L6 2.05 — L-stable)
CS:   |X| ≤ 2√(VarG·VarΦ′) ≈ 2·√(2·0.02) ≈ 0.40 > W_geo ≈ 0.32  ⇒  G″ ≳ −0.08 < 0   (yet true G″ ≈ +0.4)
RG:   VarG↓<1 under blocking ⟺ coarse stiffness ↑ at k→0 ⟺ G″(0)≥λ  = H1.
free: Var(G|0)=0  ⇒ bound vacuous  ⇒ does NOT fail the free-field acid test.
CS uses magnitude √(VarG·VarΦ′) (corr=±1) while physics is corr(a,∇S)≈0 ⇒ sign-blind.
```
**Sharpest line.** The coarse force variance is emphatically *not* small: `≈2`, equipartition-set, stable L=4→6. Cauchy–Schwarz certification needs `<1`, contradicted by data at a point where the gap is already measured positive (`G″≈0.4`); so the magnitude route cannot close, and any RG argument that *made* the variance <1 would be the emergence of IR stiffness = H1. CS is the wrong functional — bounding the cross term by the magnitude throws away the `corr(a,∇S)≈0` near-orthogonality that carries the sign. **The magnitude route does not separate from H1; read as a defect it *is* H1, read as a magnitude it is sign-blind and numerically short 2×.**

*Refs.* C4-COV-LASTMILE (§5); SUB1-COV-BOUND (gap-blind identity, `D`=gap); ABC2-CONSTRUCTIVE (§4 0.40>0.32); COMPUTE-GPP-SU3 (equipartition); PROOF-JACOBIAN-SAVES ((◆)); OTHER-GPP-KZERO ((◆)=H1); CONTINUUM-VARIANCE-DEFECT (`D/S` RG-covariant); INEQ-R7 (constancy, not edge); w7/su3_L4.log, su3_L6.log, k0_continuum_analyze.py (`V_own.G`); Brascamp–Lieb 1976; Prékopa 1971.
