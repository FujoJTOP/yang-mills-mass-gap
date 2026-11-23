# FOURIER-UNIFORMITY — species: Weyl/equidistribution on the coarse mode is sign-reversed; "strong Fourier decay" is the gapless pole and the bandwidth is the Ward Fisher information

> **判词 (verdict).** Forge: bound the `k→0` marginal `p(a)`'s characteristic function `φ(t)=E_μ e^{ita}` by a Weyl/equidistribution estimate, `sup_{|t|≥δ}|φ(t)|≤η(δ)`, and read `Var(a)=−φ″(0)≤1/λ`, so `G″=S−D=1/Var−D≥λ`. **Only #1 passes; #2–#5 fail.** The lever is *sign-reversed*: the only inequality tying decay to variance is Stam (Heisenberg), `Var(a)·J(a)≥1` (equality iff Gaussian), a **lower** bound on `Var`, with `J` the marginal's Fisher information. "Strong (fast) Fourier decay" ⇒ narrow `φ` ⇒ small `J` ⇒ a large `Var` **floor** — the gapless pole. The free massless mode `φ=e^{−c t²/2k²}` has the *strongest* decay exactly at `k→0` where it is gapless: the species is not merely vacuous on the free field, it **falsely gaps it**. Equidistribution's real output (small `|φ|` off the origin = small bandwidth) *is* the gapless signal.

## 1. The forge, stated exactly
Chart `X=log U`, coarse mode `a=⟨v,X⟩`, marginal `p(a)`, `φ(t)=E_μ e^{ita}`. `E[a]=0` (Ward/Elitzur) ⇒ `Var(a)=−φ″(0)`; and `S=1/Var(a)`, `G″=S−D` (SCBL). Proposed lever: `a` is a linear functional of a sum of link phases whose single-link phases equidistribute (Peter–Weyl), so a Weyl bound makes `|φ(t)|` small off the origin — "strong decay ⇒ `p` near-uniform ⇒ `Var` small ⇒ `G″≥λ`."

## 2. The sign is backwards — Stam/Heisenberg
The unique inequality relating characteristic-function spread to variance is
```
Var(p)·J(p) ≥ 1,   J(p)=∫(p′)²/p dx = E_μ[(∂_a log p)²] = Var_μ(E[f|a]),
```
with equality iff `p` Gaussian. It makes `Var` the **lower**-bounded side: *narrow `φ` (strong decay, small `J`) forces a large `Var` floor.* For `p=N(0,σ²)`, `φ=e^{−σ²t²/2}`: the decay rate *is* `σ=√Var`. So fast Fourier decay tracks `√Var` and is *maximal* at the gapless IR point. Equidistribution (small `|φ(t)|`, fast decay) ⇒ narrow `φ` ⇒ small `J` ⇒ large `Var` ⇒ gapless. The species predicts a gap precisely where the variance is largest — the reversed pole. (Verified exactly: Gaussian `Var·J=1.0000000000`; mixture `Var·J=3.98>1`; Plancherel `∫t²|φ|²dt=2π∫(p′)²dx` to `6×10⁻¹⁴`.)

## 3. Acid test — fails the free field in the strong sense
Massless free: `p=N(0,c/k²)`, `Var=c/k²→∞`, `φ(t)=e^{−ct²/2k²}`. As `k→0` the decay rate `√c/k→∞`: the free massless marginal has the **strongest possible Fourier decay at the IR point and is gapless.** So "strong decay ⇒ gapped" is refuted by the free field itself, not vacuously satisfied. (The brief's "heavy = slow decay" is the *smoothness* convention, `\hat p(ξ)`↔analyticity; there wide and narrow Gaussians are both analytic = fast, so it separates the phases not at all — Fourier-dual to the width reading, equally sign-blind. Only `−φ″(0)`, at the origin, carries `Var`.)

## 4. What is exact and worth banking
**(i) Dictionary.** `S=1/Var(a)=−1/φ″(0)`, `G″=S−D`; the target is `−φ″(0)≤1/λ` (mod `D`).
**(ii) Bandwidth = Fisher.** Plancherel `∫t²|φ(t)|²dt=2π∫(p′)²dx` and `∫|(√p)′|²=J/4` show the characteristic function's frequency spread is governed by the marginal's Fisher information `J=Var_μ(E[f|a])` (E2). Hence Fourier-uniformity, E2 (information) and k0-FREEPROB (`1/V`) are the same functional in three costumes.
**(iii) Equidistribution is exact but gap-blind.** Single-link marginals are Haar, `β`-independent, *identical for free and YM* (`β=0` is literally Haar). All gap information lives in cross-link correlations = the transfer spectral gap (PROOF-CONCENTRATION) — which `J` already is.

## 5. Five-point
1. **Nonperturbative — PASS** (exact Fourier/Plancherel/Peter–Weyl/Stam, all `β`).
2. **Sign-controlling `G″≥λ` — FAIL.** Stam is a `Var`-**lower** bound; fast decay ⇒ large `Var` ⇒ gapless; Weyl/equidistribution is the gapless signal.
3. **Non-circular — FAIL.** Bandwidth `=J=Var_μ(E[f|a])`; `−φ″(0)=Var=S^{-1}`; single links equidistribute gap-blindly ⇒ the only lever is cross-link = the gap.
4. **Non-SUSY + fails free field — FAIL.** Free massless is the strongest-decay object and gapless ⇒ the tool **falsely gaps** it (smoothness reading: distinguishes nothing).
5. **Uniform `≥λ` — FAIL.** A bandwidth floor `J≥λ` yields `S≥λ−(J−S)`; the Stam slack `J−S>0` is E2's non-Gaussian deficit `D`, so one needs `J−S→0` = Gaussianity = H1.

**Sharpest line.** Fourier uniformity is variance in disguise and runs the wrong pole: the free massless IR mode is the *maximally* Fourier-decaying object, so "strong decay ⇒ gap" is refuted by the free field itself; the direction that could help — a *lower* bound on the bandwidth `J` — is, via `J=Var_μ(E[f|a])` and `G″=1/Var−D`, exactly H1/E2 in a third costume. Missing step unchanged: a uniform upper bound `Var(a_{k→0})≤1/λ`, equivalently `J≥λ` saturated by the Gaussian (`G″Var=1`), which no equidistribution estimate supplies — equidistribution controls the *spread* (large-`Var`) end.

*Refs.* INVENT-VARIANCE-INEQ (`G″=S−D`, `S=1/Var`); E2-INFORMATION (`J=Var(E[f|a])`, `G″Var=1` iff Gaussian); k0-FREEPROB-AVERAGE (`1/V`, wrong normalization); k0-TRANSLATION-RIGIDITY (named-form-kernel no-go); C8-SPECIES-NOGO (sign/no-referent discipline); PROOF-COARSE-POINCARE (H1 = `Var≤1/λ`); PROOF-CONCENTRATION (correlation decay downstream); ABC-NEW-INEQ (Ward `E[af]=1`). External: Weyl 1916 / Weyl 1927 (equidistribution criterion); Davenport, *Analytic Methods for Diophantine Equations* (Weyl sums); Donoho–Stark, SIAM J. Appl. Math. **49** (1989) 906 (uncertainty); Lieb, *Inequalities* (Hausdorff–Young); Peter–Weyl 1927; Stam, *Inf. Control* **2** (1959) 101; Guo–Shamai–Verdú, IEEE IT **51** (2005) 1261 (I-MMSE/Fisher).
