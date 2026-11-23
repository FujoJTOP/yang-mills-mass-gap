# E2-INFORMATION — species E2: I-MMSE / Fisher / rate–distortion for the coarse marginal. A real exact identity and a real discriminator — but every information ceiling that bounds `Var(a)` above is the gap wearing an entropy.

> **判词 (verdict).** The I-MMSE program yields **one genuine new exact weapon and one genuine no-go.** **(Weapon)** The Ward–Stein force is *literally the score* of the coarse marginal (`∂_a log p=−m`, `m=E[f|a]`), so the coarse Fisher information is the variance of the conditional Ward force, `J:=E_p[(∂_a log p)²]=Var(E[f|a])`, and the Guo–Shamai–Verdú identity `I(snr)=½∫_0^{snr}mmse` gives `Var(a)=mmse(0)`. The three "deficits" — Fisher `S−J`, dynamical `S−G″=D`, MMSE `Var−mmse_G` — are **one object**: `G″Var(a)=1 ⟺ J=S=G″ ⟺ p Gaussian`. **(No-go)** The literal slogan *"capacity upper-bounds MMSE=variance"* holds only at `snr>0`, where it bounds the **conditional** `E[Var(a|Y_s)]`, **not** the prior `Var(a)=mmse(0)`; capacity bites at `snr>0`, variance is `snr=0`, so it gives nothing at the target. Rate–distortion duality is the correct direction and it bounds distortion **from below** ⇒ `Var(a)≥`, the **wrong sign**. And Gaussian extremality (`I(snr)≤½log(1+snr Var)`, `mmse(s)≤Var/(1+sVar)`, equality iff Gaussian) makes the free field the entropy-MMSE **maximizer**, so no ceiling separates gapped from gapless without inserting `Var` — circular. **Passes #1, #4; fails #2, #3, #5.** The acid test is met *in letter* (RHS diverges free / finite YM) and *fails in spirit* (the divergence is `Var(a)` itself).

## 1. The forge — the coarse-marginal I-MMSE / Fisher identity (exact, all β)
Chart `X=log U`, `a=⟨v,X⟩`, `f=∂_n V=⟨v,∇V⟩`, `μ∝e^{−V}`. Ward (`ABC-NEW-INEQ`): `E[af]=1`. With `p` the marginal of `a`, `m(a):=E[f|a]`, `∂_a log p=−m`, so the score is the conditional Ward force and
```
   J := E_p[(∂_a log p)²] = Var_p(m) = Var_p(E[f|a]).
```
GSV identity (exact, nonperturbative; `Y_s=√s·a+Z`): `I(s) = ½∫_0^{s}mmse(t)dt`, `mmse(t)=E Var(a|Y_t)`, `mmse(0)=Var(a)=1/S`. de Bruijn: `d/dt h(a+√tZ)=½J_t`. **Unification (new):**
```
   p(a) Gaussian  ⟺  J = S = G″  ⟺  G″·Var(a) = 1  ⟺  D = 0
   (Fisher deficit = dynamical deficit = MMSE deficit = the non-Gaussianity of the k→0 marginal).
```
So `J`, `D`, and the MMSE non-Gaussianity deficit `½log(1+sVar)−I(s)` are three measurable readings of one wall; this *is* worth banking.

## 2. Where "capacity upper-bounds variance" breaks
`I(snr)≤C` ⟹ `mmse(snr)≤2C/snr` **only for `snr>0`** — a bound on the *residual* `E[Var(a|Y_s)]`, not the *prior* `Var(a)=mmse(0)`, which sits at `snr=0` where `2C/snr→∞`: vacuous exactly at the IR mode. At small snr `I(snr)=½ snr Var(a)+O(snr²)`, so `2C/snr≈Var(a)` — tautology. Rate–distortion duality (Gaussian is worst-case source, `R(D)≤½log⁺(Var/D)`) gives `D≥Var·2^{−2R_max}`, a **lower** bound. Every correct information inequality bounds `Var(a)` **from below**; the target needs **above**.

## 3. Acid test — free field is information-extremal
Free massless `k→0`: `p` is Gaussian, `Var(a)=c/k²→∞`; at fixed noise `I(snr)=½log(1+snr·c/k²)→∞`, `mmse(snr)→c/k²`. So any ceiling `I≤C<∞` **fails for the free field** and **holds for gapped YM** (`Var→1/m²`, `I→½log(1+snr/m²)`): the bound "diverges free, finite YM" — the **letter** of the acid test is met. But the ceiling that does it is `C=C(Var(a))`: its finiteness **is** `Var(a)<∞` **is** `S≥λ` **is** H1. The discriminator is a bijection of the thing it purports to bound; the divergence is the free covariance, not an independent channel fact. `β=2.5` YM is genuinely non-Gaussian (`D/S≈2%`), so `J≠S` there — but that deficit is bounded *by* `Var`, sign again wrong.

## 4. Exact missing step
The only non-circular hope is a ceiling on the **IR capacity from short-distance physics** via strong data processing across scales: `I(a_{k=0};UV) ≤ η·I(a_{k=0};all)`. But mass generation *is* UV–IR entanglement: gapped YM has *more* coupling (`I>0`), free field has `I=0` (factorized Gaussian). An **upper** bound on `I(a;UV)` pushes toward decoupling = **toward the free field** = gapless; the SDPI contraction `η` for a Gaussian transfer operator *is* the spectral gap. So the ceiling is small ⟺ the mode is rigid ⟺ the gap: `η=1−cλ/m²`. Both endpoints collapse to H1.

## 5. Five-point verdict
1. **Nonperturbative — PASS.** Exact Ward score + GSV + de Bruijn, all β.
2. **Sign-controlling `m_eff²(k→0)≥λ` — FAIL.** Information bounds `Var` from *below* (rate–distortion duality; capacity bites only at snr>0). Wrong sign.
3. **Non-circular — FAIL.** Every ceiling that separates the phases is `C(Var(a))`; IR-capacity finiteness ≡ `S≥λ` ≡ H1; SDPI `η` = spectral gap.
4. **Non-SUSY + fails free field — PASS (in letter, not spirit).** Free `p` Gaussian, `Var→∞`, `I→∞`; bound diverges as required — but because of `Var`, not an independent channel.
5. **Uniform `≥λ` — FAIL.** No `λ`: no information quantity produces a positive floor without already containing `Var(a)`.

**Sharpest line.** E2 forges a clean exact identity — the coarse score is the Ward force, so `J=Var(E[f|a])`, and `G″Var(a)=1`, `J=S`, `D=0` are one non-Gaussianity — but the species is **sign-wrong for the target**: I-MMSE, rate–distortion and SDPI are mathematically incapable of an *upper* bound on a variance that does not already equal it, because the Gaussian (free) field is the entropy/MMSE extremizer and "gapped" = "IR mode information-poor" = "`Var` small" — the same statement twice. The exact missing step is unchanged: an upper bound on the coarse force variance `Var(G|0)` at `k→0` that is *not* derived from `Var(a)` — i.e. H1, written this time in the language of entropies.

*Refs.* `ABC-NEW-INEQ` (Ward `E[af]=1`, `G″Var(a)=1` iff Gaussian); `C7-SPECIES-CONSTRUCT` (variance-upper-bound dead by identity shape; `E[aG]=1−E[aΦ′]`); `CONTINUUM-VARIANCE-DEFECT` (`G″=S−D`, `D/S≈2%`); `OTHER-H1-POINCARE` (`Hess V_eff=G″`, blocking keeps k→0, `PᵀHP|ranP` −0.3); `INVENT-VARIANCE-INEQ` (SCBL `G″=S−D`, `S=1/Var`). External: Guo–Shamai–Verdú, IEEE Trans. IT 51 (2005) 1261 (I-MMSE); Cover–Thomas, *Elements of Information Theory* (rate–distortion duality, Gaussian worst-case); Stam 1959 / de Bruijn (Fisher–entropy); Ledoux, *Concentration of Measure* (SDPI, entropy methods); Gross 1975 (log-Sobolev contraction).
