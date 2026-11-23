# B3-TYPICAL-COERCIVITY — B2's jet has no typical-field extension: every field-adapted weight rescales *magnitudes*, the sign lives in the *correlation coefficient*, and Cauchy–Schwarz misses it by 1.8×

> **判词 (verdict).** I construct the two natural typical-field weights — log-cut/center-adapted and equipartition/fluctuation-Gaussian — on B2's jet norm `‖R‖_X=sup e^{Q/4g²}|δ²R|ℓ`. **Both collapse to one single-variable bound**, for a structural reason: a coercivity weight `w(A)=e^{φ(A)}` is a *pointwise value multiplier* bounding second-variation **magnitudes** `Var(G)`, `Var(Φ′|0)`. But the open-cell sign is the **correlation coefficient** `ρ=Cov(G,Φ′|0)/√(Var G·Var Φ′)`, and common rescaling of both variances leaves `ρ` — hence `G″` — invariant. At β=2.5 (`Var G≈2`, `Var Φ′≈0.02`, true `ρ≈−0.55`) the sharpest bound any weight yields is Cauchy–Schwarz `X≥−2√(Var G·Var Φ′)=−0.40` ⇒ `G″≥W_geo−0.40=−0.08<0`. **The null multiplier is the gap** — value-coercion cannot bound a correlation. 5-pt: **1 PASS, 2 FAIL, 3 FAIL, 4 PASS, 5 FAIL.**

## §1 Two candidate weights (constructed, then broken)
**W-A · log-cut / center.** `w_A=exp[c(‖X‖²_{geo(U,−1)}−⟨·⟩)]`, up-weighting `U` near the center `−1` where the Jacobian force `Φ′` peaks. Intent: `G∝|F|` is anti-localized there, so W-A targets the region carrying `Cov(G,Φ′|0)<0`. **Break:** under `w_A`, `Var_G` and `Var_Φ′` get O(1) factors, but C-S is *scale-covariant* — `√(Var_A G·Var_A Φ′)` rescales exactly as the variances, `ρ` untouched. Bound stays `−0.40`; a pointwise weight cannot *sign* `G·Φ′` without the correlation inequality it is meant to supply.
**W-B · equipartition / high-k Gaussian.** `w_B=exp[c Σ_{k>k_*}(k²|A_k|²/g²−1)]`, adapting to the fluctuation block (90–97% of negatives, W6-BE) — B2's coercivity *localized* to high momentum. It bounds the **diagonal** single-mode `Var(f|0)≲g²/k²`, giving `Δ_dyn=E[H|0]−Var(G|0)>0` on-shell. **Break:** it shrinks `Var(G|0)` and `Var(Φ′|0)` together (same block) ⇒ `ρ` invariant ⇒ `−0.40` again; and at the `k→0` target `k_*≳k_max>k` excludes the goal mode — null precisely there.

## §2 Why the whole class fails
`G″(0)=E[V″|0]−Var(V′|0)` is an exact **difference of a first and a second moment**; a pointwise weight multiplies the integrand but the difference survives, and its sharpest completion is C-S on the two halves. Coercivity controls each half's *size*; the open cell's sign is the *angle* between them.

## §3 Exact failing estimate
```
G″(0)|_open = W_geo + X,   X=−2Cov_μ(G,Φ′|0) ≥ −2√(Var G·Var Φ′|0)      [C-S: best any weight gives]
β=2.5:  Var G≈2,  Var Φ′≈E[Φ″]−W_geo≈0.02,  W_geo≈0.32
⇒ X ≥ −2√(2·0.02) = −0.40  ⇒  G″ ≥ 0.32−0.40 = −0.08 < 0.
True (W7):  W_full≈+0.54 ⇒ X≈+0.22 ⇒ ρ = −0.11/0.20 ≈ −0.55.
Closure needs X ≥ λ−W_geo ≈ 0.1−0.32 = −0.22
   ⟺ Cov(G,Φ′|0) ≤ 0.11  ⟺ ρ ≤ 0.55;
C-S gives ρ≤1 (Cov≤0.20): short by 0.20/0.11 ≈ 1.8 in the covariance — exactly the missing
anti-correlation. Sufficient: negative association Cov≤0 (anti-FKG) — a *correlation* inequality,
which SU(N) (no FKG) neither gives nor complements.  ≡ H1 via exact Ward G″·Var(a)=1.
```

## §4 Five-point verdict
1. **Nonperturbative — PASS.** Exact identity + Ward + Prékopa; no loops. The weights themselves are nonperturbative multipliers.
2. **Sign-controlling `m_eff²(k→0)≥λ` — FAIL.** Best bound `G″≥−0.08<0` at the cell; both weights return `−0.08`.
3. **Non-circular — FAIL at closing.** The needed `Cov(G,Φ′|0)≤0.11` is `Var(a)≤1/λ` verbatim (Ward), zero slack.
4. **Non-SUSY + fails free field — PASS.** Pure YM. Free field `g→0,F≡0 ⇒ G,Φ′→0 ⇒ G″=0`, `f₁²=0`; no spurious positivity (the log-cut weight is benign because the *force* vanishes, not the weight).
5. **Uniform `≥λ` — FAIL.** Volume-uniform *bounds* exist (C-S is), but they are sign-wrong at the typical-field cell — the one cell that decides.

## §5 What would actually close it (the species that must be forged)
Not a coercivity weight but a **correlation-sign / rearrangement inequality** on the typical cell: prove the log-cut force `Φ′` is *negatively associated* with the Wilson force `G∝|F|`. That is a 6th-species object (an NQD/anti-FKG statement for lattice SU(N) invariant measure), value-coercion cannot manufacture it, and it is the only way to beat C-S on `ρ`.

**Sharpest line.** B2's multiplier is `e^{Q/4g²}`; on the typical cell `Q≈Vκ_typ` makes it `O(1)` — not an unlucky normalization but the *equipartition of the physical measure*: no coercive penalty to spend, because value and variance balance exactly. The sign is the deviation from that balance, a correlation coefficient the weight is blind to. The null multiplier *is* the gap: `Cov(G,Φ′|0)≤0.11`, which C-S misses by 1.8×.

*Refs.* ABC2-CONSTRUCTIVE (§4 C-S `−0.40` v `W_geo 0.32`, `(★)=H1`); ABC2-SYNTHESIS (open cell = H1); ABC-NEW-INEQ (Ward `G″Var(a)=1`, cone red herring); PROOF-TRAPPING-REGION (§3 value-coercion signs no variance); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, λ_min=−g|B|); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥N/12`, free β=0 closure 0.248); SUB1-COV-BOUND (conditioning = 100% of gap); W7-Gpp0 (`G″+0.51`, `W_full≈0.54`); W6-BE-at-2.5 (21–42% neg, 90–97% UV block); Prékopa 1971; Brascamp–Lieb 1976.
