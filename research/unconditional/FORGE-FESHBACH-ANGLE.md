# FORGE-FESHBACH-ANGLE — Schur complement + graph-angle on the Savvidy complement: the cross-block *adds* (C≺0), so the angle never manufactures the positivity A lacks; exact (angle 0) only on the frozen background ⇒ vacuous, O(1) at k=0 ⇒ fatal

> 2026-11-23 · forge of a Feshbach/Schur-complement tool. Split the Hessian `H=[[A,B],[Bᵀ,C]]` with `N` the known anti-self-dual LLL (`rank r_−`, depth `−g|B|`), `P=N^⊥`, `B=PᵀHN` controlled by the graph angle `θ(P,N)`. **判词: 3/5 — not genuine; it reduces *exactly* to H1.** The decisive fact is a **sign**: `C≺0 ⇒ −BC⁻¹Bᵀ=+B|C|⁻¹Bᵀ⪰0 ⇒ Σ_P=A+B|C|⁻¹Bᵀ⪰A`. The task's proposed `a₋−‖B‖²/‖C‖` is the wrong sign for a *negative* eliminated block: **the cross-block only helps the survivor** (Lemma NS; verified `λ_min(Σ−A)=−6·10⁻¹⁶`, Haynsworth `inertia(H)=inertia(C)+inertia(Σ)`). So "small angle keeps `Σ` positive" is a category error — `Σ>0 ⇔ A>0`, and `A=PᵀHP` *is* the gap; the angle can certify an already-positive A but cannot rescue its measured min `−0.3`. Exact case (`θ=0`) = frozen constant background (vacuous, `g`-dependent); physical case `ε=sin²θ=O(1)` at `k=0` (fatal). Momentum separation does **not** rescue it — the negatives bottom at `k_∥=0`, the coarse mode's own regime.

## 1. The construction and the sign
`C=NᵀHN=−g|B|·I_{r_−}` (LLL, diagonal), `A=PᵀHP`, `B=PᵀHN`. Feshbach–Schur onto `P`:
```
Σ_P = A − B C⁻¹ Bᵀ = A + B|C|⁻¹Bᵀ ⪰ A ,   since C≺0 ⇒ C⁻¹≺0.
```
So `λ_min(Σ_P) ≥ λ_min(A)`; the correction is `+‖B‖²/g|B|`-type, **additive**. Haynsworth `inertia(H)=inertia(C)+inertia(Σ_P)`: `H` has exactly `r_−` negatives iff `Σ_P≻0`, iff `A≻0`. The graph angle enters **only** through the bound on `A`. For unit `p∈P` with overlap `ε=‖P_N p‖²=sin²θ_p` onto the true negative eigenspace, writing `λ_+` for the first positive Landau rung,
```
⟨p,Hp⟩ = λ_+(1−ε) + (−g|B|)ε = g|B|(1−2ε)   [λ_+=+g|B|]
⇒  A≻0  ⟺  ε < ½  ⟺  sin²θ < ½  (θ < 45°).          (verified)
```
This is the whole content: **the angle bound is an overlap bound, and the threshold is `ε<½`.**

## 2. Exact case: `θ=0` — but vacuous
On a *frozen constant* abelian `B`, the LLL is a genuine invariant subspace of the symmetric `H` (`[H,𝒩_Landau]=0`), so `P=N^⊥` is also invariant ⇒ `B=PᵀHN=0` **exactly**, `θ=0`, and `Σ_P=A=H|_{N^⊥}` with bottom `+g|B|` (the `n=1,s=−` rung `λ=k_∥²+g|B|` at `k_∥=0`). The Feshbach complement is *exactly* positive and exactly the spectrum. **But** this is (i) background-frozen — the measure superposes/fluctuates `F`; (ii) background-dependent, `+g|B|→0` as `g→0`; (iii) the *bare* Euclidean Hessian, blind to the measure `−Var`. As a gap statement it is the same **free-field-blind** fact as Weyl-cone-emptying (D-BBS-LOWRANK-KATO): `F≡0⇒N` empty `⇒θ=0⇒Σ=A=` free Laplacian, floor `→0` at `k=0` ⇒ holds for the massless field ⇒ vacuous.

## 3. Physical case: `θ=O(1)` at `k=0` — fatal
With non-constant `F`, `B=PᵀHN≠0`, `θ≠0`. The Savvidy band bottoms at `k_∥=0`: `λ_{0,−}=k_∥²−g|B|`, **min at zero momentum** — the *same* regime as the coarse gauge-invariant mode carrying the mass. So `P` and the true negative subspace share the soft bottom: overlap `ε=O(1)`, measured `PᵀHP` min `≈−0.3<0` at β=2.5 (FORGE5 K2) ⇒ `ε≳½`, bound fails. **Momentum separation does NOT save it** — it needs the negatives at high `k`; they live at `k_∥=0` (LOWRANK-KATO (c)). `θ=0` only on the frozen background (vacuous), `θ≳45°` at the physical `k=0` (where the gap is read): NO uniform angle exists.

## Five-point verdict
1. **Nonperturbative — PASS.** Exact Feshbach block, Haynsworth inertia, no loop expansion.
2. **Sign control — PASS, favorable-only.** `Σ_P⪰A`; on the exact `+`-eigenspace `Σ⪰+g|B|` (Lemma NS, refined with `θ=0`).
3. **Non-circular — FAIL.** Needs `λ_min(A)=λ_min(PᵀHP)>0`; `A>0⇔ε<½` is exactly the `ε`-bound `= H1`.
4. **Fails free field — FAIL.** `F≡0⇒N` empty, `θ=0,Σ=A→0` at `k=0`.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** `ε=O(1)` at `k=0`; `PᵀHP` min `−0.3`.

## Exact stuck point
**`Σ_P≻0 ⇔ A≻0 ⇔ sin²θ(k=0) < ½ ⇔ ⟨coarse|P_−|coarse⟩ < ½ ⇔ H1`.** The Schur correction `+B|C|⁻¹Bᵀ` is additive and hence cannot fix a non-positive `A`; the only independent handle, the graph angle, must satisfy `ε<½` at `k=0`, where it is `O(1)` (measured `−0.3`). The angle is provably `0` on the frozen constant background (vacuous) and provably `≳½` at the physical `k=0` (fatal). There is no momentum-separation window between them — the negatives and the coarse mode occupy the same `k=0`. Lifting `A` requires the measure/one-loop back-reaction `−Var`, i.e. **H1**; the Feshbach-angle rewrite neither adds nor removes it.

**Sharpest line.** Feshbach with the correct sign is *monotone the wrong way to help*: killing a negative-definite block only *raises* the survivor (`Σ⪰A`), so the tool is **exactly as strong as `A`'s positivity — no more, no less** — and `A=PᵀHP` is the gap in disguise. The graph angle, the only genuinely new object, is benign (`0`) exactly where the answer is free/vacuous and malignant (`O(1)`) at the `k=0` bottom where the mass is printed. The route closes iff `ε<½` at `k=0` — the ε-bound of LOWRANK-KATO, i.e. H1. Not a new tool; a faithful re-derivation of the wall.

*Refs.* D-BBS-INTEGRATE-NEG (Lemma NS `Σ⪰A`, verified); D-BBS-LOWRANK-KATO (c) (`λ_{0,−}` min at `k_∥=0`, `PᵀHP` min `−0.3`, momentum-separation fails); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, ladder `λ_{n,s}`, `λ_min=−g|B|`, `r_−∝g²B²𝓥`); FORGE5-SL-UV-CONE-EMPTY (K2 β=2.5 cone, `PᵀHP` min `−0.3`); PROOF-COARSE-POINCARE (H1, `G″=E[H]−Var`); BLIND-LOWER §5–6 (lower-bound engine needs positive effective Hessian = (H)); Haynsworth (1968, inertia additivity); Kato (1966); Reed–Simon IV Thm X.12.
