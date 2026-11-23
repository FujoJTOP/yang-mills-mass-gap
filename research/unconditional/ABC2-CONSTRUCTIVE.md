# ABC2-CONSTRUCTIVE — the coercivity-weighted jet is ONE coherent functional and it completes the region; on the typical-field high-momentum shell its multiplier is O(1), its value is the *conditioning defect of equipartition*, and that defect is (★)=H1

> **判词 (verdict).** The trapping region **completes**: one coercivity-weighted Hessian-jet functional partitions **all** of config-space, and the open cell is exactly the physical typical-field/high-momentum cell where the gap lives. But completion of the *region* is not closure of the *sign*: on that cell all four handles are simultaneously at their null — the coercivity multiplier `e^{Q/4g²}=O(1)`, the WSTC exponent `→1`, AF has `g≈1.5≫0`, and `N/12` is ~10× short of `2g|F|`. There the functional reduces to the single exact identity `G″(0)=E[V″|0]−Var(V′|0)` (`V=S_W+Φ`), whose unconditioned value is *exactly* the equipartition `E[V″]=Var(V′)`, so `G″(0)` **is** the conditioning defect of equipartition. The missing bound is that defect from below = `Var_μ(a_{k→0})≤1/λ` = **H1**. Not a completion of the sign; a completion of the *region*, with the exact failing estimate now localized to one cell and written in closed form.

**§1 One functional, exact.** `μ∝e^{−V}dX`, `V=S_W+Φ`, `Φ=−log J`. Jet `‖R‖_X=sup_{|v|=1}sup_A e^{Q/4g²}|δ²_{vv}R|ℓ`, `Q=∫TrF²` (PROOF-TRAPPING §1). On a mode `a=⟨v,X⟩`, IBP (periodic, no boundary) gives the *full-potential* identity
`G″(0)=E_μ[V″|0]−Var_μ(V′|0)=Δ_dyn+W_geo+X`, `Δ_dyn=E[H|0]−Var(G|0)`, `W_geo≥N/12`, `X=−2Cov(G,Φ′|0)`.
Unconditioned, `E[V″]=E[(V′)²]=Var(V′)+(EV′)²`; at the `P`-invariant vacuum `EV′=0`, hence `E[V″]=Var(V′)` **exactly**. So `G″(0)=[E(V″|0)−E V″]−[Var(V′|0)−Var V′]`: the gap *is* the conditioning defect of equipartition. This is the coherent functional — all four pieces are bounds on the two halves of one difference.

**§2 Three cells genuinely signed (the multiplier >1 there).** Shallow `2g|F|≤N/12`: `Φ″⪰N/12 ⇒ K^V⪰−2g|F|+N/12≥0`, BL fires, `D≤0`, `G″≥S`. High field: WSTC bounds the cone's *conditioned variance mass* `≤(cβκ_max)²e^{−cβ(κ²−κ_typ²)N_□}`, uniformly in `L`; the high-field escape — including the high-field version of `k→0` — is excluded. AF-UV tower `g_j→0`: weight `e^{−cμ_j⁴/g_j⁴}→0`. Each uses `e^{Q/4g²}>1`, `Q≥Vκ²`; each is a lower bound.

**§3 The open cell — the multiplier's null.** Typical field, high momentum: `Q≈Vκ_typ`, so `e^{Q/4g²}=O(1)` — coercivity is null. And `2g|F|≳k²` with `k` reaching the UV block (`k_*≳k_max`, K2), so *no shell is cone-free*: bare block Hessian `k_b²−2g|F_b|+N/12` is indefinite (W6-BE: 21–42% negative directions, `−1.12β`, 90–97% in the fluctuation block). WSTC's exponent `e^{−cβ(κ²−κ_typ²)}→1` at `κ≈κ_typ`; AF has `g≈1.5`; `N/12≈0.17` is 10× short of `2g|F|`. Equipartition holds *within* the cell, `Δ_dyn≈0` (measured `E[H]−Var≈0`), so `G″≈W_geo+X` — and this cell contains the `k→0` goal mode (`PᵀHP` min≈−0.3).

**§4 Why the four do not compose to a sign.** On the open cell `W_geo≥N/12` is real (Prékopa, every mode) but `X` is sign-uncontrolled. The only universal bound is Cauchy–Schwarz `|X|≤2√(Var G·Var Φ′)`; the measured `Var G≈2`, `Var Φ′≈E[Φ″]−W_geo≈0.02` give a bound `≈0.40` against slack `W_geo≈0.32`, so `G″≥0.32−0.40=−0.08<0` — short by `O(0.1)=λ`, *one variance bound short*, not a magnitude gap. (Measured `W_full≈0.54>W_geo≈0.32` ⇒ `X≈+0.22`: the sign is carried by `Cov(G,Φ′|0)<0`, the anti-localization of the Jacobian force — peaked at `F≈0`, the log-cut — against the Wilson force `∝|F|`. That is a cross-covariance rearrangement inequality; it is not rigorous, and SU(N) has no lattice FKG — which, being a *lower* bound on `Cov`, would point the wrong way regardless.) So the cell's sign is the conditioned joint statement `Δ_dyn+X≥λ−W_geo`, equivalently — via the exact Ward `G″·Var(a)=1` — `Var_μ(a_{k→0})≤1/λ`.

**§5 Five-point verdict.**
1. **Nonperturbative — PASS.** Exact identity, Prékopa/BL, WSTC mass, Savvidy spectrum; no loops.
2. **Sign-controlling = LOWER bound — FAIL at the cell.** Shallow + high-field + AF-UV are genuine lower bounds; the typical/high-momentum cell and the `k→0` mode are indefinite.
3. **Non-circular — PASS on the bulk, FAIL at closing.** The open-cell requirement is `Var(a)≤1/λ` = H1, zero slack.
4. **Non-SUSY + fails free field — PASS.** Free field `F≡0 ⇒ Q=0 ⇒` multiplier `1`, cone empty, `Φ″=0 ⇒ G″=0`, `f₁²=0`; the acid test is met (no spurious positivity).
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Exactly H1.

**Exact failing estimate (★₂).**
```
G″(0)|_{open} = E_μ[V″|0] − Var_μ(V′|0) ≥ λ,   V=S_W+Φ,
⟺ [(E[V″|0]−E V″) − (Var(V′|0)−Var V′)] ≥ λ      (conditioning defect of equipartition)
⟺ Var_μ(a_{k→0}) ≤ 1/λ = H1.
One universal bound: |X| ≤ 2√(Var G·Var Φ′) ≈ 0.40 > W_geo ≈ 0.32 ⇒ G″ ≳ −0.08.
```

**Sharpest line.** The coercivity jet is the right single functional: it signs shallow (`N/12`), high field (WSTC mass), and the AF-UV tower (`g→0`) with one multiplier `e^{Q/4g²}`. But the open cell is precisely where `Q` is *typical*, so the multiplier returns to `1`; there the functional is the exact identity `G″=E[V″|0]−Var(V′|0)`, whose unconditioned value is equipartition, and the sign is the conditioning defect — a cross-covariance no inequality certifies. The region completes; the sign does not. `(★₂)=H1`, now one cell and one difference.

*Refs.* ABC2-SYNTHESIS (plane map, one open cell); ABC-CONSTRUCTIVE §1–2; ABC-NEW-INEQ (WSTC `≤e^{−cβ(κ²−κ_typ²)N_□}`, `G″Var(a)=1`); PROOF-TRAPPING-REGION §1–4; PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥N/12`); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`); FORGE5-SL-UV-CONE-EMPTY (K2 `k_*≳k_max`); FORGE5-CONE-VARIANCE (`E[H]=E[f²]` unconditioned; conditioning = the gap); OTHER-SAVVIDY-SIGN (Var_Φ=O(1); log-cut anti-localization); SUB1-COV-BOUND (`D≤cW_geo` false, ~10×); OTHER-H1-POINCARE (`Hess V_eff=G″`); W6-BE-at-2.5 (21–42% neg, `−1.12β`, 90–97% fiber); W7-Gpp0 (`G″`+0.51, `W_full≈0.54`); Prékopa 1971; Brascamp–Lieb 1976; Holley–Stroock 1987; Savvidy PLB 71 (1977) 133.
