# INVENT-VARIANCE-INEQ — a Savvidy-corrected Brascamp–Lieb: `G″` splits exactly into a Gaussian–Schur curvature minus a non-Gaussian excess

> **Verdict.** We construct the **Savvidy-corrected Brascamp–Lieb (SCBL)** decomposition, an *exact* identity `G″(0)=S(v)−D(v)` splitting the coarse-marginal curvature into a **Gaussian–Schur term** `S(v)` (Brascamp–Lieb with the fiber Hessian `K`) and a **non-Gaussian excess** `D(v)` (the Brascamp–Lieb defect of the conditioned fiber). It needs **no log-concavity** for the identity, and it **passes the acid test**: the massless free field returns `S=D=0`, so `G″=0` *exactly* — whereas naive BL (`Var≤E[H]`) is *satisfied with equality* by the free field and is rejected as vacuous. `D≤0 ⇔` fiber log-concave (genuine BL); `D=0 ⇔` Gaussian; `D>0 ⇔` fiber Savvidy (`K=Π₀ Hess V Π₀` indefinite) — the measured β=2.5 regime. At **β=0 it is a theorem with explicit constant `G″≥S≥N/12`**; at β=2.5 it reduces `G″>0` to the single scalar inequality `D(v)<S(v)` (measured `D/S≈2%`, PROVE-Gpp), *not* provable — the residual is exactly the indefinite fiber `H_ff`, i.e. H1 renamed. A sharpened, acid-test-clean criterion; not the gap.

## (1) The object and the split
`μ∝e^{−V}dX`, `V=S_W+Φ`, `Φ=−log J`, `Φ″⪰c_Φ`, `c_Φ=N/12` (SU2 1/6, SU3 1/4). Fix the coarse mode `v` (`|v|=1`), `a=⟨v,X⟩`; condition on `a=0`. Put `h_n=vᵀHess V v=H_vv+Φ″_vv`, `f=∂_nV=⟨v,∇V⟩=G_v+Φ′`, `Π₀=I−vvᵀ`, `c=Π₀ Hess V v`, `K=Π₀ Hess V Π₀`. Differentiating the marginal twice gives exactly `G″_0=E_μ[h_n|0]−Var_μ(f|0)` (W7). Add/subtract the Schur form `⟨K^{−1}c,c⟩`: the two terms cancel, leaving

**`G″_0 = S(v) − D(v)`**,  `S:=E_μ[h_n−⟨K^{−1}c,c⟩|0]`,  `D:=Var_μ(f|0)−E_μ[⟨K^{−1}c,c⟩|0]`.

`S` is the **Gaussian–Schur curvature**: for Gaussian `V=½⟨X,AX⟩`, `S=1/((A^{−1})_{nn})=1/Var_μ(a)` exactly (Schur complement) — so `S` *is* the "essentially Gaussian" curvature the numerics find. `D` is the Brascamp–Lieb defect of the conditioned measure `μ₀=μ(·|a=0)` applied to `f=∂_nV`.

## (2) The constructed inequality (boxed)
```
SCBL (Savvidy-corrected Brascamp–Lieb):
   Var_μ(∂_nV | a=0) = E_μ[⟨K^{−1}c, c⟩ | a=0] + D(v),
   D(v) = Var_μ(∂_nV|0) − E_μ[⟨K^{−1}c,c⟩|0],
   |D(v)| ≤ Var(∂_nV) + ‖K^{−1}‖ ‖c‖²_∞ < ∞   (∂_nV periodic, group compact, K invertible);
   D ≤ 0 ⇔ K ⪰ 0 (BL on μ₀);  D = 0 ⇔ μ₀ Gaussian;  D > 0 ⇔ K indefinite (Savvidy);
   G″(0) = S(v) − D(v),  hence  G″(0) > 0  ⟺  D(v) < S(v).
```
"Corrected" because ordinary BL is exactly `Var_μ(f|0) ≤ E_μ[⟨K^{−1}c,c⟩|0]`, valid iff the *fiber* is convex (`K⪰0`); indefinite `K` injects the excess `D>0`.

## (3) Sign of `D` — Brascamp–Lieb defect lemma
On `{a=0}`, `μ₀` has Hessian `K`. If `K⪰0`, `μ₀` is log-concave and BL (1976) gives `Var_{μ₀}(∂_nV) ≤ E_{μ₀}⟨K^{−1}∇₀∂_nV,∇₀∂_nV⟩ = E_{μ₀}⟨K^{−1}c,c⟩`, i.e. `D≤0`, equality iff `μ₀` Gaussian (Schur). If `K` is indefinite, BL reverses and `D>0` is the variance injected by the negative fiber directions. The Jacobian's own fiber is provably positive: `K^Φ⪰c_Φ` ⟹ `S^Φ=1/((Φ″)^{−1})_{nn} ≥ c_Φ=N/12`.

## (4) Acid test — mandatory
**Free massless field.** `V=½⟨X,AX⟩` Gaussian ⟹ `D=0`; the k→0 mode has `(A^{−1})_{nn}=∞` ⟹ `Var(a_{k=0})=∞` ⟹ `S=1/Var(a)=0`. Hence `G″(0)=0−0=0`. The inequality `G″≥S` reads `0≥0`: **no strict positivity is produced**, so the free-field bar is met (`f₁²=0 ⇔ S=0`). Massive free field: `Var(a)<∞ ⇒ S>0, D=0 ⇒ G″>0`. **Pure Haar (β=0):** `V=Φ` convex ⟹ `K⪰c_Φ>0` ⟹ `D≤0`, and **`G″≥S≥N/12>0`** — a theorem with explicit constant.

## (5) Why it does not close β=2.5
At β=2.5 the fiber Hessian is indefinite — 90–97% of negative BE directions lie in the fiber (W6-BE) — so `K` is not `⪰0`, `D>0`, and `G″>0` requires `D(v)<S(v)`. Numerically `D/S≈2%` (the measured Gaussian-exactness `G″≈1/Var`), so it holds; but a uniform proof of `D<S` is control of the fiber `H_ff` (a fiber Poincaré / `Var(a)<∞` = H1). The determinant-zero structure makes this literal: `det K=0` ⟺ a transverse massless mode ⟺ `S=0`. The tool *localizes* the obstruction to one measurable scalar `D(v)` and one explicit competition `D<S` — sharper than "H1" — but does not evade the indefinite fiber.

**Sharpest line.** The only log-concavity the tool needs is transverse (the fiber `K`), and that is exactly what Savvidy breaks; the SCBL identity `G″=S−D` converts the entire sign question into `D(v)<S(v)`, provable with constant `N/12` at β=0, measured true at β=2.5, unprovable there because `D>0` *is* the indefinite fiber. The acid test is passed — free massless gives `G″=0`, not `>0` — so the tool is non-vacuous; it is a sharpened criterion, not the gap.

*Refs.* PROOF-JACOBIAN-SAVES (`G″` identity, `Φ″≥N/12`); PROVE-Gpp-positive (`G″≈1/Var` to 2%, marginal Gaussian, `E[H]−Var≈0`); PROOF-COARSE-POINCARNE (fiber `H_ff` indefinite); PROOF-BORN-OPPENHEIMER (Schur undefined, `H_ff` Savvidy); REPAIR-2SCALE (fiber absorption); INVENT-CONVEXIFICATION (no log-concave surrogate; periodicity no-go); W6-BE-at-2.5 (90–97% fiber negatives); W7-Gpp0 (`G″` +0.27/+0.51); Brascamp–Lieb 1976; Prékopa 1971; Schur complement; Bakry–Émery/HWI.
