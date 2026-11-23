# LOC-LOCALIZE-PROVE — the "fixed-cell oscillation" bound is not a free-standing route: on a *fixed* cell the massless free field already satisfies it (indeed `H_0→0`), so the capacity step `H_0=O(1) ⇒ Φ_eq≥ce^{−H_0}` is false without localization; the cell trichotomy is exhaustive — fixed ⇒ vacuous, natural (`λ^{−1/2}`) ⇒ `κ/λ` circular, correlation (`ξ`) ⇒ `ξ` fixed = the gap

> 2026-11-23 · Proof attempt on the C-BBS-ARRHENIUS closing reading. Target: `H_0 := sup_{|a|≤w}|V̄(a)−V̄(0)| ≤ O(1)` uniformly in `(L,a)`, at fixed `w~ξ`, claimed "strictly weaker" than `Var(a)≤1/λ`/`Hess⪰λ`, closing the λ-bootstrap to `λ_0=¼e^{−2H_0}>0` with no seed.
>
> **判词 (verdict): NOT 5/5 — 条件定理 / 精确缺失的一步.** The bound is *provable and vacuous*: Prékopa (`V̄″≥N/12`) is a **lower** bound on curvature (an upper bound on spread), never an upper bound on `V̄`-rise; `WSTC` starves the cone but the bottleneck sits at typical action. The step that turns `H_0=O(1)` into a conductance floor is where localization hides, and the **massless free field refutes that step while passing the bound** (`H_0→0`, `Φ_eq~1/L²→0`). "Strictly weaker" is wrong: on a fixed cell the statement is satisfied by the gapless theory; on the mode's natural cell it *is* `κ/λ` circular; the only other cell is `ξ=1/m`, knowing which is `m>0`. The exact missing step is **mass concentration** (`ν({|a|≤w})≥c₀`), equivalently `Var(a)≤O(1)` = H1 — not supplied by any of the three tools.

## 1 · What the local tools actually give (and their directions)
`V̄(a)=−log p(a)`, `V̄″=G″=E[H|a]−Var(f|a)` (ABC-NEW-INEQ §1), `V̄′(0)=0` since 0 is the min (`G″(0)=+0.3` measured).
- **Prékopa / Lemma B** (`Φ″⪰N/12`): on the core `{|a|≤w_c}` where Savvidy concavity is inactive, `V̄″≥c>0` ⇒ `p(a)≤p(0)e^{−ca²/2}` ⇒ **`Var(a)|_core≤1/c`**. This is a genuine local-confinement lemma — but it bounds spread *from above on the core only*; it says nothing on the thermal shell / tail, where the mass actually sits. It is a **lower** bound on curvature, hence cannot bound `V̄`-rise (`H_0`) above.
- **WSTC** (`ν(C_κ)≤e^{−cβΔκ²N_□}`): kills the cone — correct and `L`-uniform — but the sign-carriers are at `κ_typ` (ABC-NEW-INEQ §4), exactly the region it is silent on. Wrong region for `H_0`.
- **Non-upper-boundability (Lemma A).** Even a *normalized, convex* marginal can have arbitrarily large fixed-cell oscillation: `p=(K/2)e^{−K|a|}` is convex, `∫p=1`, `Var=2/K²→0` (a *good* gap!), yet `H_0=Kw→∞`. So `V̄″≥c` (convexity) provably does **not** bound `H_0` from above, and neither necessity nor sufficiency of `H_0=O(1)` for a gap holds at the level of convexity alone.

## 2 · The vacuity — the free field passes the bound (Theorem C)
Free massless coarse mode: `V̄_free=½k̂²a²`, `k̂=2π/L`. On the **fixed** cell `|a|≤w`: `H_0=½k̂²w²→0`. So the free field satisfies `H_0≤O(1)` trivially (*better* than YM). Yet its conductance is `Φ_eq^free~1/L²→0` (mode spreads over the box). Hence the capstone step
```
   H_0 = O(1)  ⟹  Φ_eq ≥ c e^{−H_0} > 0          (FALSE, counterexample = free field)
```
fails: `Φ_eq` is the `inf` over **all** cuts including the far tail `|a|~L`, where the free mode pays `~2logL`, not the fixed-cell `H_0→0`. The local capacity argument bounds the resistance across the *fixed cell* only; `Φ_eq` is a **global** inf. **Acid-test violation, decisive:** substituting the free `H_0→0` into the bootstrap gives `λ_0=¼e^{0}=¼>0` — it would "prove a mass gap for the massless free field." A gap-closer must *fail* on the free field; this bound is *satisfied* by it. The step `(C)` — not the bound `H_0=O(1)` — is where the hidden a-priori **localization** lives.

## 3 · The cell trichotomy is exhaustive (Theorem D, a no-go on the route)
The bootstrap needs `H_0` on the cell the marginal *occupies*. Three admissible cells:
1. **Fixed** `w=O(1)` (absolute): `H_0` satisfies the free field ⇒ **vacuous** (§2).
2. **Natural / variance-scaled** `w~λ^{−1/2}`: quadratic dip gives `H_0≈½|V̄″|λ^{−1}=κ/λ`; `H_0=O(1) ⟺ λ=O(1)` = H1 ⇒ **circular** (C-BBS-ARRHENIUS §4 already shows `T(λ)/λ<1/7`, drains to 0).
3. **Correlation** `w~ξ=1/m`: asserting `ξ` is fixed *is* `m>0` = **the gap**.
No fourth cell exists; each is vacuous, circular, or the gap. So "the fixed-cell barrier closes *without* a seed" cannot be true: the only non-vacuous reading coincides with the λ-reading.

## 4 · Exact missing step (honest)
A **mass-concentration / anti-spreading** bound over a fixed cell — `ν({|a|≤w})≥c₀` uniformly in `L`, equivalently `Var(a)≤O(1)` (= H1), equivalently a uniform upper bound `sup_{|a|≤w}|V̄′|≤C` (a *force*, not curvature, bound). The natural handle is the fine cancellation `V̄″=E[H]−Var(f)=O(1)`; the naive upper bound is `E[H]=O(βN_□)=O(L^d)` — **extensive and vacuous**. Prékopa bounds curvature from below, WSTC bounds the wrong region, and Lemma A shows convexity alone cannot supply it. This is precisely H1, not a weaker statement.

**Sharpest line.** The "fixed-cell oscillation bound" is a *correct* upper bound that the gapless free field already meets (`H_0→0`) while its conductance still collapses (`1/L²`); plugging it into the bootstrap manufactures a free-field gap — the acid test fires. The capacity floor needs the mode to *be* localized in the cell, and that mass concentration is H1 itself. `H_0=O(1)` is therefore neither a weakening of H1 nor a closure: it is vacuous on a fixed cell and *is* `κ/λ` on the mode's own cell. `5/5` **not** achieved; `m_eff²≥λ_0>0` still waits on `Var(a)≤O(1)`.

*Refs.* C-BBS-ARRHENIUS (λ-bootstrap, `T(λ)/λ<1/7`, fixed-vs-`λ` reading); C-BBS-FORMULATE (Cheeger, core `N/12`, `(R2)`no-thin-neck, `(R4)`k→0=H1); POTENTIAL-THEORY-CAPACITY §2,§4 (finite-volume vacuity; free `V̄≡0⇒Cap=0`); D-BBS-MEASURE-TRAP §3 (`R4` shell=whole volume at k→0); PROOF-JACOBIAN-SAVES (`W_geo≥0` sign-only, `Φ″⪰N/12`); ABC-NEW-INEQ §1,§4 (`G″=E[H]−Var`, `G″Var=1` iff Gauss, cone vs `κ_typ`); NUMERICS-AUDIT (finite-k value < free, LB≈−0.1); G2E-RIGOROUS-NU (`m₁²≥f₁²/χ`; H1 necessary-not-sufficient); G1-DRIVE §3 (object≠gap). External: Prékopa (1971); Brascamp–Lieb (1976); Cheeger (1970); Sinclair (1992); Kramers (1940); Eyring (1935).
