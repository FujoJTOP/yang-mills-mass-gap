# SOS-CERTIFICATE — the gap as a moment–SOS relaxation: `Var(a)≤1/λ` is a clean single-entry SDP with the corpus's cleanest acid test, but its degree-2 optimum is *exactly Popoviciu* (cell diameter, measure-blind), and tightening needs transcendental moments

> **判词 (verdict).** Cast `G″(k→0)≥λ` via D3 as `Var_μ(a)≤1/λ` (Gaussian-exact 2%), `a=⟨v,X⟩` **linear** in the flat chart, so `Var(a)=M₂[v,v]` is a *single entry* of the order-2 moment matrix. Forged object: a **PSD Gram matrix** `Y⪰0` certifying `1/λ−a²=σ₀+Σg_iσ_i≥0`, a finite sign-correct bound `G″≥λ` by SDP duality; acid test cleanest in the corpus — free massless ⇒ `Var(a_{k=0})=∞` ⇒ **infeasible**, massive free ⇒ feasible at `λ=m²`. **But an exact no-go closes it:** the degree-2 relaxation of a *single linear* functional is *verbatim* **Popoviciu** `Var≤(b−a)²/4` — cell diameter, *measure-blind* — giving only `λ≲1/π²≈0.10` (4× short of `G″≈0.4`, not L-uniform); tightening needs higher moments that are transcendental (A5 8192-dim integral), size-`binom(8192+d,d)`, Savvidy-slack, `d→∞` = μ = H1. **No finite rank: bottom rung information-free, top rung the gap. #1,#4 clean; #2,#3,#5 = the wall.**

## §1 The formulation (the forged object)
`μ∝e^{−S_W−Φ}dX`, flat chart `X∈su(N)^{4L⁴}`, `a=⟨v,X⟩` linear, `⟨a⟩=0` at `k→0`. BL/D3: `G″≥λ ⟺ Var_μ(a)≤1/λ`. The order-2 Lasserre/Putinar relaxation is
```
min 0  s.t.  M₂⪰0,  L_{g_i}(M₂)⪰0  (cell [−π,π]^n ∩ |a|≤ρ, chart cuts),  M₂[v,v]=Var(a);
dual certificate:   1/λ − a² = σ₀(x) + Σ_i g_i(x) σ_i(x) ≥ 0 on K,   Gram Y⪰0 .
```
`a²` is **degree 2** — the lowest non-trivial instance; `Var(a)` is one entry, the dual one PSD Gram matrix. The sign is *right*: `Var≤1/λ ⇒ G″≥λ`, unlike B4's sign-blind Cauchy–Schwarz. (Prior `defect-bootstrap/` SDPB is a *different* object — CFT defect crossing — and confirms no working PSD/SOS layer.)

## §2 Sharp no-go: the degree-2 optimum *is* Popoviciu
A measure on `[m,M]` has `Var≤(M−m)²/4`, sharp on two endpoint atoms (Popoviciu; Bhatia–Davis). This *is* the degree-2 moment-SOS optimum: the relaxation sees only `M₂⪰0` and localizers, whose extremal solutions are two-point measures, so it returns the **cell diameter** and nothing else — `λ≤1/π²≈0.10` vs `G″≈0.4`, **L-independent and measure-blind**. But `G″>0` *is* a fact about μ's **decay** (`Var(a)≈2.5<∞`), which a measure-blind certificate cannot see. Shrinking the cell to `|a|≤1/√λ` to *make* `1/λ−a²≥0` assumes the bound (**circular**); the honest cell gives the vacuous `0.10`. **The lowest SOS instance is either the trivial diameter or the assumption of the answer.**

## §3 Acid test (cleanest)
Free massless `k=0`: `Var(a_{k=0})=∞` ⇒ for any finite `1/λ`, `1/λ−a²` changes sign where the free support reaches `|a|=1/√λ` ⇒ **infeasible**. Massive free ⇒ `Var(a)=1/m²` ⇒ Gram exists at `λ=m²`. **PASS by infeasibility — Cauchy–Schwarz is merely *satisfied* by the free field.**

## §4 Higher order cannot escape
**(i) Transcendence.** μ's flat density `e^{−S_W}J` (`J` the Haar Jacobian, not polynomial); the moments `⟨x_α⟩` pinning `M_d` are A5's 8192-dim validated integral (log `U` cuts, eigen-collisions). **(ii) Size.** `dim M_d=binom(8192+d,d)`, intractable for `d≥2` at L=4, and §2 makes `d=2` content-free. **(iii) Savvidy slack.** `Var_relax≥Var_true`, uncontrolled for indefinite `Hess S_W`; `λ` swamped unless `d` huge, `d→∞` = true moments = H1. **β-continuation:** `Hess S_W` crosses 0 at the Savvidy onset while `S(β)` falls ⇒ β=0→2.5 infeasible at the crossover; β=0 itself is *analytic* `G″=W_geo≥N/12` (Prékopa+Lemma B), not SOS-representable. Hybrid (`Φ″⪰N/12` + SOS `λ_min(Hess S_W)≥−κ`) needs `κ<N/12`: false (`λ_min=−g|B|≫N/12` ⇒ `D<S` = H1).

## Five-point scorecard
1 **Nonperturbative — PASS.** Exact SDP/LP, any β; no loops/`1/N`.
2 **Sign-controlling `G″≥λ` — SHAPE PASS / CLOSURE FAIL.** Dual PSD Gram correctly signed, unlike C–S; but degree 2 = Popoviciu (`λ≲0.1`, non-uniform); tightening = H1.
3 **Non-circular — FAIL.** Self-sized cell = circular; honest cell = vacuous diameter; `d→∞` = the measure = H1.
4 **Non-SUSY + fails free field — PASS (cleanest).** Free ⇒ `Var=∞` ⇒ infeasible; non-SUSY.
5 **Uniform `≥λ` — FAIL.** Degree-2 `λ` is cell diameter (L-independent, 4× wrong); no uniform physical `λ`.

## Exact failing / constructive estimate
```
forge:  G″≥λ ⟺ Var_μ(a)≤1/λ ;  a linear ⇒ Var(a)=M₂[v,v] (single entry)   [cleanest SDP]
dual:   1/λ−a² = σ₀+Σg_iσ_i ≥0 on K ; Gram Y⪰0   ⇒ G″≥λ                  [right sign]
deg-2:  = Popoviciu  Var≤(b−a)²/4  (endpoint atoms) ⇒ measure-blind
        group cell a∈[−π,π] ⇒ λ≤1/π²≈0.10 vs true G″≈0.4   [4× short, L-indep]
        shrink cell to |a|≤1/√λ ⇒ circular; honest cell ⇒ vacuous
higher: d≥3 needs ⟨x_α⟩ = A5 8192-dim transcendental integral; dim M_d=binom(8192+d,d)
        Var_relax≥Var_true slack (Savvidy); d→∞ = true μ = G″ = H1
β=0:    G″=W_geo≥N/12 analytic (Prékopa) — NOT SOS-representable
free:   Var(a_{k=0})=∞ ⇒ 1/λ−a² sign-changes on free support ⇒ INFEASIBLE  [acid PASS]
```

**Sharpest line.** Cleanest in the corpus: the dual is a literal **PSD Gram matrix of correct lower-bound sign**, and the acid test is a genuine **infeasibility** (a free field admits no `G″≥λ>0` certificate; Cauchy–Schwarz is merely *satisfied* by it). The cleanness is the trap — for a *single linear* functional degree 2 is **exactly Popoviciu** (cell diameter, measure-blind), while `G″>0` is a statement about μ's *concentration* (`Var(a)<∞`), invisible at degree 2; reaching it needs higher moments that are transcendental (Haar Jacobian / A5) and Savvidy-slack. Bottom rung information-free, top rung H1. **A clean SDP, the cleanest acid test, a new exact coordinate on the wall (degree-2 = Popoviciu ⇒ concentration-blind) — not the gap.**

*Refs.* **Internal:** A5-RIGOROUS-NUMERICS (8192-dim validated integral, β=0 quadrature); COMPUTE-GPP-SU3 (`G″≈0.4`, W-dominance); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`); D3-SUSCEPTIBILITY-VARIATIONAL (`G″≥λ ⟺ Var(a)≤1/λ`); B4-VARIANCE-MAG (C–S sign-blind, free-field-satisfied); PROOF-JACOBIAN-SAVES (`Φ″⪰N/12`, (◆), `W_geo≥0`); C7-SPECIES-CONSTRUCT (§e); PROOF-SAVVIDY-WGEO-DOMINATION (`λ_min=−g|B|`); PROOF-P3-UNCONDITIONAL (`β<1/48`); defect-bootstrap/SDPB-DERISK (no PSD/SOS layer; different object). **External:** Lasserre, *SIAM J. Optim.* **11** (2001) 796; Parrilo, *Math. Program.* **96** (2003) 293; Putinar, *Arch. Rational Mech. Anal.* **142** (1993) 93; Handelman, *Pacific J. Math.* **132** (1988) 35; Laurent, *Emerging Appl. Alg. Geom.*, IMA **149** (2009) 157; Popoviciu (1935); Bhatia–Davis, *Amer. Math. Monthly* **107** (2000) 353; Prékopa, *Acta Sci. Math. Szeged* **34** (1973) 335; Brascamp–Lieb, *J. Funct. Anal.* **22** (1976) 366; Savvidy, *Phys. Lett. B* **71** (1977) 133.
