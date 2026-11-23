# FORGE5-SL-IR-SINGLE — the single IR `k→0` mode: can SL close it without the backward? Mode-blind terminal slice, tautological gap equation

> **2026-11-23 · Verdict (判词).** Granting the UV half closes by cone-emptying, the entire residual is the single `k→0` variance `Var_μ(a)`. All three SL shortcuts fail, honestly: **(A)** the `k→0` mode is **not** what the terminal slice "sees most favorably" — it is the *softest*, hence *largest*-variance mode `1/(m_eff²+t)`, and at `t=t*` the slice is **mode-blind** (every eigenvalue `→1/t*`); **(B)** `1/Var_μ(a)` is **not** directly computable — the only exact identity (Schwinger–Dyson `E[H]=E[f²]`) is the *unconditioned* variance (= equipartition), missing exactly the conditioning that *is* the gap; **(C)** a self-consistent SL gap equation **exists** but is empty — `m_eff²=M(t*)−t*` is the identity `m=m` in the Gaussian/RPA closure and H1 in the exact one. The `k→0` residual is genuinely **H1**.

## A — "lowest mode, seen most favorably"? Two category errors
*Transfer-matrix.* `k→0` as a *fluctuation* is the lightest **excitation** `λ₁`, **not** the PF eigenvalue (`λ_PF=e^{−aE₀}`, the vacuum); the terminal measure is a near-Gaussian on the vacuum and does not isolate the coarse excitation.
*Covariance.* BL gives `Var_{μ_t}(a)≤1/(m_eff²+t)`; the softest mode has the *largest* right side at every `t`, so the bound is **loosest** there — the opposite of favorable. As `t→∞` all eigenvalues `→1/t`: the slice loses mode information. The IR content lives in the covariance **process** `A_t`, initial value `A_0=Var_μ(a)` = the unknown.

## B — "variance directly computable"? It is the gap
H1 ⟺ `Var_μ(a)≤1/λ` (PROOF-JACOBIAN-SAVES (◆); PROOF-COARSE-POINCARE). `G″≈1/Var(a)` is exact only for a Gaussian marginal (`G″=S−D`, `S=1/Var`, `D`=BL defect, INVENT-VARIANCE-INEQ), else measured at 2–20%. The only exact "direct" identity is `E_μ[H]=E_μ[f²] ⇒ Var_μ(f)≤E_μ[H]` — the **unconditioned** variance, which at the symmetric vacuum equals `E[H]` exactly. The gap is 100% a *conditioning* effect, `G″(0)=E[H|a=0]−Var(f|a=0)+W`, and `Var(f|a=0)` is H1's object. Finiteness of `Var_μ(a)` is `m_eff>0`: the massless slice has `Var_{μ_t}=1/t`, `∫₀ dt/t=∞` at the free zero mode (PROOF-STOCHASTIC-LOCALIZATION §3). Variance-positivity bounds moments **below**; nothing upper-bounds the edge **above**.

## C — the self-consistent IR gap equation: constructed, then shown empty
Define the SL curvature `M(t):=1/Var_{μ_t}(a)`. Gaussian-exact: `M(t)=m_eff²+t`, `dM/dt=1`. The un-shift is the **gap equation**
> **(GE)**  `m_eff² = M(t*) − t*,   M(t*)=1/Var_{μ_{t*}}(a),   t*≥1.12β≈2.8.`

BL at `t*` supplies only `M(t*)≥ε`, so **(GE)** yields `m_eff²≥ε−t*<0` — vacuous. The *actual* `M(t*)` is the exact tilted marginal; in the Gaussian closure it is `1/(m_eff²+t*)`, making **(GE)** the identity `m=m` — **no nonzero branch selected**. RPA reproduces this: every `m` solves it; the nonzero branch is the dynamical self-energy `Π`, i.e. the conditioned force variance = H1 (PROOF-COARSE-POINCARE (b)). Breaking the tautology needs the *non-Gaussian* (fourth-moment/conditioning) part of the coarse marginal — exactly the data the backward would transport. **(GE) is the backward cost as a fixed point: closing it already presumes `Var_{μ_{t*}}(a)` exactly, i.e. `m_eff²`.**

## The martingale no-go — why no re-arrangement helps
`A_t=Var_{μ_t}(a)` is a nonnegative martingale with `A_∞=0` and `A_0=Var_μ(a)` **deterministic, equal to the unknown**. Optional-stopping/barrier arguments — how SL proves KLS — bound the *hitting time* `A_t→0`, the *decay rate* (a scale-invariant isoperimetric ratio); they never bound the *initial value*: `E[A_T]=A_0` is tautological, and a bound on a future `A_T` does not bound `A_0`. H1 is the **absolute magnitude** `A_0≤1/λ`; SL controls a quantity whose initial datum *is* the answer. The `+t` tilt *translates* that magnitude, invertible only through the backward `1+t*/m_eff²` (`≈10` at `m²≈0.3`; `→10⁶` massless). KLS is scale-free; the gap is scale-covariant.

## Five-point verdict
1. **Nonperturbative — PASS.** `Hess S_t=Hess S+tI` exact; no expansion.
2. **Sign at the single mode — FAIL.** Un-tilted `G″_∞>0` ⟺ `Var_μ(a)≤1/λ` = H1; SL certifies only the tilted sign.
3. **Non-circular — FAIL.** (GE) is a tautology (Gaussian/RPA) or H1 (exact); exact `M(t*)` = the gap.
4. **Non-SUSY / free-field — PASS.** Massless: `Var=∞`, `M(0)=0`, `λ=0`; the `∫dt/t` divergence pins the free zero mode.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL.** Only `≥λe^{−t*V}` (backward) or `≥0` (trivial); no uniform `λ`.

**Sharpest line.** Each shortcut collapses on inspection: the terminal slice is *mode-blind* and the softest mode is its *worst* case, not its best; the variance *is* the gap (no identity escapes the conditioning); and `m²=M(t*)−t*` is `m=m` in the Gaussian closure or H1 in the exact one. The structural reason is one line — `A_t` is a martingale whose initial value `A_0` is the unknown, so every stopping/isoperimetric argument bounds the *decay rate*, never `A_0`. SL hands H1 a sharper face (`G″=1/Var(a)`, zero hidden geometry) and proves it *KLS-orthogonal* — not a closure. **The IR residual is genuinely H1.**

*Refs.* Eldan, ECP **18** (2013) 63; Chen, GAFA **31** (2021) 34; Klartag, Mem. AMS **249** (2017); Brascamp–Lieb, JFA **22** (1976) 366; Holley–Stroock, Ann. Probab. **15** (1987) 1541; Efron–Stein (1981). Internal: HARDFORGE-SL-FULL (Step 4 = H1); FORGE5-SL-BACKWARD (`1+t*/m_eff²`); PROOF-STOCHASTIC-LOCALIZATION (free `∫dt/t`); PROOF-JACOBIAN-SAVES (`G″≈1/Var`, (◆)); INVENT-VARIANCE-INEQ (`G″=S−D`, `S=1/Var` exact Gaussian); FORGE5-CONE-VARIANCE (`E[H]=E[f²]`, conditioning=gap); INVENT-CONVEXIFICATION (Prékopa fails); PROOF-IR-H1 (H1 = bare gap).
