# IN-9 (`SPEC-RES`) — lattice gap → continuum gap

**Verdict.** The transport *implication* is true and proven; the hypothesis `SPEC-RES` is **OPEN** and is **not** supplied by the construction (`IN-0`/`IN-1`). Its content is "the physical mass gap is lower-semicontinuous under `a→0`" — a uniform-in-`a` renormalization statement of G1 difficulty, distinct from the gap itself (G2).

## 1. Precise statement
Fixed `L`; RP transfer matrices `T_a=e^{-aH_a}`, `0≤T_a≤1`, `T_a1=1`. On `1^⊥` set `H_a:=-(1/a)log T_a` (spectral calculus), `Δ_a:=inf(σ(H_a)\{0})`, `Δ_cont:=inf(σ(H)\{0})`. `SPEC-RES`: `H_a→H` in **strong resolvent sense** — `(H_a-z)^{-1}ψ→(H-z)^{-1}ψ` for all `ψ`, `Im z≠0` (⟺ `e^{-tH_a}→e^{-tH}` strongly, ∀`t>0`).

## 2. TRUE status
- `SPEC-RES ⟹ Δ_cont ≥ liminf Δ_a` — **proved**, 4 lines (Reed–Simon I, Thm VIII.24(b): `σ(H)⊆limsup_a σ(H_a)`). Lean: `gap_ge_liminf` (`research/glue-b-blocking/lean/GL2bPrime.lean`; REPORT §3.8.1).
- `SPEC-RES` — **OPEN**, not in the 16 inputs, and **not** derivable from `IN-1`/GL3. Counterexample `T_a=diag(1,e^{-c_1a},e^{-c_2a})`: `T_a→I` (every Schwinger function converges trivially) yet `Δ_a=c_1` arbitrary. The limit encodes no `Δ`; `Δ` lives in the `O(a)` term of `log T_a`, one time-derivative beyond `IN-1`.

## 3. Attempted proof — two routes, one wall
**(a) Resolvent.** Literal `SPEC-RES` on *bare* `H_a` is ill-posed (`ℋ_a≇ℋ`). Bare cylinder vectors (`u=e^{iaA}→1`) give `⟨f,H_af⟩→0`; physical ones require field-strength `Z_a`. So `SPEC-RES` presupposes a renormalized identification `ι_a:ℋ_a→ℋ` and convergence of the `O(a)`-term — not implied by OS continuity.

**(b) Correlator (no operators).** Lattice spectral measures `ν_a`, supported in `[δΛ,∞)` by `Δ_a≥δΛ`, converging weakly ⟹ `inf supp ν ≥ δΛ` (closed-support limit) ⟹ `Δ_cont≥δΛ`. But `C_a(0)=∫dν_a→∞` (UV), so `ν_a` need not converge; one needs IR tightness (no mass migration, finite `Z_a`). Same statement.

RP thus supplies only the *setting* (self-adjoint `0≤H_a`, KL/spectral representation needed to define `Δ`) — neither the convergence nor the gap.

## 4. Exact gap
`IN-9 = (i) 4-line implication [proved] + (ii) SPEC-RES [open]`. (ii) = **uniform-in-`a` renormalization of the transfer-matrix generator**: existence of `Z_a` and comparison maps `ι_a`, with `H_a`-matrix elements converging in physical normalization — a G1 obligation, isolated from the G2 half (`Δ_a≥δΛ`, input `IN-3/4`/`hcert`).

## 5. Flip condition
Flips if (A) the construction yields a uniform field-strength/tightness bound (`sup_a Z_a^{-1}C_a(t_*)<∞` at one physical `t_*>0`) **and** `H_a`-matrix-element convergence; or (B) weak convergence of `{ν_a}` is shown. Both are G1 results, open today. Note the conclusion is `Δ_cont≥liminf Δ_a` (`≥`, not `=`); equality would need norm/Hausdorff convergence (harder, unnecessary).
