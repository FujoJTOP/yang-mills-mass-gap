# D-BBS-BACKGROUND-EFFPOT — build the trap on the background effective potential, not the fluctuation Hessian: Re V_eff has a scale-free positive curvature at its own minimum, but the tachyon sits on that mode, the thermal region is concave, and the free field shares the convexity

> 2026-11-23 · constructive. Successor to C-BBS-FORMULATE / CONSTRUCTION-METHODS (c). The Hessian route asks `δ²K_j = −D²+2ad F_B ⪰ δ` pointwise; the Savvidy cone `|k|<√(2g|F_B|)` makes it false. **New mechanism:** integrate the fluctuations out and build the trap on the **background effective potential** `Re V_eff(B)` (`B` a homogeneous chromomagnetic background), whose second variation in `B` is a scalar curvature, not an operator eigenvalue.
>
> **判词 (verdict).** Three facts are 真, one fatal. **真 (b):** `Re V_eff(B)` is one-loop computable and its curvature at its Savvidy minimum `B*` is a pure positive number `V_eff″(B*)=2c`, independent of `μ,g,L` — a scale-uniform `δ` on the condensate side. **真 (d):** the imaginary part is a *power* `Im V_eff ∝ (gB)²` (an `O(g²)` effect, `~19%` at β=2.5), not an exponential. **Fatal (c):** the tachyon is the `k=0` mode *orthogonal to `B`*, so `Re V_eff` is the real part of an unstable potential whose positive curvature along `B` is orthogonal to the negative direction; `V_eff″(B)<0` on `0<B<B*/e` — the measure-dominant thermal region at β=2.5 is exactly where it is **concave**; and no inequality transfers `Re V_eff″(B)` to `G″(a₀)`. **Acid test fails:** a free photon has `V_eff(B)=B²/2`, convex, so the curvature is positive for the *free* field — it measures background-field energy, not a gap. **Route gives a 真 computable real-part curvature but no uniform `δ` for the coarse mode; residual = H1.**

## 1 · (a) The background-field effective action
Split `A=B+a`, `B` a constant Cartan chromomagnetic field (`DF_B=0`), background-gauge fix. One loop, transverse:
`Γ_1[B] = S_W[B] + ½ Tr ln(−D_B² + 2ad F_B) − Tr ln(−D_B²)  (ghosts).`
For constant `B`, the renormalized real part is the Savvidy potential
`Re V_eff(B) = B²/(2g²) + c B²[ln(B/μ²) − ½]`, `c` the one-loop coefficient `∝β₀` (SU(2): `c=11/96π²`),
and `Im V_eff(B) ∝ (gB)² > 0` from the negative rung of the Landau ladder. The proposed trap replaces "sign `H_ff` pointwise" by "lower-bound `V_eff″(B)` on background space."

## 2 · (b) Real-part curvature — positive, but only on the condensate side
`V_eff′(B)=B/g² + cB[2 ln(B/μ²)−1]`, `V_eff″(B)=1/g² + c[2 ln(B/μ²)+1]`. The stationary point `B*=μ²e^{½−1/(2cg²)}` (exponentially small, `~Λ²`) gives
`V_eff″(B*) = 2c > 0`, a pure number — `11/48π²≈0.023` for SU(2), **independent of `μ`, `g`, `L`**: a real scale-uniform `δ` at the minimum. But `V_eff″` increases in `ln B`, and `V_eff″<0` for `0<B<B*/e` — the real part is **concave near `B=0`**, `V_eff″→−∞`. At β=2.5 the typical/thermal background sits at small `B`, the concave region. So "`≥δ` on the *relevant* background space" is **false**: the positive curvature lives on the condensate side, which the `e^{−S}` measure starves. The trap holds on the wrong side of the measure.

## 3 · (c) No transfer to the coarse mode
Two obstructions, both exact.
**(i) Wrong order parameter.** `B` is the homogeneous *field strength*; the coarse mode is `a₀=⟨v,log U⟩`, the volume-averaged log-link. `k0-CONSTRAINT-EFFPOT` proves `V_eff(α)` in `α` *is* `G″(0)` — an identity. `V_eff(B)` in `B` is **not**: no inequality `a₀ ≥ f(B)` exists. Different projections of `A`; one bounds the other nowhere.
**(ii) Tachyon at `k=0`, orthogonal to `B`.** The Savvidy negative mode `λ_min=−g|B|` (rung `n=0,s=−1,k_∥=0`) is a *constant* transverse-colour fluctuation — it lives on the very homogeneous sector `V_eff(B)` controls, but **orthogonal to the `B` amplitude**. So `Re V_eff(B)` is the real part of an unstable potential whose curvature along `B` says nothing about the negative direction. The full background Hessian `Γ″[B]=H[B]+O(g²)Tr(H⁻¹H″)` carries the indefinite `H[B]` *classically*; the one-loop term is `O(g²)` against an `O(1)` indefinite classical part — the FORGE5 scaling no-go. The pointwise indefinite Hessian is **not bypassed**; it re-enters as the classical term of `Γ″[B]`.

## 4 · (d) The imaginary part — a power, not an exponent
`Im V_eff ∝ (gB)² · 𝒪(1)` is a *power*: at β=2.5 SU(3), `g²/4π≈0.19`, negative-mode weight `E[neg]/E[pos]≈0.94`. So it is **not kinematically suppressed** — it reports the homogeneous saddle is unstable (Nielsen–Olesen), with decay rate `∝(gB)²` per volume. Its **only** suppression is the `e^{−S}` weight in the full path integral — exactly the dynamical `W>Var−E[H]`, i.e. **H1**. The tachyon blocks by making the homogeneous saddle *non-controlling*, pushing the burden back onto a constrained fiber variance: same lock.

## 5 · Ledger, acid test, missing step
| leg | status |
|---|---|
| `Re V_eff″(B*)=2c>0`, scale-free | **真定理** (one-loop real part) |
| `≥δ` on the *relevant* (thermal) background space | **FALSE** — concave for `B<B*/e` |
| transfer `Re V_eff″(B)→G″(a₀)` | **no-go** — wrong order parameter; tachyon at `k=0` |
| `Im V_eff ∝ (gB)²` power, `O(g²)` | **真**, but not exponential |
| tachyon suppression | `=e^{−S}` measure `=` **H1** |
| acid test (free photon `V_eff=B²/2`) | **FAIL** — free field is convex too |

**Acid test.** Free massless photon: `c=0`, `V_eff(B)=B²/(2g²)`, `V_eff″=1/g²>0` everywhere. The route returns a *positive curvature for the free field* — it certifies "a background field costs energy," not "the zero mode is massive." Vacuity, in the project's acid-test sense.

**Exact missing step.** The uniform coarse curvature `G″(a₀)≥δ` is still `Var_ν(a)≤1/λ` — **H1**. The background effective potential does not supply it; it relocates it to the wrong variable and the wrong measure-region.

**Sharpest line.** `Re V_eff(B)` is the real part of an unstable potential: its positive curvature `2c` runs along `B` while the tachyon sits at `k=0` *orthogonal* to `B`; the measure-dominant small-`B` thermal region is concave; the free photon shares the convexity; and the imaginary part is an `O(g²)` power. The background-effective-potential trap gives a 真, computable, scale-free real-part curvature at `B*`, but **no transfer to `G″_j`** — it blocks exactly where `B` meets `a₀`, at `k=0`. Residual = H1.

*Refs.* C-BBS-FORMULATE (`G″_j=EH−VarG+W≥δ`); C-BBS-MARGINALLOG §3; D-BBS-MEASURE-TRAP; CONSTRUCTION-METHODS (c, blocker); k0-CONSTRAINT-EFFPOT (`V_eff(α)=G″`); PROOF-COARSE-POINCARE (H1); PROOF-SAVVIDY-EXACT-SIGN (`λ_min=−g|B|` at `k=0`; `E[neg]/E[pos]≈0.94`); PROOF-SAVVIDY-REFORMULATION-CONVEX (invariant inertia); FORGE5-CONSTRUCTIVE (`δ²R=O(g²)` vs cone `O(μ²)`); E1-FRG (convexity ≥0, strict floor = H1). External: Savvidy, *PLB* **71** (1977) 133; Nielsen–Olesen, *NPB* **144** (1978) 376; Leutwyler, *NPB* **179** (1981) 129; Abbott, *NPB* **185** (1981) 189.
