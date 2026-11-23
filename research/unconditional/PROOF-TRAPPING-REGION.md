# PROOF-TRAPPING-REGION — the explicit gauge BBS trap: constructible off the cone, and the cone *is* (★)

> **Verdict (判词).** **The region is now written down; its small-field half closes and its large-field half is (★).** The BBS/B–B–Y trapping region for 4D pure SU(N) is an explicit set — a coupling band, a remainder ball in a **coercivity-weighted Hessian (jet) norm**, and a cone-empty invariant inequality — on which `‖R_{j+1}‖_X ≤ e^{−c}‖R_j‖_X + Cg_j⁴` holds on the small-field half, `c` uniform in `j,a` at fixed `L`. On the large-field half coercivity `∫TrF²≥0` bounds the action's *value* (explicit tail `e^{−cΛ_{j+1}⁴/g_j⁴}`, ◇), but the contraction must transport the *second variation*, and `∫TrF²` is value-coercive, Hessian-indefinite (`−D²+2adF_B⪰−2g|F|`). Gauge invariance relocates but cannot delete it. The missing estimate is exactly **(★)**; at `k→0` it is H1. This spells out the trap PROOF-MULTISCALE-4D §2 asserts; that verdict stands.

**§1 The region — set and norm.** Chart `K_j=g_j^{−2}∫TrF²+R_j` (`θ=TrFF̃` total-derivative, drops), `R_j` on gauge-invariant dim>4 functionals, `Z_j` tracked; shell `[Λ_{j+1},Λ_j]`, `Λ_{j+1}=eΛ_j`. Fix `g*,r,θ>0` and cutoff `M`. The trap is
```
𝔗 = { (g_j,Z_j,R_j): 0<g_j≤g*,  g_{j+1}²=g_j²−2β₀g_j⁴+O(g_j⁶),
        ‖R_j‖_X ≤ r g_j²,  and  2g_j|F_{Λ_{j+1}}| < Λ_{j+1}² }.
```
`X` must be a **jet** norm, not a value norm. *Necessity.* Any `‖R‖=sup_A|R(A)|w(A)`, `w=e^{−cQ[F]/g²}≤1`, is saturated by `Q≥0` and blind to `δ²Q=−D²+2adF_B` having `λ_min=−2g|F|<0`: a functional can be value-small yet curvature-large. So take a BKAR/KP polymer norm weighted by the coercive action and dominating the shell second variation,
`‖R‖_X := sup_{|v|=1} sup_A e^{Q_{Λ_{j+1}}[F]/4g_j²}|δ²_{vv}R(A)|·ℓ(A)`,
the multiplier set by `δ²K_j=g_j^{−2}(−D²+2adF_B)+δ²R_j`. Forced; everything concentrates here.

**§2 The one-step contraction.** Split at the cone.
- **SF (cone empty, `2g_j|F|<Λ_{j+1}²`).** Retained-band rate `Λ_{j+1}²−2g_j|F| ≥ c₀Λ_{j+1}²>0`; free part contracts `e^{−tL₀}≤Ce^{−ct}` (L0, banked); the force is `∝(g_j²)²=g_j⁴`. Hence
  `‖R_{j+1}‖_X ≤ e^{−c}‖R_j‖_X + Cg_j⁴`,  `c` uniform (box gap `2π/L` a-fixed, P1-L0). **Closes.**
  *Preservation (one-loop box L1):* `g_{j+1}²<g_j²` and `Λ_{j+1}=eΛ_j`, so the cone shrinks while the shell edge grows — `𝔗` maps into itself at `j+1`.
- **LF (cone re-enters, `2g_j|F|>Λ_{j+1}²`).** Amplification `e^{t(2g_j|F|−Λ_{j+1}²)}`; then (◇), maximizing against the weight (peak at the boundary `|F|=Λ²/2g`, independent of `t`),
  `sup_{|F|>Λ²/2g} e^{−c|F|²/g²}e^{t(2g|F|−Λ²)} = e^{−cΛ_{j+1}⁴/4g_j⁴}`
  → `0` along the UV tower (`g_j→0`), an a-fixed factor at the box scale: the **value**-tail is summable.

**§3 Where Savvidy breaks it; gauge cannot escape.** (◇) bounds a value, but `e^{−c}` is a *Hessian* contraction: on the cone `δ²(g_j^{−2}∫TrF²)|_{vv}=−2g_j|F|<0` multiplies `R_j`'s gradient. The trap needs there
```
(★_LF)  Var_{μ_j}(Pᵀ∇S|a=0) ≤ E_{μ_j}[Pᵀ Hess S P|a=0] + W − δ   on {2g_j|F|>Λ_{j+1}²},
```
i.e. `G″_j≥δ` on the cone. Coercivity gives `E[H]≥0` — the *wrong* edge; `Var≥0` subtracts, so `G″≤E[H]` and no sign follows (PROOF-L3 §2; one-loop equipartition forces `Var=E[H]` typically, `G″→W`). **Gauge invariance relocates, not deletes:** the orbit-quotient transverse Hessian keeps 7–8% negative at the SU(3) crossover (W4-blocked-positivity); BRST fixes the dim-4 class `{TrF²,TrFF̃}` — counterterms, not the quadratic form (G1-BPHZ §4iii). *Structurally it must fail:* the density on the compact group is periodic, and a periodic concave function is constant (INVENT-CONVEXIFICATION) — usable convexity dies at the wrap `U≈−1`, precisely the Savvidy region. Savvidy is periodicity's face.

**§4 The missing estimate — the same (★).** `(★_LF)` **is** (★) of G1-BPHZ §3 / L3 §1: the proviso `√(2g_j|F_B|)<μ_j` is the SF condition, so the trap's open half equals the cone's open half. UV shells (`g_j→0`): (★) plausible — cone empties, Haar `W≥δ_W` carries the sign, LF value `→0` — blocked only by the unproven `Var` upper bound (Lemma A `β²N_□/λ_f`, `λ_f` open). Box/`k→0` (`Λ_{j+1}→0`): `Λ_{j+1}⁴/g_j⁴→0`, weight `→1`, `Var(G)→0`, `G″→E[H]`=bare IR curvature `⟺ m_eff²>0`=H1=G2. Trap closes at fixed volume **modulo (★)**.

**Sharpest line.** The trap is a ball about the asymptotically-free line and it is small-field-complete; every estimate but one is banked. The one — contraction on the cone — is the sign coercivity cannot give, and gauge symmetry, being periodicity in disguise (a periodic concave function is constant), can only move the cone, never empty it. The large-field bound controls what the action *is worth*, not whether it *bends the right way*; that bend is (★), and off the UV tower it is the gap.

*Refs.* PROOF-MULTISCALE-4D §1–4; W4-polchinski-target §1; P2-BBS-to-gauge §2–4; W3-marginal-BPHZ; PROOF-G1-BPHZ §3,§5; PROOF-L3 §1–2; G1-fixed-volume; GAP-G1-largefield §1,§3; A1-twosided-UV §1–4; PROVE-Gpp-positive Lemma A; PROOF-JACOBIAN-SAVES (`W≥N/12`); INVENT-CONVEXIFICATION (periodic no-go); W4-blocked-positivity (7–8% neg); Bauerschmidt–Brydges–Slade, LNM 2242 (2019); Balaban, CMP 109 (1987) 249; Savvidy, Phys. Lett. B 71 (1977) 133; Prékopa 1971.
