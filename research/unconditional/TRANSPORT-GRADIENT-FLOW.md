# TRANSPORT-GRADIENT-FLOW — can the Lüscher gradient flow be the transport machine?

> 2026-11-24 · The campaign's last gap is a **nonperturbative, sign-controlling LOWER bound** on `m₁` (or the 0⁺⁺ residue `f₁²`), `(L,a)`-uniform, that **fails on the free massless field**. This note asks whether the **Lüscher gradient flow** — the one carrier never run as a *transport machine* — supplies the missing monotone/convex **differential inequality in flow time `t`**.
>
> **Labels.** `[T]` 真定理 · `[C]` 条件定理 · `[NG]` no-go (named class) · `[O]` 精确缺失步 · `[V-here]` my own elementary checks / numerics this session · `[R]` external (cited, not re-derived). **No QED on anything unproven; Clay stays open.**

---

## Verdict up front — `[V-here]`

> ### **NOT a new species. The flow supplies exactly one monotone differential inequality — complete monotonicity — which *is* reflection positivity / the project's `B6`→`H1` positivity wall, and it integrates to an UPPER bound on `m₁`, not a lower bound. The only gapped/gapless discriminator (the correlator decay rate) is a `t`-independent spectral identity. REDUCES TO `B6` + `H1`.**

Three statements, each sharp:

1. **A monotone differential inequality does exist, but its sign is fatal.** The flowed energy `E(t)` is a Laplace transform `∫h(u)e^{−2tu}du` of a positive measure (`h≥0`), hence **completely monotone**: `(−d/dt)^k E ≥ 0` for all `k`, i.e. `E` decreasing (`dE/dt≤0`), convex (`E″≥0`), log-convex. Equivalently the **flow-time effective mass** `m_flow²(t) := −½ d/dt ln E(t)` (a tilted spectral *mean*) satisfies
   `d/dt m_flow²(t) ≤ 0`,  and  `m₁² ≤ m_flow²(t)`  for every `t>0`,  with `m_flow²(t) ↓ m₁²`.
   The inequality *marches DOWN to* `m₁²` from **above**: integrating it from `0` to `T` gives `m₁² ≤ m_flow²(T)` — an **UPPER** bound. The lower-bound direction is not produced, and the only `t` where equality holds is `t=∞`, whose uniform control *is* the gap (circular). `[T]`
2. **Every positivity monotone is vacuous on the free massless field.** For the free massless field `t²E(t) ≡ 1/(64π²) > 0` (a positive constant, `[V-here]` numerics §3). Any lower bound built from `E(t)`, `tE(t)`, `t²E(t)`, `t²dE/dt`, or the action `S(t)` is therefore **satisfied by the gapless field** ⇒ it cannot imply a gap. **Acid test FAILS** (tool is vacuous). `[NG]`
3. **The one non-vacuous object is an upper bound.** The exponential decay rate of the flowed correlator `C(t,r)` **does** distinguish gapped from gapless (`→0` on massless, `→M` on massive, `[V-here]` §3), but (i) it is `t`-independent — the flow moves only the overlaps `c_n(t)=⟨0|E(t)|n⟩`, never the masses `m_n` (spectral theorem), and (ii) its log-derivative is again a spectral mean `≥ m₁`, i.e. an upper bound. It is the spectral representation, an identity, not a lever. `[T]`

**Net.** The flow is the cleanest UV *operator* statement in the corpus (`E(t)` is a renormalized local operator, Lüscher–Weisz `[R]`) and the empty IR *transport* statement: the flow-time variable is a second Laplace variable, so "complete monotonicity in `t`" is verbatim the Källén–Lehmann positivity the theory already has (`B6`), and the `k→0` coarse curvature that sets the mass is **flow-invariant** (`H1`). **Exact missing step unchanged: the edge atom `f₁²>0` / isolation `H*(r₀)`, `(L,a)`-uniform.**

---

## 1 · Rigorous facts to anchor on

**Existence / uniqueness.**
- **Lattice (the project's object):** the flow `∂_t U_μ = −g₀² ∂_{x,μ}S[U]·U_μ` is an ODE on the compact group `G^{4L}`; its vector field is smooth and global solutions exist for **all** `t>0` (no blow-up on a compact manifold). `[T]` (Lüscher, arXiv:1006.4518, §2.)
- **Continuum:** short-time existence/uniqueness hold (quasilinear parabolic). **Long time is subtle in 4D:** the energy bound `S(t)≤S(0)` gives a priori control but does **not** by itself exclude 4D concentration/bubbling (Struwe; Råde; Feehan; Waldron). The object the project needs is the **lattice-regularized** flow (global, above) and the **renormalized flowed operator at finite `t>0`**, which is well-defined. `[R]`
  *Correction to the naive line "long-time via energy bounds": the energy bound alone does not prevent continuum 4D singularities; compactness of the lattice gauge group does. We use the latter.*

**Flowed operators are renormalized.** For every `t>0` the flowed energy density `E(t)=¼G^a_{μν}(t)G^a_{μν}(t)` is a **finite renormalized local operator** after a `t`-dependent vacuum/contact subtraction; separated-point correlators `⟨E(t,x)E(s,y)⟩_c` are finite for `x≠y`, to all orders in perturbation theory. Small-flow-time expansion: `E(t) ~ c(t)\,Tr F²_R + (c-number) + O(t)`, with `c(t)→Z>0`. `[T]` (Lüscher–Weisz, arXiv:1101.0963; Harlander–Neumann 2016; Fritzsch–Ramos 2013 `[R]`.)

**Reflection positivity is preserved.** The flow is a deterministic map compatible with Euclidean reflection; the flowed theory satisfies OS reflection positivity with the flowed field as interpolating operator, so the flowed KL measure `ν_t ≥ 0`. `[T]` (Lüscher arXiv:1006.4518; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440 `[R]`.)

**The flow is descending.** `dS/dt = −2∫ tr(D_ν F_{νμ})² ≤ 0` in the continuum (lattice: same sign), so `S(t)` is monotone non-increasing. `[T]`

---

## 2 · Candidate monotones and the differential inequality (or its absence)

### 2.1 The scalar candidates and where they point

For a free scalar with mass `M` (implicit flow `∂_tφ=(∂²−M²)φ`, the `L²`-gradient flow) the **flowed energy density** is exact:
`E_M(t) = ∫ d⁴p/(2π)⁴ \, p² e^{−2t(p²+M²)}/(p²+M²) = e^{−2tM²}∫_0^∞ ds\, e^{−sM²}/(8π²(2t+s)³)`,
with closed massless value `E_0(t)=1/(64π²t²)`. Substituting `u=p²+M²`, `E_M(t)=∫ h(u)e^{−2tu}du` with `h(u)≥0`, `supp h⊂[M²,∞)` — a **Laplace transform of a positive measure**. Hence:

| candidate | behavior (free) | direction toward `m₁` | acid test on free massless | verdict |
|---|---|---|---|---|
| `S(t)` (action) | monotone `↓`; `[T]` | none — a dimension-4 positive number, free of the atom | `S>0` on massless ⇒ **vacuous** | `[NG]` |
| `E(t)` | `↓`, convex, `[T]` | none directly | `E>0` on massless ⇒ **vacuous** | `[NG]` |
| `A(t)=t·E(t)` | `↓` to `0` (both masses) | `→0`, produces no scale | `→0`, no gap info | `[NG]` |
| `t²·E(t)` | `=1/(64π²)` const (M=0); `→0` (M>0) | **positive where massless** | `const>0` on massless ⇒ **vacuous** | `[NG]` |
| `t²·dE/dt` | `→0` (both) | `→0` | `→0` | `[NG]` |
| `m_flow²(t):=−½ d/dt ln E` | tilted mean, `↓ m₁²` | **UPPER** `m₁²≤m_flow²(t)` | finite `>0` on massless, `→0` only as `t→∞` | `[NG]` as lower bound |

The table is the core negative: **every functional that is positive on all fields is positive on the gapless free field** (`t²E≡1/(64π²)>0` is the cleanest witness), so none can lower-bound a gap. `[V-here]`

### 2.2 The one genuine differential inequality — and its direction

For the free (and linearized-flow) theory `E(t)=∫h(u)e^{−2tu}du`, `h≥0` ⇒ **complete monotonicity** `(−d/dt)^k E(t) ≥ 0` for all `k` (Bernstein 1928/29; Widder 1941 `[R]`). The first two:
`E′(t) = −2∫u\,h\,e^{−2tu} ≤ 0` (`E` decreasing), `E″E−(E′)² = 4[∫u²w·∫w−(∫uw)²] ≥ 0` (`log E` convex), `w=e^{−2tu}h`. Consequently the flow-time effective mass
`m_flow²(t) = −½ E′/E = (∫u\,w)/(∫w)` (the `t`-tilted spectral mean)
is **non-increasing** `d/dt\,m_flow²(t) ≤ 0` and obeys `m₁² ≤ m_flow²(t)` for all `t>0`; it is `+∞` as `t→0⁺` (the UV renormalon of §2.2 below) and decreases to `m₁²` as `t→∞`. `[T]` (The interacting v.e.v. is *not* of this Laplace form — see the caveat in §4.)

> **This is the monotone differential inequality the task asked for — and it is exactly wrong-way.** It is a *one-sided* statement forced by `ν≥0` alone: by Bernstein's theorem "all derivatives alternate" **⟺** `ν≥0` **⟺** reflection positivity. So it carries **no dynamical content beyond positivity**, and it integrates **downward from above**: `m₁² ≤ m_flow²(T)`. A *lower* bound needs `m_flow²(T) ≤ m₁²`, which is **false** for any gapped field at finite `t` (strict inequality unless the spectrum is a single delta). Equality holds only as `t→∞`, and controlling `E(t)` uniformly on the scale `t≳1/m₁` **is** exponential clustering = the gap (circular). This is precisely the direction theorem of `TOOL-COMPLETELY-MONOTONE.md §2`.

Two structural obstructions to "integrate from `0`":
- **`t=0` is UV-divergent.** The unflowed `E(0)=∫p²/(p²+M²)\,d⁴p` diverges (`~1/t²` as `t→0⁺` from the flowed side; exactly `1/(64π²t²)` for the free massless field). `A(0)`, `E(0)`, `E′(0)` are the `a⁻⁴` renormalon — the `L3`/marginal wall, not a finite initial condition.
- **`t=∞` is the gap.** The only time the inequality is sharp is the one limit whose uniform control is the gap itself.

**Conclusion of §2.** `[T]` that a monotone `t`-differential inequality exists (complete monotonicity / `m_flow²` non-increasing). `[NG]` that it can yield a **lower** bound: it is the positivity (`B6`) structure and points strictly the other way.

---

## 3 · Free-field tests (mandatory) — script and results

Script: `D:\Dev\math\research\unconditional\w7\transport_gradient_flow.py` (run `py -3.12`, exact 1-D momentum/`s`-integrals, no lattice job — light, single process). It computes `E_M(t)`, `m_flow²(t)`, the flowed correlator `C_M(t,r)` and its `r`-decay effective mass, and `G″(t)=1/χ(t)`.

**(A) Positivity monotone is constant on the gapless field.**
```
 t     E0(t)      t^2 E0     t E0      t^2 dE0/dt   |  E1(t)      t^2 E1
 0.1  1.5831e-1  1.5831e-3  1.5831e-2  -3.1663e-2  | 1.1144e-1  1.1144e-3
 1.0  1.5831e-3  1.5831e-3  1.5831e-3  -3.1663e-3  | 9.5411e-5  9.5411e-5
10.0  1.5831e-5  1.5831e-3  1.5831e-4  -3.1663e-4  | 2.8525e-15 2.8525e-13
30.0  1.7590e-6  1.5831e-3  5.2771e-5  -1.0554e-4  | 4.8934e-34 4.4041e-31
```
`t²E0(t) = 1/(64π²) = 1.583143e-3` at **every** `t` (gapless). Any `λ>0` lower bound on `t²E` is satisfied by `m₁=0` ⇒ **vacuous**. `[V-here]`

**(B) The direction theorem, numerically — `m_flow²(t)` decreases to `M²` from above.**
```
 t     M=0        M=1
0.1  10.00000   11.63151
1.0   1.00000    2.24560
4.0   0.25000    1.34483
16.0  0.06250    1.09120
```
`M=1`: `m_flow²(t) ≥ 1 = M²` for all `t`, decreasing `→1` (an **upper** bound on the true mass²). `M=0`: `m_flow²(t)=1/t → 0`. Decreasing confirmed for both. `[V-here]` + `[T]`

**(C) The gapped/gapless discriminator is the `t`-independent `r`-decay rate.**
```
  M=0                      M=1
 r   t=.1  t=1   t=3      r   t=.1   t=1    t=3
 2  0.966 0.229 0.081     2  1.727  0.396 0.149
 8  0.250 0.249 0.200     8  1.193  1.140 0.580
```
Massless: `m_eff → 2/r → 0` (power decay, no exponential) — **tool dies, acid test passes at the detector level**. Massive: plateaus at `M=1` for `r≫√t` (visible at `t=0.1,1`; at `t=3` the asymptote needs `r≳4Mt≈12`, beyond the tabulated range, hence still rising — honest). In both cases the rate is **`t`-independent**: the flow changes overlaps, not masses. `[V-here]`

**(D) The `k→0` curvature is the project's `H1` object.**
`G″(t)=1/χ(t)=M²e^{2tM²}` (M>0), `0` (M=0); at `t=0`: `G″(0)=M²` (massive), `0` (massless) — the exact free-field acid-test value (`G″=0` on free massless). Two distinct flow factors, to keep straight: (i) the **operator-overlap / wave-function** factor `e^{2tM²}>0` rescales `G″` at finite `t` (a positive, sign-preserving normalization — moves no direction); (ii) the **coarse-mode damping multiplier** `e^{−tλ_a(k)}→1` on `k→0` — the gradient flow does not smear the gap mode. Hence the **coarse curvature** `G″(k→0)` (the `H1` object) is flow-invariant in the only sense that matters. `[V-here]`

**Acid-test summary.** Positivity monotones (`S,E,tE,t²E,t²dE/dt`) are **positive on the gapless free field** ⇒ vacuous (`[NG]`). The rate-level effective mass **is** a genuine detector (massless→0) but is an **upper bound** and `t`-independent ⇒ supplies no lower bound.

*(2D `CP^{N-1}` was listed as optional: skipped. It would not change the verdict — the direction theorem of §2.2 is theory-level (`ν≥0 ⇒` finite-`t` upper bound), and `CP^{N-1}` has no free massless-gauge acid-test role here; only a lattice job could add it, which the present reduction does not need.)*

---

## 4 · Reduction analysis — which wall, and why non-circularity fails

**The flow time is a second Laplace variable.** `E(t)=∫h(u)e^{−2tu}du` is the same mathematical object as the Euclidean correlator `C(τ)=∫e^{−sτ}dν(s)` of `TOOL-COMPLETELY-MONOTONE.md`; "complete monotonicity in `t`" is verbatim `KL/RP positivity (ν≥0) = B6`. The flow manufactures no new cone: it is the gradient twin of the Bakry–Émery `Γ₂=Hess S` (`H1-T6-geometric-flow.md`), the same second variation, the same Savvidy-indefinite Hessian `Hess S≤−1.1β` seen as flow *instability* rather than BE *negative curvature*. `[T]`

**The gap mode is flow-invariant.** Linearizing, the flow multiplies mode `k` by `e^{−tλ_a(k)}` with `λ_a(k)=(4/a²)Σ_μ sin²(k_μa/2)`. On the coarse mode `k₁=2π/L`, `e^{−t(2π/L)²}→1` as `L→∞`: the flow is the **identity** on the one mode that sets `m₁`. The `k→0` coarse Hessian `S″_t ≡ S″` is untouched (`PROOF-GRADIENT-FLOW.md §4`, `PROOF-FLOW-CONVEXITY-TRANSFER.md §3`). So the susceptibility/`G″(0)`=`H1` object is exactly what the flow cannot move. `[T]`

**Which no-go it hits.** `B6` (complete-monotone / positivity) **plus** `H1` (the `k→0` coarse Poincaré / inverse susceptibility). The squeeze is `f₁²/χ ≤ m₁² ≤ M₀/χ` (`G2A-INFSUPP-DIRECT.md`): the flow-time machinery reaches `χ=1/G″` and thus the **upper** end `M₀/χ`; the **lower** end needs `f₁²>0`, which is the atom — not a positivity output.

**Non-circularity.** Any attempt to read the lower bound from the `t→∞` rate requires `E(t)`/`C(t)` controlled uniformly on `t≳1/m₁`, which **is** exponential clustering = the gap. So the `t→∞` route is circular; the finite-`t` route is an upper bound; and the `t→0` route is the UV renormalon (`L3`). The three exits are the banked wall, restated in flow time. `[NG]` + `[O]`

**Interacting-theory caveat (honesty).** The free-field "tracking" of `M` by `m_flow²(t)` is a free-field accident: there is a single scale. In the interacting theory `⟨E(t)⟩` at large `t` approaches the **gluon condensate** `⟨¼G²⟩` plus a perturbative tail — a dimension-4 constant, **not** `e^{−2tm₁²}`; the physical mass is read only from the **correlator's `r`-decay**, which is the `t`-independent spectral identity. So even the free-field "it tracks the mass" does not transfer to the interacting case as a bound. `[T]/[V-here]`

---

## 5 · Precise verdict and the exact missing step

- **New species?** **No.** `[V-here]` The gradient flow is (i) a rigorous UV **operator** renormalizer (`E(t)`, `[T]`), and (ii) a **positivity** structure (`B6`) whose flow-time effective mass is a **decreasing upper bound** on `m₁`. It reduces to `B6` + `H1`, the two faces of the same "single-sided variance/positivity" species the whole campaign reaches (`docs/PROOF-SKELETON.md §1`).
- **Does it lower-bound the gap?** **No.** It lower-bounds nothing; at best it *upper*-bounds `m₁²` by a computable mean, tightening downward to it only in the circular `t→∞` limit. `[T]`
- **Acid test?** Positivity monotones: **FAIL** (positive on free massless — `t²E≡1/(64π²)`). Rate-level effective mass: passes the *detector* test (massless→0) but is an upper bound, not a lower-bound source. `[V-here]`
- **Exact missing step `[O]`.** The unprotected edge residue: `∃r₀>0, f₁²=ν({m₁²})≥r₀` **and** isolation (`m₁²<` next threshold), equivalently the lower squeeze `f₁²/χ ≤ m₁²` with `f₁²≥r₀`, `(L,a)`-uniform. This is `H*(r₀)` / the isolation wall — one rung above the gap — **unchanged** by the flow. The flow *defines* `f₁²(t)` in the continuum (its genuine achievement, `CONTINUUM-GRADIENT-FLOW-DEFINE.md`) but does not **lower-bound** it.

**One-line.** The flow brings the cleanest UV statement in the corpus and the empty IR one: in flow time the theory is completely monotone for exactly the reason it is in Euclidean time (`ν≥0`), so every finite-`t` scalar is a **decreasing upper bound** on `m₁`, the `k→0` mode that sets the mass is a flow **identity**, and the free massless field keeps `t²E(t)` a positive constant. Wrong half, sharply — the transport machine is not here; the key that must be forged is still the atom/isolation `H*(r₀)`.

---

## References

**External `[R]` (cited; not re-derived this session):**
- M. Lüscher, "Properties and uses of the Wilson flow in lattice QCD," JHEP **08** (2010) 071, **arXiv:1006.4518** (lattice flow existence/uniqueness, RP, flow coupling).
- M. Lüscher & P. Weisz, "Perturbative analysis of the gradient flow in non-abelian gauge theories," JHEP **02** (2011) 051, **arXiv:1101.0963** (flowed-operator renormalization, small-flow-time expansion).
- R. Harlander & P. Neumann (2016); A. Fritzsch & R. Ramos (2013) — small-flow-time / flow-coupling coefficients.
- K. Osterwalder & E. Seiler, Ann. Phys. **110** (1978) 440 (OS reflection positivity ⇒ `ν≥0`).
- S. Bernstein (1928/29); D. V. Widder, *The Laplace Transform* (Princeton, 1941) — completely monotone ⟺ Laplace transform of a positive measure.
- M. Struwe; J. Råde; P. Feehan; J. Waldron — 4D Yang–Mills flow singularities / bubbling (long-time continuum).
- G. Savvidy, Phys. Lett. B **71** (1977) 133; N. Nielsen & P. Olesen, Nucl. Phys. B **144** (1978) 376 (indefinite Hessian sign).
- Shen–Zhu–Zhu, CMP **400** (2023) 805, **arXiv:2204.12737** (`β<1/48` mass gap / log-Sobolev).

**Internal (banked verdicts this note classifies against):**
- `docs/BRIEF-FOR-MODEL.md` §2–§3 (five hard conditions; eleven-species death map); `docs/PROOF-SKELETON.md` §1 (the single "variance/positivity" species; `H1` vs `H*`).
- `research/unconditional/PROOF-GRADIENT-FLOW.md` (spectral support-fixing; multiplier identity on `k→0`).
- `research/unconditional/PROOF-FLOW-CONVEXITY-TRANSFER.md` (`λ↦λe^{2tλ}` sign-preserving; un-smoothing re-concaves).
- `research/unconditional/CONTINUUM-GRADIENT-FLOW-DEFINE.md` (flowed operator defines `f₁²(t)`; stops at isolation wall).
- `research/unconditional/G1-L3-VIA-FLOW.md` (flow renormalizes coupling/operator, not the sharp moment bound `L3`).
- `research/unconditional/H1-T6-geometric-flow.md` (flow rate = `−Hess S`; same as BE `Γ₂`).
- `research/unconditional/TOOL-COMPLETELY-MONOTONE.md` §2 (direction theorem `m₁²≤r(t)↓`, finite-`t` upper-only; `CM⟺ν≥0`).
- `research/unconditional/G2A-INFSUPP-DIRECT.md` (squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ`).
- `research/unconditional/INVENT-NEW-MONOTONE.md` (`A5` Perelman `W`-entropy: mirror-passes free, `Γ₂≈−2.0<0` on YM).
- Witness script: `research/unconditional/w7/transport_gradient_flow.py`.

*Status: `[T]` (monotone `t`-inequality = complete monotonicity / `m_flow²↓`; flow-invariance of `k→0`; `ν_t≥0`; operator renormalization) + `[NG]` (no lower bound from positivity/finite-`t` flow data; acid test fails) + `[O]` exact missing step `f₁²≥r₀` / isolation `H*(r₀)`, `(L,a)`-uniform. **Verdict: NOT a new species — reduces to `B6` + `H1`; Clay stays open.***
