# INEQ-T1 — can Hankel / Christoffel–Darboux / Nevanlinna–Stieltjes structure LOWER-bound the edge atom `f₁²`?

> 2026-11-23 · Attack on the one open object (`f₁²=ν({m₁²})`) from the moment-matrix / orthogonal-polynomial / Stieltjes-transform side. **Verdict: No. This class is itself an "upper-bound machine." Every genuine theorem here — Christoffel–Darboux, Nevanlinna–Pick, Hankel PSD, Stieltjes inversion — bounds the atom from ABOVE, and the lower end of the admissible range is exactly 0 at every finite order.** The atom is a point-supported (singular) functional outside the span of the moment tests; no positivity datum can force it positive. No fabrication; the one valid lower bound is the *already-known* conditional level-gap lemma, which needs the non-moment datum `Λ`.

**Setup.** `ν≥0`, `supp ν⊂[m₁²,∞)`. Set `G(z)=−F(z)=∫dν(t)/(t−z)`, a Herglotz function (`Im G>0` in `ℂ⁺`). An atom `f₁²` at `m₁²` is a **pole** of `G`, and the jump/Stieltjes-inversion formula reads `f₁² = lim_{ε↓0} ε·Im G(m₁²+iε)`. The moments are the **Taylor germ at ∞**: `G(z)=−Σ M_k /z^{k+1}`. So the question is precisely: *can a residue of `G` at a finite interior point be lower-bounded by its Taylor data at ∞?*

**The one-direction theorem (Christoffel–Darboux).** Let `K_n(x,y)=Σ_{k≤n}p_k(x)p_k(y)` be the CD kernel (ON polynomials of `ν`), `λ_n(x)=1/K_n(x,x)`. Take `q(t)=K_n(t,x)/K_n(x,x)`, so `q(x)=1`. Then `∫q²dν=λ_n(x)` (reproducing property), while the atom contributes `f₁²·q(x)²=f₁²`. Hence
  **`f₁² ≤ 1/K_n(x,x) = λ_n(x)` for all `n`** — an UPPER bound; `λ_n` is built from `M_0,…,M_{2n}`. The max atom over the moment fiber is exactly this Christoffel value (SDP/duality). Same sign as every RP/PSD bound: atom ⇒ total larger. There is no companion `≥`.

**Lower end is 0 — explicit, at the physical edge.** On `[1,∞)`, edge `m₁²=1`:
- `ν_a=δ₁+δ₃`: `M_0=2`, `M_1=1+3=4`, `ν_a({1})=1`.
- `ν_b=` Lebesgue density `1` on `[1,3]`: `M_0=2`, `M_1=∫₁³s ds=4`, `ν_b({1})=0`.
Identical protected moments `(M_0,M_1)=(2,4)`, identical edge `1`, **atom `1` vs `0`**. The full range is realised: fixing the (interior) Hankel matrix leaves a convex spectrahedron of representing measures on which the atom functional attains both `0` and the CD maximum `λ_n`. So `inf = 0`, always. **No lower bound exists from moments alone.** (`ν_b` is determinate by Carleman yet shows the finite-data failure directly.)

**Nevanlinna–Pick adds nothing.** Values of `G` at finitely many points give a Pick-matrix PSD condition ⇒ an **upper** modulus bound on `|G|` ⇒ upper on the residue. Residue-`0` Herglotz interpolants (absolutely continuous near `m₁²`) are always admissible ⇒ lower end `0` again.

**Carleman/determinacy does not rescue it.** If `Σ M_k^{-1/(2k)}=∞` (plausible for YM), the FULL sequence determines `ν`, hence `f₁²`. But (i) uniqueness ≠ inequality — `f₁²` is *discontinuous* in the moments (interior vs boundary of the moment cone); (ii) the full sequence is not available from protected data; (iii) `m₁²` must be known a priori to even name the atom. So determinacy is a uniqueness remark, not a constructible bound.

**The conditional survivor (not new).** The residue–level-gap lemma (GAP-G2-edge) in `G`-language: `G(z)−f₁²/(m₁²−z)` is Herglotz with support `≥Λ`, giving `f₁² ≥ (ΛM_0−M_1)/(Λ−m₁²)`. It requires `Λ>m₁²` (a *strict* second level — a gap one rung up). Without `Λ`, "rest `≥m₁²`" is vacuous: density can pile just above `m₁²` with arbitrarily small atom and still match `M_0,M_1`. `Λ` is the extra, non-moment positivity.

**Sharpest line.** Moments are the Taylor germ of `G` at `∞`; the atom is a pole at finite distance, and extracting its residue is analytic continuation — the classically ill-posed step. Hankel/CD/Nevalinna–Pick implement only the well-posed direction, `f₁²≤λ_n`. The lower bound is `0`, it is *realised*, and the free field (`g=0`, `f₁²=0`) sits on it. **The atom weight is inaccessible to Hankel/Nevanlinna data; construction is impossible; the angle is upper-bound-only.**
