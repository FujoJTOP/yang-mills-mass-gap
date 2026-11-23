# PROVE-H2-uvgap — the UV fiber block: an explicit Poincaré/LSI gap

**Verdict.** H2 **in the block-local form the conditional theorem states and uses** (a fixed finite block; gap independent of the coarse field U′, volume L, spacing a and boundary data) is **proved**, with an explicit constant below. The two stronger demands are **open** and are pinned to the exact failing piece: **uniformity in block size b** (constant decays like e^{−2β·O(b⁴)}), and **globalization** of a block-local gap to the full fiber at β=2.5 (Zegarlinski localization needs a contraction that fails). Neither obstruction is the k→0 mode: it lives in U′=PU and never enters ∇_⊥, so H2 is genuinely easier than the gap.

## Fiber conditional
Fix U′. The fiber ξ∈G^{n_B} carries the block plaquettes not fixed by U′; each boundary plaquette sees U′ as a *fixed* link:
$$d\mu(\cdot|U')\propto e^{-\beta H_{\rm fib}(\xi;U')}\,d{\rm Haar}(\xi),\quad H_{\rm fib}=\sum_{p\subset B}\bigl(1-\tfrac1N{\rm Re\,tr}\,U_p\bigr),$$
and every summand ∈[0,2] **for every U′**.

## Theorem (proved; uniform in U′, L, a, boundary)
Reference = product Haar on the block. Since each plaquette summand lies in [0,2] regardless of U′ and ξ,
$${\rm osc}_\xi\bigl(\beta H_{\rm fib}\bigr)\le 2\beta\,N_p(B)=:C_B,$$
independent of U′, L, a. Holley–Stroock (its native LSI form, then LSI⇒Poincaré) gives, with C_G = LSI constant of normalized Haar on G (explicit, O(1)),
$${\rm Ent}_{\mu(\cdot|U')}(f^2)\le e^{C_B}C_G\,E[|\nabla_\perp f|^2],\qquad
{\rm Var}_{\mu(\cdot|U')}(f)\le e^{C_B}C_G\,E[|\nabla_\perp f|^2].$$
Hence H2 holds with
$$\boxed{\Lambda_{\rm UV}(b)=e^{-2\beta N_p(B)}/C_G>0,}$$
uniformly in U′, L, a and **all** boundary data (coarse *and* exterior-fiber), for each block. At b=2a, N_p(B) is a fixed integer; at β=2.5, Λ_UV ≥ e^{−5N_p}/C_G — finite and explicit though not optimized.

## Why the *natural* routes fail, and why they need not succeed
**Bakry–Émery Γ₂** is false pointwise: min(Ric+Hess S_W)≈−2.3 at β=2.5, Hess S_W reaching −1.12β, and 90–97% of the negative directions sit **in the fluctuation block itself** (W6-BE). **Dobrushin contraction** needs β·(plaquette-degree)<1: false at β=2.5. Neither is needed: a measure can have a gap without pointwise convexity, because the fiber is a *bounded perturbation of Haar*. That is exactly what Holley–Stroock exploits — and why H2 is provable though BE is not.

## H2 is easier than the gap — precisely
P is low-pass, so the k→0 mode is carried entirely by U′ and is **frozen** under fiber variation; ∇_⊥ annihilates it (given H3 faithfulness). The H1 wall (+0.043−0.208=−0.165, VAR-FINAL) is therefore structurally **absent** from H2. The residual difficulty is short-distance Savvidy concavity — a UV, not IR, fact.

## Honest scope
(i) Block-local gap, uniform in U′, L, boundary: **proved**. (ii) Uniform in b at β=2.5: **open** — UV obstruction. (iii) Block→global-fiber localization: **open** — Zegarlinski's one-site-to-global step fails the contraction at β=2.5; a certified b=2a interval-arithmetic gap plus an independent localization lemma would be needed. The conditional theorem's H2 is block-local, so the obligation it names is discharged.

*Refs.* Holley–Stroock 1987; Bakry–Émery; W6-BE-at-2.5; VAR-FINAL; P4.
