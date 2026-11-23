# INEQ-R8-hypercontract — can Gross hypercontractivity LOWER-bound the edge overlap `f₁²=|⟨0|O|G⟩|²`?

> 2026-11-23 · The positivity-improving-semigroup lever, attacked head-on. **Verdict: hypercontractivity is real and *explicitly available per block*, but it cannot lower-bound the overlap.** It is an `L^p→L^q` *contraction* — every matrix element it bounds is bounded **above**. Sign persistence is Perron–Frobenius order-positivity (free, sign-only, no constant); its quantitative strengthening *is* the gap. Globally at β=2.5 it collapses by the *same* two-road wall as H2: the direct BE/LSI route dies at Savvidy, the block→global route dies at inter-block coupling. Since LSI ⇒ gap, it is the *harder*, not easier, lever. No fabricated proof.

**Setup.** `P_t=e^{−tH}` heat semigroup on `L²(μ_Λ)`; `f₁²=ν({m₁²})` the edge atom. In the transfer-matrix identification the ground state is `g₀>0` (PF) and the vacuum test function is `1`, so `⟨0|O|G⟩=∫O g₀ dμ`.

**Gross's equivalence, and its direction.** LSI `Ent(f²)≤2C E|∇f|²` ⇔ hypercontractivity `‖P_t f‖_q≤‖f‖_p`, `q−1=(p−1)e^{2t/C}`. Composed with Hölder, `|⟨g,P_tf⟩|≤‖g‖_{q'}‖P_tf‖_q≤‖g‖_{q'}‖f‖_p`: **every** correlator/overlap bound is an **upper** bound. The overlap `|⟨0|O|G⟩|²` sits on the small side; hypercontractivity can only certify `≤`.

**Per block: YES — constructed.** PROVE-H2's Holley–Stroock bound `Ent≤e^{C_B}C_G E|∇_⊥|²`, `C_B=2βN_p(B)`, is an LSI with `C_HS=½e^{C_B}C_G`, uniform in `U′,L`, and boundary data. Gross then gives, explicitly, `‖P_t f‖_q≤‖f‖_p` with `q−1=(p−1)e^{2t/C_HS}`. This is genuine block-local hypercontractivity — "does it hold per-block?" is **yes**, and it is the same object as H2.

**Why persistence cannot lower-bound `f₁²`.** (i) *Sign persistence* `P_t f≥0` (strict for `t>0`) is the maximum principle/order-preservation — PF, not hypercontractive — fixing the *sign*, never a *constant*. (ii) The *quantitative* strengthening `|⟨0|O|G⟩|≥c>0` needs `g₀≳c'`, i.e. the PF gap = the mass gap. (iii) The only correlator hypercontractivity could lower-bound is `C(t)=⟨O(t)O(0)⟩≥f₁²`, which upper-bounds the atom — precisely INEQ-T2's correlator-collapse (`e^{tm₁²}C(t)↘f₁²` from above).

**Global β=2.5 collapse.** Two roads, both closed as for H2: the direct BE/LSI route fails at Savvidy (`Hess S_W≤−1.1β`, `G″≈−0.165`; W6); the block→global localization fails on the inter-block influence `C(β)>1` (Zegarlinski; H2-GLOBALIZATION). And LSI ⇒ gap makes uniform global hypercontractivity **stronger** than the gap — not a weaker lever. (Savvidy alone is intra-block and absorbed by `C_HS`; the *global* failure is the coupling contraction, exactly H2's wall.)

**Sharpest line.** Hypercontractivity is the right-*shaped* positivity tool and genuinely holds per block — but it is an `L^p→L^q` upper-bound regime, so it upper-bounds the overlap, not lower; its quantitative content is the gap; and its global form dies where H2 does. It collapses to the same **direction/regime** wall as the other residue-duals.

*Refs.* Gross 1975; Bakry–Gentil–Ledoux; Holley–Stroock 1987; Zegarlinski CMP 1989; W6-BE-at-2.5; PROVE-H2; H2-GLOBALIZATION; INEQ-T2-variational-residue; EDGE-f₁² §2; LEVER-nonconvex-edge (direction/hypothesis).
