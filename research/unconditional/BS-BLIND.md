# BS-BLIND — Lower bound on the 0⁺⁺ two-gluon mass, blind

**Target.** 4D pure SU(N) YM, defined as a Wightman QFT with transfer matrix H ≥ 0, H|Ω⟩ = 0. Let P be the orthogonal projection onto the 0⁺⁺ (J^{PC}=0^{++}) sector and set
m := inf spec( H |_{P(Hilb ⊖ Ω)} ).
Prove m > 0 uniformly in UV cutoff Λ and volume V. The 0⁺⁺ glueball is the 1⁻⁻ ⊗ 1⁻⁻ composite; the operator tr F² carries it. Free massless field: tr F² = 0 and m = 0, so every purported tool must fail there (acid test).

This note checks four routes and asks honestly whether any gives a **lower** bound. A lower bound on an eigenvalue is rare; almost everything physical gives upper bounds or sign constraints. I flag direction at every step.

---

## (a) Bethe–Salpeter / bound-state equation

The homogeneous BS equation is [K^{-1} − V]φ = 0, with K = S⊗S the free two-gluon propagator and V the exchange kernel; the binding threshold is λ_min(K^{1/2} V K^{1/2}) = 1. Positivity of m requires this eigenvalue to lie below threshold.

The defect is structural. V is built from transverse-gluon exchange ∝ 1/k², positive-definite on physical transverse modes — **but** the background quadratic form is Hess = −D² + 2 i g [F_B, ·]. At a constant chromomagnetic background F_B the magnetic-transverse block has negative eigenvalues (Savvidy). The leading kernel is therefore sign-indefinite, and ladder/rainbow summation does **not** preserve positivity; K^{1/2}VK^{1/2} is not a positive operator. A negative quadratic mode is the one-loop **instability of the perturbative vacuum** toward magnetic condensation — not a physical tachyon. Controlling the sign of the resummed kernel to all orders *is* the mass-gap problem. **No lower bound. Direction: undetermined sign.**

## (b) Temple / Weinstein variational lower bound

Temple: for normalized ψ, E = ⟨H⟩, and a number β with β ≤ E₁ (next eigenvalue) and β > E,
E₀ ≥ E − (⟨H²⟩ − E²)/(β − E).
This *is* a genuine lower bound on the ground-state energy of a sector — but only once (i) β is a **lower bound to the next eigenvalue** and (ii) the variance is controlled.

In the 0⁺⁺ channel the "next eigenvalue / threshold" is the two-constituent continuum edge ≥ 2m₁, so β demands a **lower bound on one-particle masses** — circular. Near Savvidy the variance is uncontrollable, the quadratic form changes sign, β−E → 0, and the bound degenerates to −∞. **Direction: correct (lower), but vacuous without β − E > 0.** It *reduces*, not proves, the gap.

## (c) Reflection positivity of tr F²

RP gives ⟨θ A A⟩ ≥ 0, hence the Källén–Lehmann measure ρ(ds) on the 0⁺⁺ channel is **nonnegative**, and C(t) = ⟨tr F²(t) tr F²(0)⟩_c = ∫ e^{−s t} ρ(ds) is completely monotone (log-convex, alternating-derivative signs). These are *upper* bounds on C(t) and force monotone decay. They are fully compatible with supp ρ ∋ 0: a positive measure may pile up at s = 0. **RP yields no lower bound on inf supp ρ.** Concretely, the free field obeys RP with m = 0 — the acid test refutes any "RP ⇒ gap" claim. **Direction: constraints/upper bounds only.**

## (d) Confining string / area law

Area law ⇒ static potential V(R) = σR ⇒ flux-tube picture with glueball mass ~ c√σ. A rigorous lower bound m_{0⁺⁺} ≥ c′√σ would need a **rigorous flux-tube excitation-gap bound**. The effective-string (Lüscher-term) expansion is asymptotic order-by-order, not a bound; and the 4D area law itself is unproven for non-Abelian YM. **Direction: would be positive — but conditional on two open inputs.**

---

## Savvidy region: does any route cross it?

The indefinite Hessian is real and gauge-quotient-stable (negative transverse-magnetic mode survives). It is a statement about the quadratic approximation to the **effective action in a background**, i.e. about the instability of the perturbative vacuum toward condensation. It is *not* a statement about the transfer-matrix spectrum of gauge-invariant states. The physical m is a transfer-matrix eigenvalue; it is positive iff the nonperturbative theory confines with a gap. None of (a)–(d) lower-bounds it across Savvidy: (a) loses sign, (b) needs the threshold oracle, (c) is sign-blind, (d) is doubly conditional.

## Reduction to one inequality

All four routes bottom out at the **same** statement:

> **(I)** inf spec( H |_{P(Hilb⊖Ω)} ) ≥ δ > 0, uniform in Λ and V;
> equivalently, the 0⁺⁺ spectral measure has support bounded away from 0;
> equivalently, ⟨tr F²(t) tr F²(0)⟩ ≤ C e^{−δ t} with δ > 0 uniform.

- (a) needs (I) as the BS threshold; (b) needs it as β − E > 0 plus variance control it cannot supply near Savvidy; (c) never yields it — RP allows δ = 0; (d) yields it only granted area law **and** a rigorous flux-tube gap.

So the obstruction collapses to a **single** lower bound: the exponential decay rate of the 0⁺⁺ correlator is bounded below — exactly the Clay mass gap in the 0⁺⁺ channel.

## Verdict

| Route | Kind | Gives genuine lower bound across Savvidy? |
|---|---|---|
| (a) Bethe–Salpeter | open | No — kernel sign-indefinite |
| (b) Temple/Weinstein | **conditional** on (I) | Only if threshold β given |
| (c) RP / spectral measure | proved (for constraints) | No — sign-blind, δ=0 allowed |
| (d) area law / string | open | Conditional on 2 open inputs |

**Proved:** RP/positivity of the spectral measure; Källén–Lehmann existence. **Conditional:** Temple/Weinstein lower bound *given* the threshold inequality (I). **Open:** the inequality (I) itself, hence m_{0⁺⁺} > 0 unconditional, across the Savvidy indefinite-Hessian region.

No route here is a proof. The honest statement: the two-gluon 0⁺⁺ positivity is **equivalent** to the one inequality (I); every standard lower-bound arrow (Temple/Weinstein, positivity, string) either assumes it, is direction-blind, or is itself open. The missing step is precise: a uniform lower bound on the 0⁺⁺ spectral threshold that survives the sign change of Hess = −D² + 2 i g[F_B,·].
