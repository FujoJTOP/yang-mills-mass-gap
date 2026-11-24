# AREALAW-GAP — does "4D area law ⟹ 0⁺⁺ mass gap" give G″(k→0)≥λ?

> 2026-11-23 · **Verdict: passes the acid test, dies on the channel test.** The area law is a genuine confining object (free field is perimeter-law, so it *fails* the area law ⇒ non-vacuous), but it is a **1-form** statement about the Wilson loop, |⟨W(C)⟩|≲K^{|C|}e^{−σ·Area}, while the target G″ is the **0-form neutral** plaquette correlator — the 0⁺⁺ gap. No theorem bridges them. The "implication" is a *co-derivation* where true (strong coupling) and **open** where needed (β=2.5); and "area law at every β below the transition" does not exist beyond strong coupling for non-abelian SU(N).

## What is actually true
**Strong coupling (β<β₀): unconditional, but co-derived.** Osterwalder–Seiler '78 — refined by PROOF-P3-UNCONDITIONAL (f₁²>0 for β<β₁=O(1/N)) — draw from **one** KP cluster expansion the area law, confinement (V(R)≥σR), and exponential clustering *together*, V/L-uniform. Both true; neither obtained *from* the other. The species is not the mechanism even here.

**Small λ / large N ('t Hooft regime): unconditional area law.** CNS (arXiv:2505.16585, Thm 1.2; arXiv:2509.04688, Thm 1.6): |⟨W(γ)⟩|≤K^{|γ|}e^{−c·Area(γ)} for λ in the strong-coupling regime λ<λ₀(d) (**⚠ 更正 2026-11-24**：此前 `λ₀=1/(8(d−1))` 未在 CNS 源文核实，真实阈值远小于此，见下). The gap there rides on the Durhuus–Fröhlich mass-gap condition — again an *independent* input, not a corollary. At β=2.5 the physical point is far beyond this strong-coupling regime — **not ~75× but ~40 orders of magnitude**（重核阈值量级 `β₀(d)≈10^{−10d}/d`，SU(3) d=4 ≈ `10^{−41}`；确切常数须读 Thm 1.2 正文，见 `TRANSPORT-RIGOROUS-NUMERICS.md`）. CNS still name "small N, large λ, d=4" — this point — as open.

## Peierls / contour / screening beyond strong coupling? No.
Peierls needs a **positive** contour/surface fugacity. For non-abelian SU(N) the inverse (area-law) series is **Weingarten**, whose local weights are **signed** (Euler-characteristic sign for SU(2)); positive c_r≥0 live only in the forward series. There is no FKG (compact ⇒ torsion ⇒ no lattice order), no order parameter, no 4D duality. CPS–CNS tame the signed surface sum only at small λ/large N and flag the positive regrouping as the missing ingredient. "Screening" is the physical *content* of the area law, not a tool that proves it.

## Channel mismatch — why the implication fails
The area law's transfer-matrix content is a gap in the **string / flux-tube (1-form) sector**: Ψ_{W(C)} has energy V(R)≈σR, and m_string=√σ is a 1-form excitation. The 0⁺⁺ glueball is a **0-form** neutral plaquette state. The measured m_{0⁺⁺}=c√σ (c≈3.5–4, Morningstar–Peardon) is a **physical identity, not a theorem** — the "area law ⇏ gap" already banked in PROOF-STRING-TENSION-STABILIZER / PROOF-CONFINEMENT-RANDOM-SURFACE. Area law bounds the charged flux sector; G″(k→0) is the neutral 0⁺⁺ correlator. **For the SAME channel the area law yields √σ — a different, string-sector state — so m_eff²≥λ does not follow.**

## 5-point scorecard
1. **Nonperturbative — PASS** (surface/cluster; no loop/1/N/g), but only β<β₀ or small λ.
2. **Sign-controlling LOWER bound — FAIL.** Bounds |⟨W⟩|, not G″; wrong channel. The only bridge to a G″ bound (m_{0⁺⁺}≥c√σ) is the missing step.
3. **Non-circular — FAIL.** Reaching G″ needs either the transmutation identity (circular: assumes glueball–string) or an OS-positive Savvidy-crossing RG (the G2/C6 wall).
4. **Non-SUSY + FAILS free field — PASS.** Free field: perimeter law, no area law, no gap ⇒ genuine discriminator. Necessary, not sufficient.
5. **Uniform m_eff²(k→0)≥λ — FAIL.** No λ for the neutral channel; strong-coupling σ>0 is the string gap √σ, not G″; β=2.5 open.

## Acid test
Free field has a perimeter law (no area law) and no gap, so the area-law object **fails for the free field** — it is *not* vacuous and stays silent/wrong-free there. It **passes the free-field acid test as a discriminator**. But that is not the target: the species needs the area law to *lower-bound G″*, and it does not.

## Exact missing step
A **proven neutral-channel bound** m_{0⁺⁺}≥c√σ (c>0), or an **OS-positive cluster/surface expansion crossing Savvidy (β≈2.47) yielding area law ∧ 0⁺⁺ gap simultaneously** — the same negative-correlation object G2 and C6 show does not yet exist.

**Sharpest line.** The area law survives the free-field acid test (free field is perimeter-law, gameless — it fails the area law) but dies on the channel test: it measures the electric flux tube (1-form, |⟨W⟩|, gap √σ) while G″ is the neutral 0⁺⁺ plaquette correlator (0-form). Every rigorous area law — Osterwalder–Seiler at β<β₀, CNS at small λ/large N — arrives *with* the gap from a shared expansion, never *implying* it, in a different sector; at the finite-N crossover β=2.5 the authors name the point open, and the bridge m_{0⁺⁺}=c√σ is a measured transmutation. The species turns a genuine object into a theorem only where the gap is independently true, never into the needed implication.

*Refs.* Osterwalder–Seiler, Ann. Phys. 110 (1978) 440; Cao–Nissim–Sheffield, arXiv:2505.16585 (Thm 1.2), arXiv:2509.04688 (Thm 1.6); Cao–Park–Sheffield, arXiv:2307.06790; Durhuus–Fröhlich, CMP 75 (1980) 103; Seiler, LNP 159 (1982); Morningstar–Peardon, PRD 60 (1999) 034509. Internal: PROOF-CONFINEMENT-RANDOM-SURFACE; PROOF-P3-UNCONDITIONAL; PROOF-STRING-TENSION-STABILIZER; GAP-G2-magnetic; C6-DEPENDENCE-CLASS.
