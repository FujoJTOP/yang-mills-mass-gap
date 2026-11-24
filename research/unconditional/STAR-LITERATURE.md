# STAR-LITERATURE — what is actually known about (★)

**(★)** = *a gapless interacting 4D Wilsonian IR attractor of an asymptotically-free gauge theory is
conformal.* This file surveys the literature **precisely**; citations verified by web search (URLs at bottom).
Labels: [T] theorem · [C] conditional/proved-under-assumptions · [O] open · [NG] no-go/false.

---

## §0 · Why (★) matters for the campaign

If (★) holds for **pure** YM, then the pure-YM β-function (no Banks–Zaks zero, §3) would force pure YM to be
**gapped** (a gapless attractor would have to be conformal = a BZ/fixed point, which pure YM has none).
So (★) is a non-measure route to the gap. Its status is the crux.

---

## §1 · scale ⟹ conformal in 4D — NOT fully proved, but strongly argued

- **Polchinski, *Scale and conformal invariance in QFT*, Nucl. Phys. B303 (1988) 226–236 [C/partial].**
  Completes the **2D** argument (Zamodolchikov: rigid scale ⟹ local/conformal under broad conditions). In
  **3+ dimensions he explicitly could find NEITHER a general proof NOR a counterexample** — i.e. 4D
  scale ⟹ conformal was left **OPEN** by the canonical reference.
- **Dymarsky–Komargodski–Schwimmer–Theisen, *On Scale and Conformal Invariance in Four Dimensions*,
  arXiv:1309.2921, JHEP 10 (2015) 171 [C].** Under unitarity they show infinitely many matrix elements vanish
  (a nontrivial *necessary* condition for conformality) and argue it is also *sufficient*, i.e. scale ⟹
  conformal in unitary 4D — but they **"discuss possible exceptions."** Not a clean theorem.
- **Dymarsky–Farnsworth–Komargodski–Luty–Prilepina, arXiv:1402.6322, JHEP 02 (2016) 099 [C].** Follow-up on
  scale/conformality and generalized free fields.

⇒ The "scale ⟹ conformal" half of (★) is **argued [C], not proved [T]**, for unitary 4D, with caveats.

---

## §2 · scale WITHOUT conformal — REAL, via limit cycles; this makes (★) FALSE in general

- **Fortin–Grinstein–Stergiou, *Scale without Conformal Invariance: An Example*, arXiv:1106.2540,
  Phys. Lett. B704 (2011) 74–80 [T-model].** An explicit **unitary** model in `D=4−ε` that is **scale- but
  not conformally-invariant**, with finite correlators — realized by a **limit-cycle** RG trajectory.
- **Fortin–Grinstein–Stergiou, *Limit Cycles in Four Dimensions*, arXiv:1206.2921, JHEP 2012 [T-model].** A
  limit cycle in a **unitary 4D gauge theory** (three-loop); proves **4D gauge β-functions are NOT gradient
  flows**.

⇒ **A gapless IR can be a scale-invariant limit cycle, not a CFT.** Therefore **(★) is FALSE as a general
statement** about asymptotically-free gauge theories: an AF theory with the right (multi-coupling) structure
can be gapless and scale-invariant without being conformal.

---

## §3 · BUT limit cycles require ≥2 couplings — pure YM is excluded

- A **single real coupling** `g` gives a **1D** β-flow `dg/dln μ = β(g)`; a 1D autonomous flow has **no closed
  recurrent orbit** (no rotation structure) ⇒ **no limit cycle**. Limit cycles need ≥2 couplings (a 2D+
  vector field with a rotation).
- **Pure YM has ONLY the gauge coupling** `g` (no matter, no tunable second coupling). So the
  Fortin–Grinstein–Stergiou escape **does not apply to pure YM**.
- Pure-YM 2-loop β: `β(g) = −β₀g³ − β₁g⁵`, `β₀ = 11N/(48π²)`, `β₁ = 17N²/(384π⁴)`, **both > 0 for N≥2** ⇒
  `β(g) < 0` for **all** `g>0` ⇒ **no Banks–Zaks zero, no conformal window** (no matter to tune).
  (witness `w7/star_literature.py` confirms for SU(2),SU(3): no positive root of β.)

⇒ For **pure YM**, (★) reduces to: *a one-coupling AF flow with `β<0` for all `g`, IF gapless in the IR, must
reach a scale-invariant attractor (a fixed point, possibly `g=∞` asymptotic safety) — hence a CFT.*

---

## §4 · The answer to the posed question

**Has anyone proved "a gapless AF 4D gauge theory is gapped or CFT"? — NO [O].**
- (★) is **FALSE** for AF gauge theories **with matter/multiple couplings** (limit cycles, §2).
- The "scale ⟹ conformal" half is **[C]** (Dymarsky et al.), not **[T]**.
- For **pure YM** the limit-cycle escape is closed by the **single-coupling no-cycle** fact (§3), but
  no-limit-cycle + scale⟹conformal(4D, caveated) is **not a proof** of (★) for pure YM.
- **Closest rigorous theorem: the a-theorem** (Komargodski–Schwimmer, JHEP 12 (2011) 099, arXiv:1107.3987
  [T]) — `a_UV ≥ a_IR`, via the dilaton NG boson + trace anomaly + S-matrix unitarity. It **constrains the flow
  but does not force the gap.**

---

## §5 · The other named topics

- **'t Hooft anomaly matching for scale anomalies:** the trace anomaly `T^μ_μ = aE₄ − cW²` is **not** an 't
  Hooft anomaly (no global symmetry to match); the a-theorem (KS) is the rigorous analogue. There is **no**
  standard "scale-anomaly 't Hooft matching" that forces a massless sector.
- **Dilaton theorem / no spontaneous scale breaking without a dilaton:** the dilaton **is** the NG boson of
  spontaneously broken *conformal* symmetry (KS build their proof on it). But scale symmetry is **anomalous**,
  so spontaneous breaking of the *anomalous* scale symmetry is subtle; and a massless dilaton is itself a
  gapless mode (it is **not** the `0⁺⁺` glueball). So "gapless via a dilaton" is not a gap.
- **Conformal-window absence for pure YM:** trivial — pure YM has no matter, hence no Banks–Zaks window
  (§3): `β<0` for all `g`.

---

## §6 · Verdict for the (★) attack

1. **(★) as stated is FALSE** (limit cycles) — it must be **restricted to pure YM** (single coupling).
2. **Restricted to pure YM, (★) is UNPROVED [O]**, but the only general escape (limit cycles) is **closed by
   the single-coupling no-cycle fact** [C], and the "scale ⟹ conformal" half is **[C]** (Dymarsky et al.).
   So the honest target is: *prove that a one-coupling AF flow with `β<0` for all `g`, if gapless, has a
   scale-invariant (hence conformal) attractor* — combining §1 [C] with §3 [C].
3. **No literature proves the gap via (★).** The nearest results are the a-theorem [T] (flow constraint) and
   the 4D scale⟹conformal argument [C]. **(★) for pure YM is a genuine, sharp, unproved target — not a known
   theorem and not a known counterexample.**

---

## §7 · Citations (verified)

- J. Polchinski, *Scale and conformal invariance in quantum field theory*, Nucl. Phys. B303 (1988) 226–236.
- A. Dymarsky, Z. Komargodski, A. Schwimmer, S. Theisen, *On Scale and Conformal Invariance in Four
  Dimensions*, arXiv:1309.2921, JHEP 10 (2015) 171.
- A. Dymarsky, K. Farnsworth, Z. Komargodski, M. Luty, V. Prilepina, *Scale invariance, conformality, and
  generalized free fields*, arXiv:1402.6322, JHEP 02 (2016) 099.
- J.-F. Fortin, B. Grinstein, A. Stergiou, *Scale without Conformal Invariance: An Example*, arXiv:1106.2540,
  Phys. Lett. B704 (2011) 74–80.
- J.-F. Fortin, B. Grinstein, A. Stergiou, *Limit Cycles in Four Dimensions*, arXiv:1206.2921, JHEP 2012.
- Z. Komargodski, A. Schwimmer, *On Renormalization Group Flows in Four Dimensions*, arXiv:1107.3987,
  JHEP 12 (2011) 099.
- (Critique) *On ε-conjecture in a-theorem*, arXiv:1110.2586 (notes KS relies on dilaton decoupling).

*URLs (web-search verified):* arXiv.org/abs/1309.2921 ; arxiv.org/abs/1106.2540 ; arxiv.org/abs/1206.2921 ;
arxiv.org/abs/1107.3987 ; arxiv.org/abs/1402.6322.

*Labels:* [T] a-theorem, Fortin–Grinstein–Stergiou models; [C] scale⟹conformal 4D, single-coupling no-cycle;
[O] (★) for pure YM. **No QED; Clay 0%.**
