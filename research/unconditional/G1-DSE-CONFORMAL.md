# G1-DSE-CONFORMAL — the pure-YM (N_f=0) DSE **does admit** a conformal scaling solution (κ≈0.595), so it is **not inconsistent** at N_f=0 — **but that solution is a gauge-fixed, truncated object about unphysical propagators; it is neither a physical pure-glue CFT nor a refutation of (G1)**

> **VERDICT (decisive, no QED).** Taking the Landau-gauge gluon + ghost Dyson–Schwinger equations (DSE) with the
> scaling ansatz `D_A(k) ∼ k^{−2+2κ}`, `D_c(k) ∼ k^{−2−κ}`, **the N_f=0 DSE ADMITS a scale-covariant ("conformal")
> solution** — the *scaling solution* with IR exponent `κ ≈ 0.595` (von Smekal–Hauck–Alkofer 1997;
> Fischer–Alkofer 2003). So the propagator DSE is **not inconsistent** with a gapless-looking IR, and this looks
> superficially like "a pure-glue CFT, contradicting (G1)". **It is not.** The scaling solution is `[C]`
> (truncation-dependent), refers to the **gauge-fixed, unphysical** gluon/ghost propagators, and the scaling
> gluon `D_A(0)=0` is a **positivity-violating / confinement signal, not a massless particle**. The physical
> gauge-invariant (glueball) spectrum is a *separate* object the propagator DSE does not compute. Hence the DSE
> **cannot settle (G1)**, and does **not** refute it: (G1) — "no physical strongly-coupled IR fixed point in pure
> YM" — is untouched. The honest status of (G1) is `[O]`, and the DSE scaling solution is a red herring for it.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric. Here the key labels: scaling solution `[C]`; physical interpretation `[O]`; (G1) refutation `[NG]`.

---

## §1 · Setup and the scaling ansatz

Landau-gauge gluon `D^{ab}_{μν}(k)=δ^{ab}(δ_{μν}−k_μk_ν/k²)Z(k)/k²` and ghost `D^{ab}(k)=−δ^{ab}G(k)/k²`. A
"conformal" IR means scale-covariant dressings `Z(k)∼k^{2κ}`, `G(k)∼k^{−κ}`, i.e. propagators
`D_A∼k^{−2+2κ}`, `D_c∼k^{−2−κ}`, with IR exponent `κ>0` (the ghost is "enhanced", the gluon "suppressed"). A free
CFT would have `κ=0`; the scaling solution has `κ≠0` — an **interacting** CFT-like IR.

## §2 · The scaling solution EXISTS at N_f=0 `[C]`

The coupled ghost + gluon DSEs with the standard "1+κ" ghost-gluon-vertex truncation reduce in the IR to a
transcendental for `κ`. It has a positive root for pure YM (`N_f=0`), `κ≈0.595` (SU(3), 4D;
von Smekal–Hauck–Alkofer 1997; Fischer–Alkofer, *Prog. Part. Nucl. Phys.* 51 (2003) 581; reviewed in
Huber, *Phys. Rept.* 879 (2020) 1). The associated gluon exponent is `κ_A=2κ`, ghost `κ_G=κ`, and the exponents
satisfy the IR power-counting + the (ghost-gluon) Slavnov–Taylor constraint. **So, at the level of the truncated
DSE, pure YM admits a conformal scaling solution — the DSE is NOT inconsistent at `N_f=0`.** (Witness solves the
representative leading IR balance and finds the positive root; the exact value is truncation-dependent.)

## §3 · Why this is NOT a physical pure-glue CFT, and does NOT contradict (G1)

1. **Gauge-fixed, unphysical.** The gluon and ghost are *not* in the physical (gauge-invariant, BRST-cohomology)
   Hilbert space. The scaling solution is a statement about Landau-gauge *propagators*, not the physical spectrum.
2. **`D_A(0)=0` is a confining signal, not a massless particle.** A vanishing gluon propagator in the IR is the
   hallmark of **positivity violation** (the gluon is not an asymptotic state) — the standard "scaling" picture of
   *confinement*, not of a gapless free/CFT excitation. The ghost's `G(k)∼k^{−κ}` divergence is likewise
   unphysical (no asymptotic ghost).
3. **The physical gap is a separate object.** The `0⁺⁺` glueball mass is a property of the gauge-invariant
   correlators; the propagator DSE never computes it. A scale-covariant `D_A,D_c` is fully consistent with a
   *gapped* glueball spectrum (the gluon is confined, not the lightest state).
4. **The scaling-vs-decoupling debate is truncation/scheme-dependent `[O]`.** The *lattice* (Aguilar–Binosi–
   Papavassiliou; Bogolubsky et al.) finds the **decoupling** solution `D_A(0)>0` (a "massive" gluon), not the
   scaling solution. Whether scaling or decoupling is the true IR is *not settled* (it depends on the vertex
   truncation, the Gribov–Zwanziger horizon, and the renormalization scheme). So even as a statement about
   propagators, "conformal scaling" is contested.
5. **Truncation caveat.** The IR exponent `κ` and even the *existence* of the scaling solution depend on the
   vertex ansatz (ghost-gluon vertex IR behavior, three-gluon vertex) and on imposing the Gribov–Zwanziger
   horizon condition. The scaling solution can be an artifact of a truncation that neglects the decoupling
   mechanism.

## §4 · Verdict on (G1)

**(G1): "pure YM has no physical nonperturbative IR fixed point."**
- The DSE **admits** a conformal *scaling* solution at `N_f=0` `[C]` — but it is gauge-fixed, truncation-dependent,
  and about **unphysical** propagators. It is **not** a physical gapless CFT.
- Therefore the DSE **does not refute (G1)** and **does not prove it**: it is simply blind to the physical
  spectrum. (G1) remains `[O]`.
- The DSE "conformal solution" is the **gauge-fixed cousin** of the "gapless" worry: it looks conformal in the
  propagators while the physical theory is confined. It is a red herring for the physical mass gap. Consistent
  with `STAR-NO-BZ` (no nonperturbative β sign theorem) and `ONLY-EXIT`.

*Lesson:* the DSE scaling solution is exactly the kind of object that *looks* like a pure-glue CFT and is not —
because it lives on unphysical, gauge-fixed, truncated propagators. It neither contradicts (G1) nor settles it.

*Refs:* von Smekal–Hauck–Alkofer, PRL 79 (1997) 3591; Fischer–Alkofer, *Prog. Part. Nucl. Phys.* 51 (2003) 581;
Huber, *Phys. Rept.* 879 (2020) 1; Aguilar–Binosi–Papavassiliou, PRD 78 (2008) 025010; Zwanziger, PRD 76
(2007) 125014; `STAR-SYNTHESIS.md`, `G1-NONPERT-BETA.md`.
