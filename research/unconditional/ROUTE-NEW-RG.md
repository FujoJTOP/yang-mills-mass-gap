# ROUTE-NEW-RG — modern rigorous RG at the crossover: **no modern result gives the gap at β≈2.5**; the Savvidy-indefinite signed measure is the wall, and the exact missing estimate is the uniform large-field bound = H1

> **VERDICT (no QED).** Every modern rigorous RG for 4D gauge is a **small-parameter** method —
> strong coupling (`u(β)` small), weak coupling / UV (`g²` small), or subcriticality — and at the
> crossover `α_s~O(1)` **no such parameter exists** (witness: `u(2.5)≈0.51`, `g²(2.5)=1.6`, both `O(1)`).
> The best modern rigorous result at β=2.5 is **none that yields the mass gap**; the closest objects are the
> strong-coupling gap (`β<1/48`) and Balaban's `β→∞` UV control, and neither reaches 2.5. The **exact
> missing estimate** is a uniform large-field bound on the Polchinski flow with the **signed indefinite
> (Savvidy) effective action** — the `D`-term sub-criticality `κ_bare > Σ`, i.e. `χ≤1/λ` = `H1`.

---

## §1 · What reaches where

- **Strong coupling `β<1/48` `[T]`:** Osterwalder–Seiler 1978; Shen–Zhu–Zhu, *CMP* 400 (2023) (cluster/LSI, `L`-uniform gap). Rigorous — but `β<1/48`, ~25× in `u` below 2.5.
- **UV / `β→∞` `[T]`:** Balaban (Comm. Math. Phys., 1980s) — UV stability / renormalized trajectory of 4D lattice gauge. UV control, **not** an IR gap.
- **Modern constructive RG `[T, small-parameter]`:** Bauerschmidt–Brydges–Slade (RG method; weak coupling / small field); Bauerschmidt–Bodineau–Dagallier (Polchinski / SPDE for `φ⁴`). All require a small parameter.
- **Singular SPDE / regularity structures `[T, subcritical]`:** Chandra–Chevyrev–Hairer–Shen (3D YM–Higgs stochastic quantization). 4D gauge is **critical** (`[g²]=0`, asymptotic freedom); regularity structures need **subcriticality** ⇒ no direct control; the counterterms are the `u=2` renormalon.
- **Area law `[T]`** (Aizenman–Duminil-Copin for finite groups; Cao–Nissim–Sheffield per corpus): area law, **not** the transfer-matrix `0⁺⁺` gap.

## §2 · The best modern rigorous result at β=2.5

**There is none that gives the gap.** The two rigorous ends are SZZ (`β<1/48` gap) and Balaban (`β→∞` UV); the modern RG (BBD/BBS) reaches only where a small parameter exists. At `α_s~O(1)` the small parameter vanishes (witness), so no constructive-RG convergence is available.

## §3 · The exact missing estimate (= H1)

In the Polchinski-flow / BBS framework the controlling step is a **uniform large-field bound** on the effective action. The positive / ferromagnetic large-field control (BBS) **fails** for Yang–Mills because the Savvidy Hessian is **indefinite** (`λ_min=−g|B|`). Precisely, one needs a rigorous bound
> `κ_bare > Σ` (the `D`-term sub-criticality), uniformly in `(β,L,a)` at the crossover — i.e. `χ≤1/λ` `[O]` = `H1`.

So a "modern rigorous RG at 2.5" is blocked by exactly the same uniform large-field/Savvidy estimate the campaign pinned elsewhere: **H1 in the RG-flow coordinate.**

*Labels:* `[T]` survey pieces; `[NG]` no modern RG reaches 2.5 with a gap; `[O]` exact step = the uniform large-field/Savvidy bound = `H1`. Witness: `w7/route_new_rg.py`.
