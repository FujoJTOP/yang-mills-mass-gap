"""
pnt_largen_loop_witness.py  --  PNT-LARGEN-LOOP structural witness

Question: is the 4D Makeenko-Migdal (MM) loop equation CLOSED in W(C) (no DSE n-point
tower), and does it give a closed scalar gap equation for m_0++?  Where exactly does it
fail, and is the obstruction the same infinite tower as SLOT-DSE or structurally better?

Three checks (one exact numerical, two structural/demonstrative -- NOT Yang-Mills MC):

 (A) D=2 [T exact]: the signed area of a closed curve is a BOUNDARY functional (Stokes),
     (1/2)*oint(x dy - y dx), hence additive at a self-intersection; so
     W(C) = exp(-lambda*A(C)/2) is multiplicative, W(fig8) = W(sq1)*W(sq2), and the MM
     self-intersection factorization closes.  (Singer's 2D master field realizes this.)

 (B) D=4 [NG, demonstrated]: the area of a SPANNING surface of a 1D boundary is NOT
     boundary-determined (unlike the 2D signed area): two spanning surfaces of the same
     unit-square boundary -- flat (area 1) and a tent (area 2.236 at height 1) -- have
     different areas.  So "A(C)" is not a functional of C in D=4; W(C) is a SUM over
     spanning surfaces (a string functional), not exp(-sigma*A).  A single scalar sigma
     cannot be the closed solution.

 (C) structural [T]/[O]: the MM RHS is W(C1)*W(C2) -- a PRODUCT of the SAME functional --
     so the loop equation is closed in W and has NO n-point tower (contrast: DSE is an
     infinite tower).  But W(C) lives on infinite-dimensional loop space, the D=4 master
     field is not constructed, and the large-loop decay of W gives sigma (1-form /
     confinement), NOT m_0++ (0-form); m_0++ needs linearizing around the UNKNOWN master
     field.  => structurally better than DSE, but not a closed scalar gap equation.

Run:  py -3.12 research/unconditional/w7/pnt_largen_loop_witness.py
"""
import numpy as np


def signed_area_2d(pts):
    p = np.asarray(pts, float)
    x, y = p[:, 0], p[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def tri_area(a, b, c):
    a, b, c = np.array(a, float), np.array(b, float), np.array(c, float)
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a))


def main():
    print("=" * 74)
    print("PNT-LARGEN-LOOP witness")
    print("=" * 74)

    # ---------------- (A) D=2 exact closure ----------------
    print("\n(A) D=2 closure  [T exact]")
    sq1 = [(0, 0), (1, 0), (1, 1), (0, 1)]          # unit square, CCW
    sq2 = [(0, 0), (-1, 0), (-1, -1), (0, -1)]      # unit square, CCW, shares origin
    fig8 = sq1 + sq2 + [sq1[0]]                      # single closed figure-eight curve
    A1, A2 = signed_area_2d(sq1), signed_area_2d(sq2)
    A_fig = signed_area_2d(fig8)
    lam = 0.7
    W = lambda A: np.exp(-lam * A / 2.0)
    print(f"  A(sq1)={A1:.3f}  A(sq2)={A2:.3f}  A(fig8)={A_fig:.3f}  A1+A2={A1 + A2:.3f}")
    print(f"  area additive at self-intersection : {np.isclose(A_fig, A1 + A2)}")
    print(f"  W(fig8)={W(A_fig):.6f}   W(sq1)*W(sq2)={W(A1) * W(A2):.6f}")
    print(f"  W multiplicative (MM factorizes)   : {np.isclose(W(A_fig), W(A1) * W(A2))}")
    print("  => D=2: signed area is boundary-determined (Stokes) => exp(-lam A/2) solves MM.")

    # ---------------- (B) D=4 surface non-uniqueness ----------------
    print("\n(B) D=4 obstruction  [NG, demonstrated]")
    corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]   # unit-square boundary, z=0
    flat_area = 1.0
    for h in (0.5, 1.0, 2.0):
        apex = (0.5, 0.5, h)
        tent = sum(tri_area(corners[i], corners[(i + 1) % 4], apex) for i in range(4))
        print(f"  tent (apex height h={h:>3}): area = {tent:.4f}   vs flat = {flat_area:.4f}")
    print("  same boundary, many spanning surfaces with DIFFERENT area")
    print("  => in D=4 'A(C)' is NOT a boundary functional (contra D=2 signed area);")
    print("     W(C) = sum over spanning surfaces (string functional), not exp(-sigma*A).")

    # ---------------- (C) structural ----------------
    print("\n(C) structural  [T]/[O]")
    print("  MM RHS = W(C1)*W(C2): a PRODUCT of the same functional => closed in W, NO n-point tower.")
    print("  Contrast DSE: infinite tower of n-point functions.  => large-N loop eq structurally better.")
    print("  BUT: W(C) on infinite-dim loop space; D=4 master field not constructed (Singer: 2D only).")
    print("  Loop eq large-loop decay -> sigma (1-form/confinement); m_0++ (0-form) needs")
    print("  linearizing around the UNKNOWN master field.  => no closed scalar gap equation.")


if __name__ == "__main__":
    main()
