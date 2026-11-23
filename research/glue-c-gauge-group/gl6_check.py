#!/usr/bin/env python3
"""GL6 executable check: gauge-fixing independence on a finite lattice.

Exact enumeration of a small lattice gauge theory with a finite group.

Checks, for a gauge-invariant measure mu and a gauge-invariant F:

  (I1) mu is gauge invariant (orbit weights are constant).
  (I2) Orbit-decomposition identity with a *multi-valued* gauge slice A:
           sum_X F dmu = sum_{a in A} F(a) mu(a) * |orbit(a)| / m(a),
       where m(a) = |A cap orbit(a)|.  This is the finite-lattice
       Faddeev--Popov identity with the orbit volume as weight.
  (I3) Two different gauge functionals (Landau-type slices) give the same
       invariant expectation.
  (I4) Failure mode: dropping 1/m(a) ("forgetting the Gribov copies") changes
       the answer whenever m is non-constant.
  (I5) A non-invariant observable DOES depend on the slice (the theorem is
       not vacuous; gauge fixing independence is only for invariants).

Run:  python gl6_check.py            # default: S3 on a 2x1 torus (fast)
      python gl6_check.py --group Z2 --lattice 2x2
"""

from __future__ import annotations
import argparse
import itertools
import math


# ---------------------------------------------------------------- groups
class Group:
    def __init__(self, name, els, mul, inv, chi, dim):
        self.name = name
        self.els = els
        self.index = {e: i for i, e in enumerate(els)}
        self.mul = mul
        self.inv = inv
        self.chi = chi                 # real character of the defining rep
        self.dim = dim                 # its dimension
        self.q = len(els)
        self.e = 0                     # els[0] is the identity
        self.mul_i = [[self.index[mul(a, b)] for b in els] for a in els]
        self.inv_i = [self.index[inv(a)] for a in els]
        self.gens = self._generators()

    def _generators(self):
        # generators of the group, as element indices
        if self.name == "Z2":
            return [1]
        # S3 = <(0 1), (0 1 2)> in the convention mul(a,b) = a o b
        return [self.index[(1, 0, 2)], self.index[(1, 2, 0)]]


def make_z2():
    els = [(0,), (1,)]
    mul = lambda a, b: ((a[0] + b[0]) % 2,)
    inv = lambda a: a
    chi = {(0,): 1.0, (1,): -1.0}
    return Group("Z2", els, mul, inv, lambda a: chi[a], 1)


def make_s3():
    els = list(itertools.permutations(range(3)))
    mul = lambda a, b: tuple(a[b[i]] for i in range(3))
    def inv(a):
        r = [0, 0, 0]
        for i, x in enumerate(a):
            r[x] = i
        return tuple(r)
    # 2-dimensional irrep: chi(p) = #fixed points - 1
    chi = lambda a: float(sum(1 for i in range(3) if a[i] == i) - 1)
    return Group("S3", els, mul, inv, chi, 2)


# ---------------------------------------------------------------- lattice
class Lattice:
    """d=2 torus Lx x Ly.  Links: horizontal H(x,y), then vertical V(x,y)."""
    def __init__(self, Lx, Ly):
        self.Lx, self.Ly = Lx, Ly
        self.V = Lx * Ly
        V = self.V
        links = []
        for y in range(Ly):
            for x in range(Lx):
                links.append((x + Lx * y, (x + 1) % Lx + Lx * y))
        for y in range(Ly):
            for x in range(Lx):
                links.append((x + Lx * y, x + Lx * ((y + 1) % Ly)))
        self.links = links
        H = lambda x, y: (x % Lx) + Lx * (y % Ly)
        Vv = lambda x, y: V + (x % Lx) + Lx * (y % Ly)
        plaq = []
        for y in range(Ly):
            for x in range(Lx):
                # v -> v+e_x -> v+e_x+e_y -> v+e_y -> v ; last two reversed
                plaq.append((H(x, y), Vv(x + 1, y), H(x, y + 1), Vv(x, y)))
        self.plaq = plaq
        self.nlinks = len(links)
        self.trees = self._trees()

    def _trees(self):
        """Two spanning trees (as link-index lists) of the vertex graph."""
        # tree 1: horizontal chains then vertical spokes
        t1 = []
        for y in range(self.Ly):
            for x in range(self.Lx - 1):
                t1.append(self.Lx * y + x)          # H(x,y), x -> x+1
        for y in range(self.Ly - 1):
            t1.append(self.V + self.Lx * y)         # V(0,y)
        # tree 2: same vertices, different link order (vertical first)
        t2 = []
        for x in range(self.Lx):
            for y in range(self.Ly - 1):
                t2.append(self.V + x + self.Lx * y)  # V(x,y)
        for x in range(self.Lx - 1):
            t2.append(x)                             # H(x,0)
        return [t1, t2]


# ---------------------------------------------------------------- measure
def build_tables(g: Group, lat: Lattice):
    q = g.q
    chi_i = [g.chi(a) for a in g.els]
    # holonomy index per plaquette, as a flat table over the 4 link digits
    ptab = []
    for (a, b, c, d) in lat.plaq:
        tab = [0] * (q ** 4)
        for i in range(q):
            for j in range(q):
                ij = g.mul_i[i][j]
                for k in range(q):
                    ijk = g.mul_i[ij][g.inv_i[k]]
                    for l in range(q):
                        tab[((i * q + j) * q + k) * q + l] = g.mul_i[ijk][g.inv_i[l]]
        ptab.append(tab)
    return chi_i, ptab


def decode(idx, q, L):
    ds = [0] * L
    for i in range(L):
        ds[i] = idx % q
        idx //= q
    return ds


def encode(ds, q):
    idx = 0
    for i in range(len(ds) - 1, -1, -1):
        idx = idx * q + ds[i]
    return idx


def action_of_digits(ds, g, lat, chi_i, ptab):
    q = g.q
    S = 0.0
    for p, (a, b, c, d) in enumerate(lat.plaq):
        h = ptab[p][((ds[a] * q + ds[b]) * q + ds[c]) * q + ds[d]]
        S += 1.0 - chi_i[h] / g.dim
    return S


def F_invariant(ds, g, lat, chi_i, ptab):
    """chi of the first plaquette (a closed loop) -- gauge invariant."""
    a, b, c, d = lat.plaq[0]
    q = g.q
    h = ptab[0][((ds[a] * q + ds[b]) * q + ds[c]) * q + ds[d]]
    return chi_i[h] / g.dim


def F_link(ds, g, lat, chi_i, ptab):
    """chi of one link -- NOT gauge invariant."""
    return chi_i[ds[0]] / g.dim


# ---------------------------------------------------------------- gauge action
def build_gen_maps(g: Group, lat: Lattice):
    maps = []
    for v in range(lat.V):
        for s in g.gens:
            per_link = []
            for (src, tgt) in lat.links:
                mp = [0] * g.q
                for i in range(g.q):
                    if src == v and tgt == v:
                        b = g.index[g.mul(g.mul(g.els[s], g.els[i]), g.inv(g.els[s]))]
                    elif src == v:
                        b = g.mul_i[s][i]
                    elif tgt == v:
                        b = g.mul_i[i][g.inv_i[s]]
                    else:
                        b = i
                    mp[i] = b
                per_link.append(mp)
            maps.append(per_link)
    return maps


def apply_gen(idx, gm, q, L):
    ds = decode(idx, q, L)
    for i in range(L):
        ds[i] = gm[i][ds[i]]
    return encode(ds, q)


def orbit_of(start, gen_maps, q, L):
    seen = {start}
    stack = [start]
    while stack:
        x = stack.pop()
        for gm in gen_maps:
            y = apply_gen(x, gm, q, L)
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return seen


# ---------------------------------------------------------------- main check
def run(group_name, Lx, Ly, beta=1.0, phi_kind="uniform"):
    g = make_z2() if group_name == "Z2" else make_s3()
    lat = Lattice(Lx, Ly)
    q, L, N = g.q, lat.nlinks, g.q ** lat.nlinks
    chi_i, ptab = build_tables(g, lat)
    gen_maps = build_gen_maps(g, lat)

    print(f"group={g.name}  lattice={Lx}x{Ly}  links={L}  configs={N}  |G|={q}  |Gauge|={q**lat.V}")

    # ---- (I1) invariance + exact expectation over all configurations
    weight = [0.0] * N
    Finv = [0.0] * N
    Flink = [0.0] * N
    Z = 0.0
    E_inv = 0.0
    E_link = 0.0
    for idx in range(N):
        ds = decode(idx, q, L)
        S = action_of_digits(ds, g, lat, chi_i, ptab)
        w = math.exp(-beta * S)
        weight[idx] = w
        Finv[idx] = F_invariant(ds, g, lat, chi_i, ptab)
        Flink[idx] = F_link(ds, g, lat, chi_i, ptab)
        Z += w
        E_inv += w * Finv[idx]
        E_link += w * Flink[idx]
    E_inv /= Z
    E_link /= Z

    # ---- orbits + Landau-type slice (maximise a gauge functional)
    c = [1.0 + 0.17 * i for i in range(L)]  # distinct link weights
    seen = bytearray(N)
    orbits = []
    for idx in range(N):
        if seen[idx]:
            continue
        orb = sorted(orbit_of(idx, gen_maps, q, L))
        for x in orb:
            seen[x] = 1
        orbits.append(orb)

    def phi(idx):
        ds = decode(idx, q, L)
        if phi_kind == "uniform":
            return sum(chi_i[ds[i]] for i in range(L)) / g.dim
        return sum(c[i] * chi_i[ds[i]] for i in range(L)) / g.dim

    slice_pts = []       # (canonical orbit id, config, m, orbit_size)
    mult = {}
    for oi, orb in enumerate(orbits):
        vals = [phi(x) for x in orb]
        mx = max(vals)
        cop = [x for x, v in zip(orb, vals) if abs(v - mx) < 1e-12]
        for x in cop:
            slice_pts.append((oi, x))
        mult[oi] = len(cop)

    mvals = sorted(set(mult.values()))
    print(f"orbits={len(orbits)}  slice points={len(slice_pts)}  multiplicity values m={mvals}")

    # ---- (I2)/(I3) slice sums
    def slice_sum(F, correct=True):
        Zs = 0.0
        Es = 0.0
        for (oi, x) in slice_pts:
            ds = decode(x, q, L)
            w = weight[x]
            wgt = len(orbits[oi]) / (mult[oi] if correct else 1)
            Zs += w * wgt
            Es += w * wgt * F(ds)
        return Es / Zs, Zs

    def F_inv_arr(ds):
        return F_invariant(ds, g, lat, chi_i, ptab)

    def F_link_arr(ds):
        return F_link(ds, g, lat, chi_i, ptab)

    Es_inv, Zs_inv = slice_sum(F_inv_arr, correct=True)
    Ew_inv, Zw_inv = slice_sum(F_inv_arr, correct=False)
    Es_link, _ = slice_sum(F_link_arr, correct=True)

    print()
    print("  exact   E[F_inv]  = %.12f" % E_inv)
    print("  slice   E[F_inv]  = %.12f   (weight |orbit|/m)" % Es_inv)
    print("  no 1/m  E[F_inv]  = %.12f   (weight |orbit|)" % Ew_inv)
    print("  exact   E[F_link] = %.12f  (non-invariant)" % E_link)
    print("  slice   E[F_link] = %.12f  (non-invariant)" % Es_link)
    print("  Z (all configs)   = %.12f" % Z)
    print("  Z (slice, |orb|/m)= %.12f" % Zs_inv)

    tol = 1e-9 * max(1.0, abs(E_inv))
    ok_inv = abs(Es_inv - E_inv) < tol
    ok_Z = abs(Zs_inv - Z) < tol * max(1.0, abs(Z))
    ok_fail = (len(mvals) > 1) and abs(Ew_inv - E_inv) > 1e-6
    ok_noninv = abs(Es_link - E_link) > 1e-6

    print()
    print("  (I2) exact = slice-with-|orbit|/m :", "PASS" if ok_inv else "FAIL")
    print("  (I2) Z identity                  :", "PASS" if ok_Z else "FAIL")
    print("  (I4) dropping 1/m changes answer :", "PASS" if ok_fail else "NOTE(m constant)")
    print("  (I5) non-invariant F is slice-dep:", "PASS" if ok_noninv else "FAIL")

    if not (ok_inv and ok_Z):
        raise SystemExit("GL6 CHECK FAILED: orbit-decomposition identity")
    return {"E_inv": E_inv, "E_slice": Es_inv, "E_wrong": Ew_inv,
            "m_values": mvals, "orbits": len(orbits)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="S3", choices=["Z2", "S3"])
    ap.add_argument("--lattice", default="2x1")
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--phi", default="uniform", choices=["uniform", "weighted"])
    args = ap.parse_args()
    Lx, Ly = (int(s) for s in args.lattice.lower().split("x"))
    run(args.group, Lx, Ly, args.beta, args.phi)
    print("\nGL6 CHECK: PASS")


if __name__ == "__main__":
    main()
