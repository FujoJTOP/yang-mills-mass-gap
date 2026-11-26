"""
UR1-fast — n=6 layer only, convention (A):  B_6 = sum_{Gamma} |mu(G_Gamma)|
Fast paths: tree (m=n-1) and unicyclic (m=n, mu=(-1)^(n-1)(L-1), L=cycle length);
DP only for cyclomatic number >= 2.  Reuses the validated lattice code from ur1_enum.py.
"""
import time, sys, json
import ur1_enum as E

n_target = int(sys.argv[1]) if len(sys.argv) > 1 else 6
t0 = time.time()
levels = E.animals(n_target)
S_list = levels[n_target - 1]
n = n_target
print(f"enumeration {n}: {len(S_list)} animals in {time.time()-t0:.1f}s", flush=True)

sgn = (-1) ** (n - 1)
t1 = time.time()
B = 0
n_tree = 0
n_uni = 0
n_dp = 0
maxabs = 0
hist = {}
step = max(1, len(S_list) // 10)
for idx, S in enumerate(S_list):
    # --- adjacency by link grouping ---
    groups = {}
    for i, p in enumerate(S):
        x, (mu, nu) = E.dec(p)
        for lk in ((E.sc(x) << 2) | mu,
                   (E.sc(E.add(x, E.DIRS[mu])) << 2) | nu,
                   (E.sc(E.add(x, E.DIRS[nu])) << 2) | mu,
                   (E.sc(x) << 2) | nu):
            g = groups.get(lk)
            if g is None:
                groups[lk] = [i]
            else:
                g.append(i)
    adj = [0] * n
    deg = [0] * n
    mE = 0
    for g in groups.values():
        for a in range(len(g)):
            for b in range(a + 1, len(g)):
                i, j = g[a], g[b]
                if not (adj[i] >> j) & 1:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
                    deg[i] += 1
                    deg[j] += 1
                    mE += 1
    if mE == n - 1:
        mu = sgn
        n_tree += 1
    elif mE == n:
        # peel degree-1 vertices -> remaining cycle, length L
        d = deg[:]
        alive = [True] * n
        stack = [i for i in range(n) if d[i] == 1]
        removed = 0
        while stack:
            i = stack.pop()
            if not alive[i] or d[i] != 1:
                continue
            alive[i] = False
            removed += 1
            for j in range(n):
                if alive[j] and (adj[i] >> j) & 1:
                    d[j] -= 1
                    if d[j] == 1:
                        stack.append(j)
        L = n - removed
        mu = sgn * (L - 1)
        n_uni += 1
    else:
        mu, _ = E.conn_sum(adj, n, 1)
        n_dp += 1
    B += abs(mu)
    maxabs = max(maxabs, abs(mu))
    hist[mu] = hist.get(mu, 0) + 1
    if idx % step == 0:
        print(f"  {idx}/{len(S_list)}  B={B}  [{time.time()-t1:.0f}s]", flush=True)

print(f"\nn={n}: count={len(S_list)} trees={n_tree} unicyclic={n_uni} DP={n_dp} "
      f"max|mu|={maxabs}")
print(f"B_{n} = {B}")
print(f"lambda_{n} = B^(1/n) = {B ** (1.0/n):.6f}   u*lambda = {E.U * B ** (1.0/n):.6f}")
print(f"m_{n}^(A) at u={E.U} = {E.U**n * B:.6e}")
print(f"m^{{1/n}} = {(E.U**n * B) ** (1.0/n):.6f}   ubar = {E.EDELTA*E.U:.10f}")
print(f"sign histogram (top 12): {sorted(hist.items(), key=lambda kv: -kv[1])[:12]}")
print(f"total {time.time()-t0:.1f}s")
json.dump(dict(n=n, count=len(S_list), trees=n_tree, unicyclic=n_uni, dp=n_dp,
               B=B, lam=B ** (1.0 / n), u=E.U, max_abs_mu=maxabs,
               hist={str(k): v for k, v in sorted(hist.items())}),
          open(f"ur1_enum_n{n}.json", "w"), indent=1)
