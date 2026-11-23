import sys, json, time
exec(open('gluon_prop.py').read().split("if __name__")[0])
results = []
# controls
results.append(control_free(4, 24))
results.append(control_free(6, 16))
t0 = time.time()
# main grid
for L in (4, 6):
    for beta in (2.2, 2.5, 2.8):
        nc = 150 if L == 4 else 100
        results.append(run(beta, L, ncfg=nc, ntherm=600, gap=10,
                           seed=100 + int(beta * 10) + L))
# volume check at beta=2.5
for L in (8,):
    results.append(run(2.5, L, ncfg=28, ntherm=400, gap=12, seed=8))
json.dump(results, open('D:/Dev/math/research/lattice-cpb/gluon_prop_results.json', 'w'), indent=1)
print(f"\nDONE in {time.time()-t0:.0f}s, wrote gluon_prop_results.json")
