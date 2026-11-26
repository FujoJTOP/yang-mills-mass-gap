import io, re

p = "r3_l4_volume.py"
s = io.open(p, encoding="utf-8").read().replace("\r\n", "\n")

old_print = s[s.index('    r = rows[L]\n    say(f"  L={L}'):s.index('dump("S2", rows)')]
new_print = '''    r = rows[L]
    npq = 6 * L ** 4
    say(f"  L={L}  E={E:5d}  plaquettes={npq:5d}   [{(time.time() - t0):.1f}s, n={nmeas}, acc={acc:.2f}]")
    say(f"     A = <|grad S_W|^2>_mu      = {A:11.3f} +- {Aerr:7.3f}   per link {A / E:7.3f}"
        f"   per plaque {A / npq:7.3f}")
    say(f"     B = sup_phi lower bound    = {B:11.3f}             per link {B / E:7.3f}"
        f"   per plaque {B / npq:7.3f}   (B/A = {B / A:.3f})")
    say(f"     C = Haar random field      = {C:11.3f} +- {Cerr:7.3f}   per link {C / E:7.3f}"
        f"   per plaque {C / npq:7.3f}")
    say(f"     per plaquette S_W (mu)     = {r['S_mean']:.4f}   (Haar/infinite-T value: "
        f"{S_W(haar(L)) / (6 * L ** 4):.4f})")
    say(f"     tr Cov(grad S)  = {trC_modes:11.3f} (modes)  {trC_direct:11.3f} (direct)   per link {trC_direct / E:8.4f}")
    say(f"     lambda_max Cov  = {lmax:+11.4f}    lambda(k=0) = {lok:+11.4f}"
        f"    median mode = {r['lam_med']:+.4f}")
'''
s = s.replace(old_print, new_print)

tail_start = s.index('def fit(ys, Ls=(2, 4, 6)):')
new_tail = '''def fit(ys):
    Ls = sorted(ys)
    y = np.log([ys[L] for L in Ls])
    x = np.log(Ls)
    sl = np.polyfit(x, y, 1)[0]
    loc = [float(np.log(ys[b] / ys[a]) / np.log(b / a)) for a, b in zip(Ls[:-1], Ls[1:])]
    return float(sl), loc


say("-" * 108)
say("S2b  L-exponent fits (log-log through L=2,4,6,8; local slopes between neighbours)")
for name, key in (("A  <|grad S|^2>_mu", "A"), ("B  sup_phi lower bound", "B"),
                  ("C  Haar random field", "C"), ("tr Cov(grad S_W)", "trC_direct"),
                  ("lambda_max Cov", "lam_max"), ("lambda(k=0) Cov", "lam_k0")):
    ys = {L: rows[L][key] for L in rows}
    sl, loc = fit(ys)
    say(f"    {name:24s} : exponent {sl:6.3f}   (local " + " / ".join(f"{v:6.3f}" for v in loc) + ")")
dump("S2fits", {k: fit({L: rows[L][k] for L in rows}) for k in
                ("A", "B", "C", "trC_direct", "lam_max", "lam_k0")})
say("=" * 108)
LOG.close()
print("\\nwrote r3_l4.json / r3_l4.log")
'''
s = s[:tail_start] + new_tail
io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("patched", len(s))
