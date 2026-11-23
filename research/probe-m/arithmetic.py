"""Reproducible arithmetic for research/probe-m/REPORT.md.

Checks, from handoff-proof REPORT 1.5 / RB2:
  * two-loop t_up(labm), q0 = (X/Lambda) exp(-t_up)
  * the UV floor xi/a >= 1/q1 (cutoff 1/a) or 1/(pi q1) (cutoff pi/a)
  * the IR floor L >= 4 xi/a  (project's L/xi >= 4 evidence standard)
  * cost model N_eff = 4000 (xi/0.93)^2 (8.2/percent)^2  (lattice-cpb 10.3)
"""
import math

C = 24 * math.pi**2 / 11          # C_RG
A_2L = 51 / 121


def t_up(labm):
    return C / labm - A_2L * math.log(1 + 121 * C / (51 * labm))


def q0(X_over_Lambda, labm=6.0):
    return X_over_Lambda * math.exp(-t_up(labm))


def n_eff(xi, pct):
    return 4000 * (xi / 0.93) ** 2 * (8.2 / pct) ** 2


def main():
    t = t_up(6.0)
    print(f"t_up^2L(6) = {t:.4f}")
    for X, name in ((8.0, "m/Lambda"), (2.0, "sqrt(sigma)/Lambda")):
        q = q0(X)
        print(f"  {name} = {X}: q0 = {q:.4f}  -> 1/q0 = {1/q:.3f}")

    print("\nfloors and requirements (q1 = q0 at threshold):")
    for name, q1, xi_at_b6 in (("m", q0(8.0), 1.0), ("sigma", q0(2.0), 4.4)):
        print(f"  {name}: q1={q1:.4f}")
        print(f"    UV floor, cutoff 1/a : xi/a >= {1/q1:.2f}")
        print(f"    UV floor, cutoff pi/a: xi/a >= {1/(math.pi*q1):.2f}")
        print(f"    IR floor L >= 4*xi/a at that xi:")
        for xi in (1/q1, 1/(math.pi*q1)):
            print(f"       xi/a={xi:.2f} -> L >= {4*xi:.1f}")

    print("\ncost model (L=12**3x24 for m, 32**4 for sigma; 50us/site NumPy):")
    for name, xi, vol, L in (("m", 1.75, 12**3 * 24, 12), ("sigma", 7.0, 32**4, 32)):
        N = n_eff(xi, 5.0)
        sweeps = 20 * N          # 2*tau_int, tau=10
        sec = sweeps * vol * 50e-6
        print(f"  {name}: xi/a={xi}, V={vol}, N_eff={N:.3g}, sweeps={sweeps:.3g}")
        print(f"       NumPy ~ {sec/3600:.0f} h ; @7x ~ {sec/3600/7:.0f} h "
              f"; @50x ~ {sec/3600/50:.0f} h")

    # quenched SU(3) Wilson, r0/a = 5.368, 7.380, 9.660 (Necco-Sommer),
    # r0 = 0.5 fm -> a^-1 = 2.12, 2.91, 3.81 GeV; m(0++) r0 = 4.33 -> 1.71 GeV
    print("\nxi/a(0++) from r0/a and m(0++)=1.71 GeV:")
    ainv = {6.0: 2.12, 6.2: 2.91, 6.4: 3.81}
    for beta, a in ainv.items():
        print(f"  beta={beta}: am={1.71/a:.3f}  xi/a={a/1.71:.2f}  "
              f"xi_sigma/a={a/0.44:.2f}")
    target = 1.75
    b = 6.0 + (target - ainv[6.0] / 1.71) / ((ainv[6.2] - ainv[6.0]) / 1.71) * 0.2
    print(f"  xi/a = {target} at beta ~ {b:.2f}")

    print("\nUV floor is certificate-independent: a <= q1/X = exp(-t_up)/Lambda")
    print(f"  exp(-t_up)/Lambda = {math.exp(-t_up(6.0)):.4f}/Lambda ; "
          f"for Lambda=0.2 GeV: a <= {math.exp(-t_up(6.0))/0.2:.3f} GeV^-1 "
          f"-> a^-1 >= {0.2/math.exp(-t_up(6.0)):.2f} GeV (beta ~ 6.2)")


if __name__ == "__main__":
    main()
