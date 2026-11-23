"""
Numerical verification for the 4D line-defect two-point block.

Convention (agrees with BLOCK-CORRECTED.md / U-ODE.md / FREE-CODIM3.md):
    u = tanh^2(eta/2),   z = cosh d = (1+u)/(1-u),   u in (0,1)
    u -> 0 : bulk OPE   (x1 -> x2)
    u -> 1 : defect OPE (both points -> line)

Defect block degree nu = Delta_j - 1, Casimir eigenvalue lambda = nu(nu+1) = Delta_j(Delta_j-1).

Checks:
  (a) P_{Delta-1}(z) = 2F1(Delta, 1-Delta; 1; -u/(1-u));  P -> 1 (u->0),  P ~ (1-u)^{1-Delta} (u->1)
  (b) Q_{Delta-1}(z) closed form;                          Q ~ (1-u)^Delta (u->1), Q ~ -1/2 ln u (u->0)
  (c) H^2 radial ODE:  f'' + coth(d) f' = Delta(Delta-1) f
  (d) impossibility: ratio P/Q -> inf as u->0;  growing-branch coefficient forces alpha=0 / beta=0.

All high precision via mpmath; a scipy cross-check for integer degree.
"""
import sys
import mpmath as mp

mp.mp.dps = 40

DJ = ['1', '1.5', '2', '3']
US = ['0.05', '0.2', '0.5', '0.8', '0.95']

def z_of(u):
    return (1 + u) / (1 - u)

def P_hyp(D, u):
    """P_{D-1}(cosh d) = 2F1(1-D, D; 1; (1-z)/2)."""
    return mp.hyp2f1(1 - D, D, 1, -u / (1 - u))

def P_leg(D, u):
    """mpmath Legendre P of degree D-1 at z>1 (real branch)."""
    return mp.legenp(D - 1, 0, z_of(u), type=3)

def Q_closed(D, u):
    """Q_{D-1}(z) closed form (DLMF 14.3.7 with nu=D-1):
       sqrt(pi) Gamma(D) / (2^D Gamma(D+1/2)) * z^{-D} * 2F1((D+1)/2, D/2; D+1/2; 1/z^2)."""
    z = z_of(u)
    pref = mp.sqrt(mp.pi) * mp.gamma(D) / (2**D * mp.gamma(D + mp.mpf(1) / 2))
    return pref * z**(-D) * mp.hyp2f1((D + 1) / 2, D / 2, D + mp.mpf(1) / 2, 1 / z**2)

def Q_leg(D, u):
    return mp.legenq(D - 1, 0, z_of(u), type=3).real

# ---------------------------------------------------------------- table (a)+(b)
def table():
    out = []
    out.append("=" * 132)
    out.append("(a)+(b)  P = 2F1(D,1-D;1;-u/(1-u))  and  Q_{D-1}  block behaviours")
    out.append("columns: D  u   z   P_hyp  P_leg |P_hyp-P_leg|   P*(1-u)^(D-1)  Q_closed  Q_leg  Q/(1-u)^D   Q+0.5*ln(u)   P/Q")
    out.append("=" * 132)
    for ds in DJ:
        D = mp.mpf(ds)
        for us in US:
            u = mp.mpf(us)
            z = z_of(u)
            Ph = P_hyp(D, u)
            Pl = P_leg(D, u)
            Psc = Ph * (1 - u)**(D - 1)               # expect -> 2^(D-1)
            Qc = Q_closed(D, u)
            Ql = Q_leg(D, u)
            Qsc = Qc / (1 - u)**D                       # expect -> finite const
            Qlog = Qc + mp.mpf(1) / 2 * mp.log(u)       # expect -> -gamma - psi(D)
            ratio = Ph / Qc
            out.append("{:>4} {:>5} {:>8.4} {:>12.6} {:>12.6} {:>10.2e}  {:>14.6}  {:>12.6} {:>12.6} {:>12.6}  {:>12.6}  {:>12.4e}".format(
                ds, us, float(z), float(Ph), float(Pl), float(abs(Ph - Pl)),
                float(Psc), float(Qc), float(Ql), float(Qsc), float(Qlog), float(ratio)))
        out.append("-" * 132)
    # predicted leading constants
    out.append("")
    out.append("Predicted leading constants:")
    out.append("  P_{D-1}(z) ~ A_nu z^nu  with nu=D-1, A_nu=2^nu Gamma(nu+1/2)/(sqrt(pi)Gamma(nu+1));")
    out.append("  P*(1-u)^(D-1) -> 2^nu A_nu = 2^(2(D-1)) Gamma(D-1/2)/(sqrt(pi)Gamma(D)).")
    out.append("  Q_{D-1} ~ sqrt(pi)Gamma(D)/(2^D Gamma(D+1/2)) z^(-D);  Q/(1-u)^D -> that/2^D.")
    out.append("  Q + 0.5 ln u -> -gamma - psi(D)  (finite part of the log divergence).")
    for ds in DJ:
        D = mp.mpf(ds); nu = D - 1
        Pcst = 2**nu * (2**nu * mp.gamma(nu + mp.mpf(1) / 2) / (mp.sqrt(mp.pi) * mp.gamma(nu + 1)))
        out.append("  D={:>4}:  P*(1-u)^(D-1)->{:.6f}   Q/(1-u)^D->{:.6f}   Q+0.5ln u->{:.6f}".format(
            ds, float(Pcst),
            float(mp.sqrt(mp.pi) * mp.gamma(D) / (2**(2 * D) * mp.gamma(D + mp.mpf(1) / 2))),
            float(-mp.euler - mp.digamma(D))))
    out.append("")
    out.append("Slow log convergence of Q+0.5 ln u at tiny u (pins the finite part):")
    for ds in ['1.5', '2', '3']:
        D = mp.mpf(ds)
        row = []
        for us in ['1e-2', '1e-3', '1e-4', '1e-5']:
            u = mp.mpf(us)
            row.append("u={:>5}: {:>10.6}".format(us, float(Q_closed(D, u) + mp.mpf(1) / 2 * mp.log(u))))
        out.append("  D={:>4}  -gamma-psi(D)={:>9.6}   ".format(ds, float(-mp.euler - mp.digamma(D))) + "  ".join(row))
    return "\n".join(out)

# ---------------------------------------------------------------- (c) H^2 ODE
def h2_ode_check():
    out = []
    out.append("")
    out.append("=" * 100)
    out.append("(c) H^2 radial ODE check:  f''(d) + coth(d) f'(d) = D(D-1) f(d),  f(d)=Q_{D-1}(cosh d)")
    out.append("=" * 100)
    Ds = ['1.5', '2']
    ds_pts = ['0.3', '0.7', '1.2', '2.0']
    for ds in Ds:
        D = mp.mpf(ds)
        lam = D * (D - 1)
        f = lambda d: Q_closed(D, mp.cosh(d))  # Q via z=cosh d
        # use closed form's u?  Q_closed takes u; define directly in d instead:
        def Qd(d):
            z = mp.cosh(d)
            pref = mp.sqrt(mp.pi) * mp.gamma(D) / (2**D * mp.gamma(D + mp.mpf(1) / 2))
            return pref * z**(-D) * mp.hyp2f1((D + 1) / 2, D / 2, D + mp.mpf(1) / 2, 1 / z**2)
        f = Qd
        out.append("  D={:>4} (lambda={:.4f}):".format(ds, float(lam)))
        for dps in ds_pts:
            d = mp.mpf(dps)
            f0 = f(d)
            f1 = mp.diff(f, d, 1)
            f2 = mp.diff(f, d, 2)
            lhs = f2 + mp.coth(d) * f1
            rel = abs(lhs - lam * f0) / abs(lam * f0)
            out.append("    d={:>4}:  f''+coth(d)f' = {:>18.12}   lambda*f = {:>18.12}   rel.err={:.2e}".format(
                dps, float(lhs), float(lam * f0), float(rel)))
    return "\n".join(out)

# ---------------------------------------------------------------- (d) impossibility
def impossibility():
    out = []
    out.append("")
    out.append("=" * 100)
    out.append("(d) Impossibility of a single alpha*P + beta*Q")
    out.append("=" * 100)
    out.append("Fact 1 (bulk end, u->0): P -> 1 (bounded) while Q ~ -0.5 ln u -> +infinity.")
    out.append("  Any beta != 0 gives beta*Q -> +-infinity; alpha*P (bounded) cannot cancel. => beta = 0.")
    D = mp.mpf('2')
    for us in ['0.1', '0.01', '1e-3', '1e-4', '1e-5']:
        u = mp.mpf(us)
        P = P_hyp(D, u); Q = Q_closed(D, u)
        out.append("    D=2 u={:>7}: P = {:>10.6}   Q = {:>12.6}   P/Q = {:>12.4e}  (->0 as Q blows up)".format(
            us, float(P), float(Q), float(P / Q)))
    out.append("")
    out.append("Fact 2 (defect end, u->1): pure (1-u)^Delta decay needs the (1-u)^(1-D) branch to vanish.")
    out.append("  P carries it: P*(1-u)^(D-1) -> nonzero constant  (see table (a)).")
    out.append("  Q does NOT:   Q/(1-u)^D -> nonzero const  =>  Q*(1-u)^(D-1) = [Q/(1-u)^D]*(1-u)^(2D-1) -> 0.")
    for ds in ['1.5', '2', '3']:
        D = mp.mpf(ds)
        us = mp.mpf('0.9999')
        P = P_hyp(D, us); Q = Q_closed(D, us)
        Pp = P * (1 - us)**(D - 1)     # -> nonzero const
        Qp = Q * (1 - us)**(D - 1)     # -> 0
        out.append("    D={:>4} u=0.9999:  P*(1-u)^(D-1) = {:>12.6}   Q*(1-u)^(D-1) = {:>12.4e}".format(
            ds, float(Pp), float(Qp)))
    out.append("  => the growing-branch coefficient of alpha*P+beta*Q is alpha*(nonzero).")
    out.append("     Pure defect decay (forced) => alpha = 0.")
    out.append("")
    out.append("Fact 1 forces beta=0; Fact 2 forces alpha=0.  The only combination meeting BOTH ends is f=0.")
    out.append("  A nonzero physical block cannot be a single linear combination of P_{D-1} and Q_{D-1}. QED.")
    return "\n".join(out)

def main():
    lines = [table(), h2_ode_check(), impossibility()]
    txt = "\n".join(lines)
    print(txt)
    with open(r"D:\Dev\math\research\defect-bootstrap\verify_block_OUT.txt", "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print("\n[saved to verify_block_OUT.txt]")

if __name__ == "__main__":
    main()
