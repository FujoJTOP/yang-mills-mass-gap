#!/usr/bin/env python3
r"""
cov_2dym.py -- COV-T, model (3): the 2D YM character tower, and the FLOW-CHOICE test.

Same single-link SU(2)=S^3 machinery as cov_link.py, but the flow kernel K_t is a free input:

  FLOW A (the BBD/BD choice, = cov_link.py):  K_t = sum_l (l+1) e^{-l(l+2) t} chi_l
                                              = the heat semigroup of the group.  lam_1 = 3.
  FLOW B (the 2D YM transfer matrix):         K_t = sum_l (l+1) e^{-E_l t} chi_l,
        E_l = -log(c_{l/2}/c_0),   c_j(beta) = 2 I_{2j+1}(2 beta)/(2 beta)   (r1_2dym.py convention)
        => E_1 = -log(c_{1/2}/c_0) = the *physical* 2D Yang-Mills mass gap.

Everything else is identical:
    B_t(x)     = e^{-beta} sum_l c_l e^{-E_l t} chi_l(psi_x)      (class-function convolution)
    mu_{t,x}   propto K_t(x^-1 w) e^{-S_W(w)} dw
    Cov        = (beta^2/4) Cov_mu(q_vec)
    Hess V_t   = (beta/4) E_mu[q_0] 1_3 - (beta^2/4) Cov_mu(q_vec)
    q(t)       = ||(Hess V_t)_-||_op

The point of the model: the DECAY RATE of q(t) is the smallest nonzero E_l, i.e. it is a property
of the chosen flow, not of the theory.  Flow A: rate = 3 = one-link group gap (beta-independent).
Flow B: rate = E_1 = 2D YM mass gap (beta-dependent)  =>  the "integrability of C(t) needs m>0" step
is circular *for that flow*, and only for that flow.

Run: py -3.12 -I cov_2dym.py    numpy/scipy only, deterministic, single process.
Writes cov_2dym.log.
"""
import numpy as np
from scipy.special import iv

LOG = []
LOGF = "cov_2dym.log"


def say(s=""):
    print(s)
    LOG.append(s)


def flush():
    with open(LOGF, "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


def chi(l, psi):
    s = np.sin(psi)
    out = np.empty_like(np.asarray(psi, dtype=float))
    small = np.abs(s) < 1e-9
    out[~small] = np.sin((l + 1) * psi[~small]) / s[~small]
    out[small] = (l + 1) * np.where(np.cos(psi[small]) > 0, 1.0, (-1.0) ** l)
    return out


def c_l(l, beta):
    return iv(l, beta) - iv(l + 2, beta)


def c_2dym(j, beta):
    """2D YM character coefficient (r1_2dym.py): c_j(beta) = 2 I_{2j+1}(2 beta)/(2 beta)."""
    return 2.0 * iv(2 * j + 1, 2.0 * beta) / (2.0 * beta)


def s3_grid(npsi=80, nth=80, nphi=80):
    xg, wg = np.polynomial.legendre.leggauss(npsi)
    psi = 0.5 * np.pi * (xg + 1.0)
    wpsi = 0.5 * np.pi * wg
    xt, wt = np.polynomial.legendre.leggauss(nth)
    theta = np.arccos(xt)
    phi = 2.0 * np.pi * np.arange(nphi) / nphi
    PS, TH, PH = np.meshgrid(psi, theta, phi, indexing="ij")
    W = (wpsi[:, None, None] * wt[None, :, None] * (2.0 * np.pi / nphi)
         * np.sin(PS) ** 2 * np.ones_like(PH)) / (2.0 * np.pi ** 2)
    q0 = np.cos(PS).ravel()
    st = np.sin(PS).ravel()
    qv = np.stack([st * np.sin(TH).ravel() * np.cos(PH).ravel(),
                   st * np.sin(TH).ravel() * np.sin(PH).ravel(),
                   st * np.cos(TH).ravel()], axis=1)
    return q0, qv, W.ravel()


class Tower:
    """single SU(2) link;  flow kernel K_t = sum_l (l+1) e^{-E_l t} chi_l."""

    def __init__(self, beta, E, lmax=240, ndelta=8192, npsi=80, nth=80, nphi=80):
        self.beta = beta
        self.E = E
        self.lmax = lmax
        self.nd = ndelta
        self.q0, self.qv, self.w = s3_grid(npsi, nth, nphi)
        self.c = np.array([c_l(l, beta) for l in range(lmax + 1)])

    def B_char(self, psi_x, t):
        tot = 0.0
        ps = np.array([psi_x])
        for l in range(self.lmax + 1):
            tot += self.c[l] * np.exp(-self.E[l] * t) * float(chi(l, ps)[0])
        return np.exp(-self.beta) * tot

    def kernel_lut(self, t):
        d = np.linspace(0.0, np.pi, self.nd)
        tot = np.zeros_like(d)
        for l in range(self.lmax + 1):
            co = (l + 1) * np.exp(-self.E[l] * t)
            if co < 1e-18 and l > 0:
                break
            tot += co * chi(l, d)
        return d, tot

    def moments(self, psi_x, t):
        beta = self.beta
        cx, sx = np.cos(psi_x), np.sin(psi_x)
        ip = cx * self.q0 + sx * self.qv[:, 0]
        delta = np.arccos(np.clip(ip, -1.0, 1.0))
        d, k = self.kernel_lut(t)
        Kt = np.interp(delta, d, k)
        rho = self.w * Kt * np.exp(beta * self.q0)
        Z = rho.sum() * np.exp(-beta)
        wn = rho / rho.sum()
        return Z, float(wn @ self.q0), wn @ self.qv, (self.qv * wn[:, None]).T @ self.qv

    def sup_x(self, t, nalpha=9):
        """sup over x of ||Cov||_op and sup over x of ||(Hess V_t)_-||_op (separately)."""
        Cs, qs, xs, lms, Zs = [], [], [], [], []
        for a in np.linspace(0.0, np.pi, nalpha):
            Z, m0, mv, Mv = self.moments(a, t)
            C = Mv - np.outer(mv, mv)
            Cs.append(np.linalg.norm((self.beta ** 2 / 4) * C, 2))
            HV = (self.beta / 4) * m0 * np.eye(3) - (self.beta ** 2 / 4) * C
            ev = np.linalg.eigvalsh(HV)
            qs.append(max(0.0, -ev.min()))
            xs.append(a)
            lms.append(ev.min())
            Zs.append(Z)
        i = int(np.argmax(Cs))
        j = int(np.argmax(qs))
        return dict(t=t, C=Cs[i], x=xs[i], q=qs[j], xq=xs[j], lmin=lms[j], Z=Zs[j])


def spectrum_A(lmax):
    return np.array([l * (l + 2.0) for l in range(lmax + 1)])


def spectrum_B(beta, lmax):
    c0 = c_2dym(0.0, beta)
    E = []
    for l in range(lmax + 1):
        r = c_2dym(0.5 * l, beta) / c0
        E.append(-np.log(r) if r > 1e-300 else 1e3)
    return np.array(E)


def curve(mod, taus, lmin=3.0):
    rows = []
    for tau in taus:
        t = tau / lmin
        r = mod.sup_x(t)
        rows.append(r)
    return rows


def report(tag, mod, taus, lmin, unit):
    say(f"  [{tag}]  spectrum: E_0 = {mod.E[0]:.3e}, E_1 = {mod.E[1]:.6f}, E_2 = {mod.E[2]:.6f}, "
        f"E_3 = {mod.E[3]:.6f}   (lam_1 used for the time axis = {lmin})")
    say(f"  {'lam1*t':>8} {'C_abs':>12} {'q(t)':>13} {'lam_min HV':>13} {'local rate':>11} {'B_quad/B_char':>14}")
    prev = None
    loc = []
    for r in curve(mod, taus, lmin):
        lr = np.nan
        if prev is not None and prev["q"] > 1e-13 and r["q"] > 1e-15:
            lr = -np.log(r["q"] / prev["q"]) / (r["t"] - prev["t"])
        loc.append(lr)
        Bc = mod.B_char(r["x"], r["t"])
        say(f"  {r['t']*lmin:8.4g} {r['C']:12.6f} {r['q']:13.6e} {r['lmin']:13.6f} "
            f"{lr:11.5f} {r['Z']/Bc:14.2e}")
        prev = r
    return loc


def main():
    beta = 2.5
    say("=" * 100)
    say("MODEL (3)  2D YM CHARACTER TOWER on one link  --  and the FLOW-CHOICE test")
    say("=" * 100)
    say(f"  beta = {beta};  S_W = beta(1-q_0);  s3 quadrature 80^3 nodes.")
    say("")
    taus = [0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0]

    say("  ---- FLOW A: the BBD/BD heat-kernel flow  (must reproduce model (1)) ----")
    mA = Tower(beta, spectrum_A(240))
    locA = report("A", mA, taus, 3.0, "lam1=3")
    say(f"  => rate(t->inf) ~ {locA[-3]:.6f} ... {locA[-1]:.6f}   (lam_1 = 3 ; KERNEL flow)")
    say(f"     sup_t C_abs = {max(r['C'] for r in curve(mA, taus, 3.0)):.6f}")
    say("")

    say("  ---- FLOW B: the 2D YM TRANSFER-MATRIX flow  (E_1 = the physical 2D YM gap) ----")
    mB = Tower(beta, spectrum_B(beta, 240))
    E1 = mB.E[1]
    locB = report("B", mB, [t * E1 for t in taus], E1, "E_1")
    say(f"  => rate(t->inf) ~ {locB[-3]:.6f} ... {locB[-1]:.6f}   vs  E_1 = {E1:.6f}  (the 2D YM gap)")
    say(f"     sup_t C_abs = {max(r['C'] for r in curve(mB, [t*E1 for t in taus], E1)):.6f}")
    say("")

    say("  ---- beta-dependence of the two rates (each fitted in units of its own candidate gap) ----")
    say(f"  {'beta':>7} {'E_1=2DYM gap':>13} {'rate_B(abs)':>12} {'rate_B/E_1':>11} "
        f"{'rate_A(abs)':>12} {'rate_A/lam_1':>13}")
    taur = np.array([1.0, 1.6, 2.5, 4.0, 6.3, 10.0])
    for b in (1.0, 2.0, 2.5, 4.0, 8.0):
        mb = Tower(b, spectrum_B(b, 200), lmax=200)
        E1b = mb.E[1]
        qB = np.array([mb.sup_x(t / E1b)["q"] for t in taur])
        ok = qB > 1e-15
        rB = -np.polyfit(taur[ok], np.log(qB[ok]), 1)[0] if ok.sum() >= 3 else np.nan
        ma = Tower(b, spectrum_A(200), lmax=200)
        qA = np.array([ma.sup_x(t / 3.0)["q"] for t in taur])
        oka = qA > 1e-15
        rA = -np.polyfit(taur[oka], np.log(qA[oka]), 1)[0] if oka.sum() >= 3 else np.nan
        say(f"  {b:7.3f} {E1b:13.6f} {rB*E1b:12.6f} {rB:11.6f} {rA*3.0:12.6f} {rA:13.6f}")
    say("")
    say("  => FLOW B's rate tracks E_1 = the 2D YM mass gap (a dynamical, beta-dependent number):")
    say("     'integrability of C requires m>0' is circular EXACTLY for that flow.")
    say("     FLOW A's rate is 3 = lam_1(SU(2)) for every beta: the one-link group gap, NOT dynamical.")
    say("     BBD sec.5 fixes the flow to be P_t = the group heat semigroup  =>  FLOW A.")
    say("")
    flush()


if __name__ == "__main__":
    main()
