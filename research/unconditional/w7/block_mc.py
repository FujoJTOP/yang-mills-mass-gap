"""
块内 q 场的 SU(2) Wilson MC — 输出整块 q 场（供 F2 / RTC 检验）
d=4, 周期边界, Metropolis + (μ, 宇称) 棋盘并行
SU(2) 四元数 (a,b) ∈ C²: U = [[a, b], [-conj(b), conj(a)]], tr U = 2 Re a
U 表示为 tuple of 4 pairs: U[mu] = (a,b), 每个 shape (L,L,L,L)
"""
import argparse, time, sys
import numpy as np

def qmul(x, y):
    a1, b1 = x; a2, b2 = y
    return (a1*a2 - b1*np.conj(b2), a1*b2 + b1*np.conj(a2))

def qdag(x):
    a, b = x
    return (np.conj(a), -b)

def qtr2(x):                       # tr U = 2 Re a
    return 2.0*x[0].real

def sh(L_):
    return L_

def shift(Umu, mu, k):             # 沿 mu 平移 k 步（周期）：out[x] = in[x+k]
    return (np.roll(Umu[0], k, axis=mu), np.roll(Umu[1], k, axis=mu))

def init_u(L, rng):
    out = []
    for _ in range(4):
        a = rng.normal(size=(L,L,L,L)) + 1j*rng.normal(size=(L,L,L,L))
        b = rng.normal(size=(L,L,L,L)) + 1j*rng.normal(size=(L,L,L,L))
        n = np.sqrt(np.abs(a)**2 + np.abs(b)**2)
        out.append((a/n, b/n))
    return tuple(out)

def staple(U, mu):
    """Σ(x,μ) = Σ_{ν≠μ}[ U_ν(x+μ) U_μ(x+ν)† U_ν(x)† + U_ν(x+μ−ν)† U_μ(x−ν)† U_ν(x−ν) ]"""
    s = None
    for nu in range(4):
        if nu == mu: continue
        unu  = U[nu]; umu = U[mu]
        t1 = qmul(qmul(shift(unu, mu, -1), qdag(shift(umu, nu, -1))), qdag(unu))
        u_nu_xmu_mnu = shift(shift(unu, mu, -1), nu, +1)
        u_mu_xmnu    = shift(umu, nu, +1)
        u_nu_xmnu    = shift(unu, nu, +1)
        t2 = qmul(qmul(qdag(u_nu_xmu_mnu), qdag(u_mu_xmnu)), u_nu_xmnu)
        s = (t1[0]+t2[0], t1[1]+t2[1]) if s is None else (s[0]+t1[0]+t2[0], s[1]+t1[1]+t2[1])
    return s

def sum6(U, mu):
    """Re tr(U_μ Σ) = Σ_{6 plaquettes 含 (x,μ)} tr(U_□)"""
    return qtr2(qmul(U[mu], staple(U, mu)))

def propo(eps, rng, shape):
    n = rng.normal(size=shape+(3,))
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    c, s = np.cos(eps), np.sin(eps)
    return (c + 1j*s*n[...,0], s*(n[...,1] + 1j*n[...,2]))

def sweep(U, beta, L, eps, rng, masks):
    """一次 sweep。优化依据：宇称棋盘 pass 只改 U_mu(x) (x 在 mask 上)，
       而 staple(x,mu) 只依赖 U_mu(x±nu) 与 U_nu(·)——x±nu 宇称相反、U_nu 不同方向
       => 本 pass 内 staple 前后【完全相同】。故只算一次。"""
    at = ac = 0
    hb = beta/2.0
    for mu in range(4):
        A = U[mu][0]; B = U[mu][1]
        for par in (0,1):
            mask = masks[par]
            S = staple(U, mu)                      # 前后相同，只算一次
            sa, sb = S
            dU = propo(eps, rng, (L,L,L,L))
            # Un = dU * U
            a0, b0 = A, B
            ua = dU[0]*a0 - dU[1]*np.conj(b0)
            ub = dU[0]*b0 + dU[1]*np.conj(a0)
            # Re tr(U S) = 2 Re(a*s - b*conj(sb))，逐点
            def f2(a,b):
                return 2.0*np.real(a*sa - b*np.conj(sb))
            dS = -hb*(f2(ua,ub) - f2(a0,b0))
            acc = (np.log(rng.random(size=(L,L,L,L))) < -dS) & mask
            na = int(acc.sum())
            if na:
                A[acc] = ua[acc]; B[acc] = ub[acc]   # 原地写，不做整块 where
            at += na; ac += int(mask.sum())
    return U, at, ac

def qfield(U):
    L = U[0][0].shape[0]
    tot = np.zeros((L,L,L,L)); cnt = 0
    for mu in range(4):
        for nu in range(mu+1,4):
            P = qmul(qmul(qmul(U[mu], shift(U[nu], mu, -1)), qdag(shift(U[mu], nu, -1))), qdag(U[nu]))
            tot += 1.0 - 0.5*qtr2(P); cnt += 1
    return tot/cnt

def total_plaq_re(U):
    s = 0.0
    for mu in range(4):
        for nu in range(mu+1,4):
            P = qmul(qmul(qmul(U[mu], shift(U[nu], mu, -1)), qdag(shift(U[mu], nu, -1))), qdag(U[nu]))
            s += qtr2(P).sum()
    return s

def selftest(L=4, beta=2.5, seed=7, ntests=6):
    """逐 link 核 ΔS：staple 形式 vs 全格重算"""
    rng = np.random.default_rng(seed)
    U = init_u(L, rng)
    worst = 0.0
    for _ in range(ntests):
        mu = int(rng.integers(0,4)); x = tuple(int(rng.integers(0,L)) for _ in range(4))
        eps = float(rng.uniform(0.1,0.6))
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        c, s = np.cos(eps), np.sin(eps)
        d0 = c + 1j*s*n[0]; d1 = s*(n[1] + 1j*n[2])          # SU(2) 元素，(a,b)
        a0 = U[mu][0][x]; b0 = U[mu][1][x]
        ua = d0*a0 - d1*np.conj(b0); ub = d0*b0 + d1*np.conj(a0)   # U' = d*U
        Ub = list(U); A = U[mu][0].copy(); B = U[mu][1].copy()
        A[x] = ua; B[x] = ub
        Ub[mu] = (A,B); Ub = tuple(Ub)
        d1_ = -(beta/2.0)*(total_plaq_re(Ub) - total_plaq_re(U))
        S_ = staple(U, mu); sa = S_[0][x]; sb = S_[1][x]
        f = lambda aa, bb: 2.0*np.real(aa*sa - bb*np.conj(sb))
        d2_ = -(beta/2.0)*(f(ua,ub) - f(a0,b0))
        worst = max(worst, abs(d1_-d2_))
        print('  mu=%d x=%s eps=%.2f  dS_brute=%+.10f  dS_staple=%+.10f  diff=%.2e'
              % (mu,x,eps,d1_,d2_,abs(d1_-d2_)))
    print('SELFTEST worst = %.3e  => %s' % (worst, 'PASS' if worst<1e-9 else 'FAIL'))
    return worst < 1e-9

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--sweeps', type=int, default=20000)
    ap.add_argument('--therm', type=int, default=3000)
    ap.add_argument('--every', type=int, default=1)
    ap.add_argument('--eps', type=float, default=0.35)
    ap.add_argument('--seed', type=int, default=2026)
    ap.add_argument('--out', type=str, default='block_mc_out.npz')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--store-fields', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    rng = np.random.default_rng(a.seed)
    L, beta = a.L, a.beta
    idxs = tuple(np.meshgrid(*[np.arange(L)]*4, indexing='ij'))
    masks = [((sum(idxs) % 2) == p) for p in (0,1)]
    U = init_u(L, rng)
    t0 = time.time(); at = ac = 0
    N = L**4
    # 累加器（不预先减均值：k=0 的抵消必须在【平均之后】做，否则 χ_q 被抹掉）
    acc_qhat = np.zeros((L,L,L,L))                 # Σ |FFT(q)|²/N
    acc_q = 0.0; acc_q2 = 0.0; acc_qb = 0.0; acc_qb2 = 0.0
    qs, sws, pls = [], [], []
    M = 0
    for sw in range(a.therm + a.sweeps):
        U, x_, y_ = sweep(U, beta, L, a.eps, rng, masks); at += x_; ac += y_
        if sw >= a.therm and (sw - a.therm) % a.every == 0:
            qf = qfield(U)
            if a.store_fields: qs.append(qf.astype(np.float32))
            sws.append(sw); pls.append(float(qf.mean()))
            qh = np.fft.fftn(qf)
            acc_qhat += (np.abs(qh)**2)/N
            qb = qf.mean()
            acc_q += qf.sum(); acc_q2 += (qf**2).sum()
            acc_qb += qb; acc_qb2 += qb*qb
            M += 1
        if sw % max(1,(a.therm+a.sweeps)//20) == 0:
            print('[%s] sw=%d/%d plaq=%.6f acc=%.3f %.1fs'
                  % (time.strftime('%H:%M:%S'), sw, a.therm+a.sweeps,
                     float(qfield(U).mean()), at/max(ac,1), time.time()-t0), flush=True)
    qs = np.array(qs) if qs else np.zeros((0,L,L,L,L),np.float32); M = max(M,1)
    qbar = acc_qb/M
    Ck = acc_qhat/M
    Ck[0,0,0,0] -= N*qbar*qbar                      # 连通化：只在 k=0 扣 N<q>²
    Cz = np.real(np.fft.ifftn(Ck))
    varq   = acc_q2/(M*N) - (acc_q/(M*N))**2        # 单点 Var(q)
    varbar = acc_qb2/M - qbar**2                    # Var(q̄)
    chi_q  = N*varbar                                # = Ĉ(0) = Σ_z C_q(z)
    np.savez_compressed(a.out, L=L, beta=beta, q_fields=qs, sweeps=np.array(sws),
                        plaq=np.array(pls), Ck=Ck, Cz=Cz,
                        varq=varq, varbar=varbar, chi_q=chi_q, M=M, N=N,
                        acc=at/max(ac,1), eps=a.eps, seed=a.seed)
    print('DONE L=%d M=%d cfgs  q_fields=%s' % (L, M, qs.shape))
    print('  plaq mean   = %.8f' % float(np.mean(pls)))
    print('  Var(q)      = %.8e' % varq)
    print('  Var(q_bar)  = %.8e' % varbar)
    print('  chi_q = N*Var(q_bar) = %.8e   (should equal Ck[0])' % chi_q)
    print('  Ck[0,0,0,0] = %.8e   max_k Ck = %.8e   Cz(0) = %.8e'
          % (Ck[0,0,0,0], Ck.max(), Cz[0,0,0,0]))
    print('  check: sum_z Cz / Ck[0] = %.9f' % (Cz.sum()/Ck[0,0,0,0] if abs(Ck[0,0,0,0])>0 else float('nan')))
    print('  time = %.1fs -> %s' % (time.time()-t0, a.out))

