"""GP7 -- figures + numbers for the standalone note on vacuum subtraction in the 0++ GEVP."""
import numpy as np, json, os, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT = r'D:\Dev\math\docs\figs\gp6'
os.makedirs(OUT, exist_ok=True)
OBS = [0,3,4]; NAM = ['W11','W12','W22']; NB = 20
d = np.load('gh3_probe.npz'); S, Lt, N = d['S'], int(d['Lt']), int(d['n'])
beta, Ls = float(d['beta']), int(d['Ls'])
A = np.ascontiguousarray(S[:,0][:,OBS,:]).astype(np.float64); del S
Abar = A.mean(axis=(0,2)); Ac = A - Abar[None,:,None]
blocks = np.array_split(np.arange(N), NB)
def build(X):
    Cs = np.empty((Lt, N, 3, 3))
    for t in range(Lt):
        Cs[t] = np.einsum('nct,ndt->ncd', X, np.roll(X,-t,axis=2), optimize=True)/Lt
    Pb = np.stack([Cs[:,b].sum(1) for b in blocks],0); tot = Cs.sum(1)
    return tot, Pb
def gevp(Ct):
    C0 = Ct[0]; w,V = np.linalg.eigh(0.5*(C0+C0.T))
    k = w > 1e-12*max(w.max(),1e-300); L = V[:,k]*np.sqrt(w[k])
    P = np.linalg.solve(L, np.eye(L.shape[1]))
    return np.array([np.linalg.eigvalsh(0.5*((P@Ct[t]@P.T)+(P@Ct[t]@P.T).T))[-1] for t in range(Lt)])
def meff_err(lt, jk):
    out=[]
    for t in range(1, Lt-1):
        if lt[t]<=0 or lt[t+1]<=0: out.append((np.nan,np.nan)); continue
        me = np.log(lt[t]/lt[t+1])
        v = np.array([np.log(j[t]/j[t+1]) if j[t]>0 and j[t+1]>0 else np.nan for j in jk])
        v = v[~np.isnan(v)]
        e = np.sqrt((len(v)-1)/len(v)*((v-v.mean())**2).sum()) if len(v)>2 else np.nan
        out.append((me,e))
    return out
t1=time.time()
tot_u, Pb_u = build(A);  lt_u = gevp(tot_u); jk_u = np.array([gevp(tot_u-Pb_u[b]) for b in range(NB)])
tot_c, Pb_c = build(Ac); lt_c = gevp(tot_c); jk_c = np.array([gevp(tot_c-Pb_c[b]) for b in range(NB)])
eu = meff_err(lt_u, jk_u); ec = meff_err(lt_c, jk_c)
diag = {NAM[i]: [float(tot_c[t,i,i]/tot_c[0,i,i]) for t in range(1,7)] for i in range(3)}
ts = np.arange(1, Lt-1)
# ---- Fig 1 ----
fig, ax = plt.subplots(1,2, figsize=(11,4.2))
ax[0].plot(np.arange(Lt), lt_u, 'o-', color='crimson', label='un-subtracted (vacuum kept)')
ax[0].axhline(1.0, ls=':', c='gray'); ax[0].set_ylim(0.999, 1.0005)
ax[0].set_xlabel('t'); ax[0].set_ylabel(r'$\lambda_{\max}(t)$'); ax[0].set_title('(a) un-subtracted: $\lambda\\to 1$, $m_{eff}\\equiv 0$')
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].semilogy(np.arange(Lt), np.abs(lt_c), 'o-', color='navy', label='connected (vacuum subtracted)')
ax[1].set_xlabel('t'); ax[1].set_ylabel(r'$|\lambda_{\max}(t)|$'); ax[1].set_title('(b) connected: genuine decay')
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which='both')
fig.suptitle(fr'SU(2) $\beta={beta}$, ${Ls}^3\times{Lt}$, N={N}: the vacuum constant dominates the GEVP', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig1_lambda.png'), dpi=150); plt.close(fig)
# ---- Fig 2 ----
fig, ax = plt.subplots(figsize=(7.2,4.6))
m = np.array([x[0] for x in ec]); e = np.array([x[1] for x in ec])
ax.errorbar(ts, m, yerr=e, fmt='o-', color='navy', capsize=3, label='connected GEVP $\\{W_{11},W_{12},W_{22}\\}$')
for i,nm in enumerate(NAM):
    r = np.array(diag[nm]); mm = np.log(r[:-1]/r[1:])
    ax.plot(np.arange(1,6), mm, 's--', ms=4, alpha=.65, label=f'naive diag {nm}')
ax.axhline(1.108, ls=':', c='green', lw=1.4, label='GH3 $L_t$=12 value (1.108)')
ax.axvspan(4, Lt-2, color='gray', alpha=.12)
ax.text(4.4, 1.95, 'S/N < 2', fontsize=9, color='gray')
ax.set_xlabel('t'); ax.set_ylabel(r'$m_{eff}(t)=\ln[\lambda(t)/\lambda(t+1)]$')
ax.set_title(fr'connected GEVP at $\beta={beta}$, ${Ls}^3\times{Lt}$, N={N} (20-block JK)')
ax.legend(fontsize=7.5); ax.grid(alpha=.3); ax.set_ylim(-0.9, 2.6)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig2_meff.png'), dpi=150); plt.close(fig)
# ---- Fig 3 ----
fig, ax = plt.subplots(figsize=(7.0,4.5))
for i,nm in enumerate(NAM):
    r = np.array(diag[nm]); ax.semilogy(np.arange(1,7), np.abs(r), 'o-', label=nm)
ax.set_xlabel('t'); ax.set_ylabel('$C(t)/C(0)$  (connected, diagonal)')
ax.set_title('connected diagonal correlators $\\beta=2.5$')
ax.legend(); ax.grid(alpha=.3, which='both')
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig3_ct.png'), dpi=150); plt.close(fig)
# ---- Fig 4 ----
fig, ax = plt.subplots(figsize=(7.0,4.0))
x = np.arange(3); sig = np.sqrt(np.abs([tot_c[0,i,i] for i in range(3)]))
ax.bar(x-0.2, np.abs(Abar), 0.38, label=r'$|\bar O_i|$ (vacuum)', color='crimson')
ax.bar(x+0.2, sig, 0.38, label=r'$\sqrt{C_{ii}(0)}$ (connected signal)', color='navy')
ax.set_yscale('log'); ax.set_xticks(x); ax.set_xticklabels(NAM)
ax.set_title('why the subtraction matters: signal is $10^{-2}$ of the constant')
ax.legend(); ax.grid(alpha=.3, axis='y', which='both')
for i in range(3): ax.text(x[i]-0.2, abs(Abar[i])*1.15, f'{abs(Abar[i]):.0f}', ha='center', fontsize=8)
for i in range(3): ax.text(x[i]+0.2, sig[i]*1.15, f'{sig[i]:.2f}', ha='center', fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig4_vacuum.png'), dpi=150); plt.close(fig)
json.dump(dict(beta=beta, Ls=Ls, Lt=Lt, N=N, vac=[float(v) for v in Abar],
               lam_unsub=[float(v) for v in lt_u], lam_conn=[float(v) for v in lt_c],
               meff_conn=[[float(a), float(b)] for a,b in ec],
               meff_unsub=[[float(a), float(b)] for a,b in eu],
               diag=diag, wall=time.time()-t1), open(os.path.join(OUT,'numbers.json'),'w'), indent=1)
print('figs ->', OUT); print('wall %.1fs'%(time.time()-t1))
for t,(a,b) in enumerate(ec[:6],1): print(f' t={t}: m_eff={a:.4f} +- {b:.4f}')
