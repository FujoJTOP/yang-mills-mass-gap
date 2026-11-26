"""GP8 -- fixed figures (half-period only, S/N panel, no noise tail)."""
import numpy as np, json, os, time
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT = r'D:\Dev\math\docs\figs\gp6'; os.makedirs(OUT, exist_ok=True)
OBS=[0,3,4]; NAM=['W11','W12','W22']; NB=20
d=np.load('gh3_probe.npz'); S,Lt,N=d['S'],int(d['Lt']),int(d['n']); beta,Ls=float(d['beta']),int(d['Ls'])
A=np.ascontiguousarray(S[:,0][:,OBS,:]).astype(np.float64); del S
Abar=A.mean(axis=(0,2)); Ac=A-Abar[None,:,None]
blocks=np.array_split(np.arange(N),NB)
def build(X):
    Cs=np.empty((Lt,N,3,3))
    for t in range(Lt): Cs[t]=np.einsum('nct,ndt->ncd',X,np.roll(X,-t,axis=2),optimize=True)/Lt
    return Cs.sum(1), np.stack([Cs[:,b].sum(1) for b in blocks],0)
def gevp(Ct):
    C0=Ct[0]; w,V=np.linalg.eigh(0.5*(C0+C0.T)); k=w>1e-12*max(w.max(),1e-300)
    L=V[:,k]*np.sqrt(w[k]); P=np.linalg.solve(L,np.eye(L.shape[1]))
    return np.array([np.linalg.eigvalsh(0.5*((P@Ct[t]@P.T)+(P@Ct[t]@P.T).T))[-1] for t in range(Lt)])
tot_u,Pb_u=build(A); lt_u=gevp(tot_u); jk_u=np.array([gevp(tot_u-Pb_u[b]) for b in range(NB)])
tot_c,Pb_c=build(Ac); lt_c=gevp(tot_c); jk_c=np.array([gevp(tot_c-Pb_c[b]) for b in range(NB)])
TM=5    # t=1..4 only: S/N>=2 window (t>=5 is noise)
ts=np.arange(1,TM)
me,er,sn=[],[],[]
for t in ts:
    m=np.log(lt_c[t]/lt_c[t+1]); v=np.array([np.log(j[t]/j[t+1]) if j[t]>0 and j[t+1]>0 else np.nan for j in jk_c]); v=v[~np.isnan(v)]
    e=np.sqrt((len(v)-1)/len(v)*((v-v.mean())**2).sum()) if len(v)>2 else np.nan
    me.append(m); er.append(e); sn.append(m/e if e and e>0 else np.nan)
me,er,sn=map(np.array,(me,er,sn))
diag={NAM[i]:np.array([tot_c[t,i,i]/tot_c[0,i,i] for t in range(1,7)]) for i in range(3)}
# Fig1
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
ax[0].plot(np.arange(TM+1),lt_u[:TM+1],'o-',color='crimson',ms=5,label='un-subtracted')
ax[0].axhline(1.0,ls=':',c='gray'); ax[0].set_ylim(0.9997,1.00005)
ax[0].set_xlabel('t'); ax[0].set_ylabel(r'$\lambda_{\max}(t)$')
ax[0].set_title(r'(a) un-subtracted: $\lambda\to1$, $m_{eff}\equiv0$'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].semilogy(np.arange(TM+1),np.abs(lt_c[:TM+1]),'o-',color='navy',ms=5,label='connected')
ax[1].set_xlabel('t'); ax[1].set_ylabel(r'$|\lambda_{\max}(t)|$')
ax[1].set_title('(b) connected: genuine decay'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3,which='both')
fig.suptitle('SU(2) beta=%g, %d^3x%d, N=%d  (only t<=Lt/2 shown: C(t)=C(Lt-t))'%(beta,Ls,Lt,N),fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig1_lambda.png'),dpi=150); plt.close(fig)
# Fig2 two panels
fig,(a1,a2)=plt.subplots(2,1,figsize=(7.4,6.4),sharex=True,gridspec_kw={'height_ratios':[2,1]})
a1.errorbar(ts,me,yerr=er,fmt='o-',color='navy',capsize=3,label=r'connected GEVP $\{W_{11},W_{12},W_{22}\}$')
for nm in NAM:
    r=diag[nm]; a1.plot(np.arange(1,5),np.log(r[:4]/r[1:5]),'--',lw=1,alpha=.45,label=f'diag {nm}')
a1.axhline(1.108,ls=':',c='green',lw=1.3,label='$L_t$=12 reference (1.108)')
a1.set_ylabel(r'$m_{eff}(t)$'); a1.set_ylim(0.55,1.75); a1.legend(fontsize=7.5); a1.grid(alpha=.3)
a1.set_title(fr'connected GEVP, $\beta={beta}$, ${Ls}^3\times{Lt}$, N={N} (20-block JK)')
a2.axhspan(-1,5,color='gray',alpha=.15); a2.axhline(5,ls='--',c='crimson',lw=1,label='S/N = 5')
a2.plot(ts,sn,'s-',color='darkorange',ms=4); a2.set_yscale('log'); a2.set_ylim(1.2,150)
a2.set_xlabel('t'); a2.set_ylabel('S/N'); a2.legend(fontsize=8); a2.grid(alpha=.3,which='both')
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig2_meff.png'),dpi=150); plt.close(fig)
# Fig3
fig,ax=plt.subplots(figsize=(7.0,4.4))
for nm in NAM: ax.semilogy(np.arange(1,7),np.abs(diag[nm]),'o-',label=nm)
ax.set_xlabel('t'); ax.set_ylabel('$C(t)/C(0)$'); ax.set_title(r'connected diagonal correlators, $\beta=2.5$')
ax.legend(); ax.grid(alpha=.3,which='both'); fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig3_ct.png'),dpi=150); plt.close(fig)
# Fig4
fig,ax=plt.subplots(figsize=(7.0,3.9))
x=np.arange(3); sig=np.sqrt(np.abs([tot_c[0,i,i] for i in range(3)]))
ax.bar(x-0.2,np.abs(Abar),0.38,label=r'$|\bar O_i|$ (vacuum)',color='crimson')
ax.bar(x+0.2,sig,0.38,label=r'$\sqrt{C_{ii}(0)}$ (connected signal)',color='navy')
ax.set_yscale('log'); ax.set_xticks(x); ax.set_xticklabels(NAM)
ax.set_title('signal is $\\sim10^{-2}$ of the constant: why subtraction is mandatory')
ax.legend(); ax.grid(alpha=.3,axis='y',which='both')
for i in range(3): ax.text(x[i]-0.2,abs(Abar[i])*1.2,f'{abs(Abar[i]):.0f}',ha='center',fontsize=8)
for i in range(3): ax.text(x[i]+0.2,sig[i]*1.2,f'{sig[i]:.2f}',ha='center',fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig4_vacuum.png'),dpi=150); plt.close(fig)
json.dump(dict(meff=[[float(a),float(b)] for a,b in zip(me,er)],sn=[float(v) for v in sn],
               t=list(map(int,ts)),vac=[float(v) for v in Abar],lam_conn=[float(v) for v in lt_c],
               lam_unsub=[float(v) for v in lt_u]),open(os.path.join(OUT,'numbers.json'),'w'),indent=1)
print('figs regenerated (t<=%d)'%TM)
for t,m,e,s in zip(ts,me,er,sn): print(f' t={t}: m_eff={m:+.4f} +- {e:.4f}  S/N={s:6.2f}')


