"""sos3 — SOS(H-c = sum A_i^* A_i, 形式核上) 撞不撞五类 + 真障碍。
对象: Savvidy 算子 H = -D^2 + 2 s g|B| (s=±1 螺旋度), f02 §1 (1) 式。
A: 谱; B: H-c>=0 的可行性 ⟺ c<=lam_min; C: SOS⟺PSD (k=1 够) + 证书核=负模;
D: 平移 H+g|B|. 输出 w7/sos3_obstruction.log。py -3.12。
"""
import numpy as np

out = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); out.append(s)

# ---------- A. Savvidy 算子 (截断振子基 x 纵向周期 FD, dx=1) ----------
gB, N, kT = 1.0, 60, 16
a = np.diag(np.sqrt(np.arange(1, N)), 1)
Nop = a.T @ a
T1 = np.diag(-2.0 * np.ones(kT)) + np.diag(np.ones(kT - 1), 1) + np.diag(np.ones(kT - 1), -1)
T1[0, -1] = T1[-1, 0] = 1.0
L1 = -T1                                          # -d^2/dx^2, 周期
I_N, I_T = np.eye(N), np.eye(kT)
H0 = np.kron(I_T, 2 * gB * (Nop + 0.5 * I_N)) + np.kron(L1, I_N)   # -D^2
dim = N * kT
Hm = H0 - 2 * gB * np.eye(dim)                     # s=-1
ev = np.linalg.eigvalsh(Hm)
P("== A. Savvidy 谱 (gB=1, s=-1) ==")
P(f"dim={dim}  lambda_min = {ev[0]:.12f}   (解析 -g|B| = -1)   [T]+[V]")
P(f"前 6 个: {np.round(ev[:6],12)}")
P(f"0 是否在谱里: min|lambda| = {np.min(np.abs(ev)):.6f}  >0  => 0 落在谱隙内"
  f" (本格点值; 连续极限 = g|B| = 1)")
P(f"负本征值个数 = {int(np.sum(ev < -1e-12))}  (LLL 横向磁支 x k_par=0 邻域)")
P(f"s=+1 支: lambda_min = {np.linalg.eigvalsh(H0 + 2*gB*np.eye(dim))[0]:.12f}  (解析 +3g|B|)")

# ---------- B. 可行性: H - c >= 0 <==> c <= lambda_min ----------
P("\n== B. H-c >= 0 可行性扫描 ==")
P(f"{'c':>8} {'lam_min(H-c)':>16} {'>=0':>5}")
best = None
for c in [-2.0, -1.5, -1.0, -0.5, 0.0, 0.01, 0.5, 1.0]:
    lm = np.linalg.eigvalsh(Hm - c * np.eye(dim))[0]
    ok = lm >= -1e-10
    best = c if ok else best
    P(f"{c:8.2f} {lm:16.10f} {str(ok):>5}")
P(f"最大可行 c = {best:.12f} = lambda_min(H) = -g|B| < 0  =>  c>0 全部不可行 (含 c=0.01)")

# ---------- C. 稠密 A: SOS <==> PSD ----------
P("\n== C. 稠密 A: sum A_i^*A_i = H-c 的可行集 ==")
c0 = ev[0]
T = Hm - c0 * np.eye(dim)
w, V = np.linalg.eigh(T)
tol = 1e-8 * w[-1]
P(f"c = lambda_min: lam_min(T) = {w[0]:.3e}  rank_1e-8 = {int(np.sum(w > tol))}/{dim}")
Lc = np.linalg.cholesky(T + 1e-13 * np.eye(dim))
P(f"k=dim: Cholesky 残差 ||LL^T-T||/||T|| = {np.linalg.norm(Lc@Lc.T - T)/np.linalg.norm(T):.3e}")
A1 = V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.T
P(f"k=1: A=sqrt(T), ||AA^T-T||/||T|| = {np.linalg.norm(A1@A1.T - T)/np.linalg.norm(T):.3e}"
  f"  => k 不是约束, 一个稠密 A 永远够")
P(f"A psi_min: ||A psi_min|| = {np.linalg.norm(A1@V[:,0]):.2e}  => ker A = 负模方向 (证书=该本征向量)")
P(f"c = lam_min+0.5: lam_min(H-c) = {np.linalg.eigvalsh(Hm-(c0+0.5)*np.eye(dim))[0]:.10f} < 0 => 无任何因子")

# ---------- D. 平移逃逸 ----------
P("\n== D. H -> H+|gB| ==")
es = np.linalg.eigvalsh(Hm + gB * np.eye(dim))
P(f"lam_min(H+|gB|) = {es[0]:.3e} => 最佳 c'=0, 得到 H >= -|gB| (非 m>0)")
P(f"平移后谱前 3: {np.round(es[:3],6)}  (0 之上的间隔 0.1522 = 盒长 k_par 间隔; 连续极限 2g|B|)")

open("w7/sos3_obstruction.log", "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n[log] w7/sos3_obstruction.log")
