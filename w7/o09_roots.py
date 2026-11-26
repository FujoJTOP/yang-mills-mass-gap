import numpy as np
src = open(r"D:\Dev\math\w7\o09_phicurv.py",encoding="utf-8").read()
head = src.split('say("=" * 100)\nsay("[1]')[0]
ns = {}
exec(head, ns)
body = src[src.index("def slice_a"):src.index("for th, tag, ref in")]
body = body.replace('say(','pass #(').replace('LOG.','pass # LOG.')
ns.update(RMAX=2*np.pi-1e-3, NS=2000, HH=1e-4)
exec(body, ns)
slice_a = ns["slice_a"]
def root_bs(f, lo, hi, it=60):
    flo=f(lo)
    for _ in range(it):
        mid=0.5*(lo+hi)
        if f(mid)*flo>0: lo=mid; flo=f(mid)
        else: hi=mid
    return 0.5*(lo+hi)
f = lambda a: slice_a(a, 2.5, np.pi)["Gpp"]
r1 = root_bs(f, 0.5, 0.75); r2 = root_bs(f, 2.9, 3.4)
L=[]
L.append("G''(a)=0 at beta=2.5, theta=pi:  a1 = %.5f   a2 = %.5f" % (r1, r2))
for a in (0.25,0.5,0.75,1.0,1.5,2.0):
    r = slice_a(a, 2.5, np.pi)
    L.append("a=%.2f  Phi''=%.5f  Var(force)=%.5f  E[H]=%+.5f  G''=%+.5f  Phi''-Var=%+.5f  Phi''/Var=%.2f"
             % (a, r["Phi2"], r["varF"], r["H"], r["Gpp"], r["Phi2"]-r["varF"], r["Phi2"]/r["varF"]))
dev = []
for a in (0.25,0.5,1.0,1.5,2.0):
    r = slice_a(a, 2.5, np.pi)
    dev.append(abs(r["Gpp"]-(r["Phi2"]+r["H"]-r["varF"])))
L.append("identity G''=Phi''+E[H_vv]-Var(force): max |dev| = %.2e" % max(dev))
print("\n".join(L))
open(r"D:\Dev\math\w7\o09_roots.log","w",encoding="utf-8").write("\n".join(L)+"\n")
