import json
p=r"D:\Dev\math\research\v3\species\w7\repro\analysis.json"
a=json.load(open(p,encoding="utf-8"))
tot_hi=sum(r["cited_hi"] for r in a); tot_ok=sum(r["ok_hi"] for r in a)
print("TOTAL scripts",len(a),"hi cited",tot_hi,"matched",tot_ok)
for r in a:
    print("%-22s hi %d/%d miss %d notes=%d" % (r["script"], r["ok_hi"], r["cited_hi"], len(r["miss_hi"]), len(r["notes"])))
