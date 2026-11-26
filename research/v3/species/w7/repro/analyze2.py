#!/usr/bin/env python3
"""analyze2.py -- judge note-cited numbers against fresh stdout captures.

fresh output  = repro/out/<group>__<stem>.out.txt   (stdout of *my* run only)
note evidence = lines in research/v3/species/**/*.md that either name `<stem>.py`
                or carry a "[F]" fact stamp

A cited number is OK when some number in the fresh output agrees to rtol<=1e-3
(tolerances are then re-checked by hand against each note's own printed precision).
High-signal = mantissa with >=4 significant digits, or scientific notation.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SPECIES = os.path.dirname(ROOT)
OUT = os.path.join(HERE, "out")

NUM = re.compile(r"[-−+]?\d+(?:\.\d+)?(?:[eE][-+−]?\d+)?")
RTOL = 1e-3


def fval(tok):
    return float(tok.replace("−", "-"))


def nsig(tok):
    t = tok.replace("−", "-").lstrip("+-")
    if "e" in t or "E" in t:
        return 99
    mant = t
    return len([c for c in mant.split(".")[0].lstrip("0")]) + (
        len(mant.split(".")[1]) if "." in mant else 0)


def nums(text):
    out = []
    for m in NUM.finditer(text):
        try:
            out.append((m.group(0), fval(m.group(0))))
        except ValueError:
            pass
    return out


def fresh(stem, group):
    p = os.path.join(OUT, f"{group}__{stem}.out.txt")
    if os.path.exists(p):
        return open(p, encoding="utf-8", errors="replace").read()
    return ""


TOKRE = re.compile(r"[-+]?(\d+)(?:\.(\d+))?(?:[eE]([-+]?\d+))?$")


def tol_of(tok):
    """half-ulp of the precision the note itself printed"""
    m = TOKRE.match(tok.replace("−", "-"))
    if not m:
        return 5e-3
    dec = len(m.group(2)) if m.group(2) else 0
    ex = int(m.group(3)) if m.group(3) else 0
    return 0.5 * (10.0 ** -dec) * (10.0 ** ex)


def best_match(v, cands, tok):
    tol = tol_of(tok)
    b = None
    for ctok, u in cands:
        if u == v:
            return (ctok, 0.0, tol)
        d = abs(u - v)
        if b is None or d < b[1]:
            b = (ctok, d, tol)
    if b is not None:
        # agreement at the coarser of (note precision, printed output precision)
        b = (b[0], b[1], max(b[2], tol_of(b[0])))
    return b


def note_lines():
    ev = {}
    for dp, _, fns in os.walk(SPECIES):
        for f in fns:
            if not f.endswith(".md"):
                continue
            p = os.path.join(dp, f)
            txt = open(p, encoding="utf-8", errors="replace").read()
            ev[f] = txt.splitlines()
    return ev


def main():
    nl = note_lines()
    scripts = []
    for f in sorted(os.listdir(ROOT)):
        if f.endswith(".py"):
            scripts.append(("w7", f[:-3]))
    ad = os.path.join(ROOT, "audit")
    for f in sorted(os.listdir(ad)):
        if f.endswith(".py"):
            scripts.append(("audit", f[:-3]))

    report = []
    for group, stem in scripts:
        txt = fresh(stem, group)
        cands = nums(txt)
        pat = re.compile(re.escape(stem) + r"\.py")
        per_note = {}
        for nname, lines in nl.items():
            if not any(pat.search(l) for l in lines):
                continue
            hits = []
            for i, l in enumerate(lines):
                if pat.search(l) or l.strip().startswith("[F]"):
                    for tok, v in nums(l):
                        if not (-1e-6 < v < 1e-6 and "." not in tok):
                            hits.append((tok, v, i + 1, nsig(tok)))
            per_note[nname] = hits
        rows = []
        seen = set()
        for nname, hits in per_note.items():
            for tok, v, ln, sig in hits:
                k = round(v, 10)
                if k in seen:
                    continue
                seen.add(k)
                m = best_match(v, cands, tok) if cands else None
                ok = bool(m and m[1] <= m[2])
                rows.append(dict(note=nname, line=ln, tok=tok, val=v, sig=sig,
                                 ok=ok, match=(m[0] if m else None),
                                 absdiff=(m[1] if m else None),
                                 tol=(m[2] if m else None)))
        hi = [r for r in rows if r["sig"] >= 4]
        report.append(dict(script=stem, group=group, notes=sorted(per_note),
                           fresh_bytes=len(txt), fresh_nums=len(cands),
                           cited=len(rows), cited_hi=len(hi),
                           ok_hi=sum(1 for r in hi if r["ok"]),
                           miss_hi=[r for r in hi if not r["ok"]],
                           rows=rows))
    with open(os.path.join(HERE, "analysis.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=1, ensure_ascii=False)

    tot_hi = sum(r["cited_hi"] for r in report)
    tot_ok = sum(r["ok_hi"] for r in report)
    print(f"scripts={len(report)}  high-signal cited numbers={tot_hi}  matched={tot_ok}")
    for r in report:
        if not r["notes"]:
            print(f"\n-- {r['script']:22s} [{r['group']}]  no citing note; stdout {r['fresh_bytes']}B")
            continue
        flag = "" if not r["miss_hi"] else f"  <-- {len(r['miss_hi'])} unmatched"
        print(f"\n== {r['script']:22s} [{r['group']}] notes={','.join(r['notes'])}{flag}")
        print(f"   hi {r['ok_hi']}/{r['cited_hi']}  stdout={r['fresh_bytes']}B")
        for x in r["miss_hi"][:10]:
            print(f"     MISS {x['tok']:>16}  ({x['note']}:{x['line']}) nearest={x['match']}")


if __name__ == "__main__":
    main()
