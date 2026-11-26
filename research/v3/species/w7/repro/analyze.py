#!/usr/bin/env python3
"""analyze.py -- compare numbers cited in notes against fresh repro outputs.

For each w7/*.py (and w7/audit/*.py) script:
  * collect the fresh output text (stdout capture + snapshot of the .log it wrote)
  * find every note under research/v3/species that names the script
  * pull the numeric tokens from the note lines that name the script, plus lines
    carrying a "[F]" fact stamp
  * try to match each note number against the fresh text (exact string, then
    numeric |rel err| <= tol, tol taken from the number's own printed precision)

Output: repro/analysis.json + a readable dump on stdout.
"""
import json
import os
import re
import sys

W7 = os.path.dirname(os.path.abspath(__file__))          # .../w7/repro
ROOT = os.path.dirname(W7)                               # .../w7
SPECIES = os.path.dirname(ROOT)                          # .../species
OUT = os.path.join(W7, "out")
SNAP = os.path.join(W7, "snap")

NUM = re.compile(r"[-−+]?\d+(?:\.\d+)?(?:[eE][-+−]?\d+)?")


def as_float(tok):
    t = tok.replace("−", "-").replace("+", "")
    return float(t)


def sig_digits(tok):
    """digits after the decimal point in the mantissa -> implied precision"""
    t = tok.replace("−", "-").lstrip("+-")
    mant = t.split("e")[0].split("E")[0]
    return len(mant.split(".")[1]) if "." in mant else 0


def numbers_in(text):
    out = []
    for m in NUM.finditer(text):
        tok = m.group(0)
        try:
            out.append((tok, as_float(tok), m.start(), m.end()))
        except ValueError:
            pass
    return out


def fresh_text(stem, group):
    parts = []
    for name in (f"{group}__{stem}.out.txt",):
        p = os.path.join(OUT, name)
        if os.path.exists(p):
            parts.append(open(p, encoding="utf-8", errors="replace").read())
    sd = os.path.join(SNAP, f"{group}__{stem}")
    if os.path.isdir(sd):
        for f in sorted(os.listdir(sd)):
            if f.endswith((".log", ".json")):
                p = os.path.join(sd, f)
                parts.append(open(p, encoding="utf-8", errors="replace").read())
    return "\n".join(parts)


def notes():
    out = {}
    for dp, _, fns in os.walk(SPECIES):
        for f in fns:
            if f.endswith(".md"):
                p = os.path.join(dp, f)
                out[f] = open(p, encoding="utf-8", errors="replace").read()
    return out


def script_stems():
    res = []
    for f in sorted(os.listdir(ROOT)):
        if f.endswith(".py"):
            res.append(("w7", f[:-3]))
    ad = os.path.join(ROOT, "audit")
    if os.path.isdir(ad):
        for f in sorted(os.listdir(ad)):
            if f.endswith(".py"):
                res.append(("audit", f[:-3]))
    return res


def match_number(val, txt_nums, tol_default=5e-3):
    best = None
    for tok, v, _, _ in txt_nums:
        if v == val:
            return tok, 0.0
        denom = max(abs(v), abs(val), 1e-300)
        rel = abs(v - val) / denom
        if best is None or rel < best[1]:
            best = (tok, rel)
    if best and best[1] <= tol_default:
        return best
    return None


def main():
    NT = notes()
    report = []
    for group, stem in script_stems():
        fresh = fresh_text(stem, group)
        fnums = numbers_in(fresh)
        cited_notes, cited = [], []
        pat = re.escape(stem) + r"\.py"
        for nname, ntext in NT.items():
            if not re.search(pat, ntext):
                continue
            cited_notes.append(nname)
            for line in ntext.splitlines():
                if re.search(pat, line) or line.strip().startswith("[F]"):
                    cited.append((nname, line.strip()))
        # numbers on the citing lines only (keeps the set small and relevant)
        rows = []
        seen = set()
        for nname, line in cited:
            for tok, val, _, _ in numbers_in(line):
                key = (round(val, 12),)
                if key in seen:
                    continue
                seen.add(key)
                m = match_number(val, fnums)
                sd = sig_digits(tok)
                rows.append({
                    "note": nname, "note_tok": tok, "note_val": val,
                    "frac_digits": sd,
                    "match_tok": m[0] if m else None,
                    "rel": m[1] if m else None,
                    "ok": bool(m),
                })
        report.append({
            "script": stem, "group": group,
            "notes": cited_notes,
            "fresh_chars": len(fresh),
            "fresh_numbers": len(fnums),
            "cited_numbers": len(rows),
            "matched": sum(1 for r in rows if r["ok"]),
            "rows": rows,
        })
    with open(os.path.join(W7, "analysis.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=1, ensure_ascii=False)
    # readable dump
    for r in report:
        miss = [x for x in r["rows"] if not x["ok"]]
        print(f"\n### {r['script']} [{r['group']}]  notes={','.join(r['notes']) or '-'}")
        print(f"    fresh {r['fresh_chars']} chars / {r['fresh_numbers']} numbers; "
              f"cited {r['cited_numbers']}, matched {r['matched']}, unmatched {len(miss)}")
        for x in miss[:12]:
            print(f"      MISS note={x['note_tok']:>14}  (from {x['note']})")
    print("\nwrote analysis.json")


if __name__ == "__main__":
    main()
