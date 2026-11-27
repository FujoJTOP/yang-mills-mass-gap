"""HR4 inventory -- headers only (mmap), zero MC, zero heavy load.

Question: which npz under research/ carry a *usable* C(t) sequence (t = time separation
on a T^3 x R / T^4 cylinder), and for which (Ls, Lt, N, seed)?

A file qualifies if it has either
  (a) 'S' with ndim>=2 whose LAST axis is the time axis (per-slice observable),
  (b) a 1-D 'C'/'corr'/'Ct' of length >= 4 (pre-computed correlator).
Nothing is loaded beyond the header; the arrays are never read.

Run: py -3.12 hr4_inventory.py
"""
import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))  # D:\Dev\math

OUT = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


def scalar(d, *names):
    for n in names:
        if n in d.files:
            a = d[n]
            if a.shape == ():
                return a.item()
    return None


def main():
    rows = []
    nfile = 0
    for dirpath, _, files in os.walk(os.path.join(ROOT, 'research')):
        for fn in sorted(files):
            if not fn.endswith('.npz'):
                continue
            fp = os.path.join(dirpath, fn)
            nfile += 1
            try:
                d = np.load(fp, mmap_mode='r')
            except Exception as e:                                # noqa: BLE001
                rows.append(dict(file=os.path.relpath(fp, ROOT), err=repr(e)))
                continue
            with d:
                keys = list(d.files)
                kshape = {k: (str(d[k].dtype), tuple(d[k].shape)) for k in keys}
                has_t = None
                if 'S' in keys and d['S'].ndim >= 2:
                    has_t = 'S'
                else:
                    for k in ('C', 'Ct', 'corr', 'Cfull', 'C_t'):
                        if k in keys and d[k].ndim == 1 and d[k].shape[0] >= 4:
                            has_t = k
                            break
                if has_t is None:
                    rows.append(dict(file=os.path.relpath(fp, ROOT), tseq=None,
                                     keys=kshape))
                    continue
                Ls = scalar(d, 'Ls', 'L')
                Lt = scalar(d, 'Lt')
                N = scalar(d, 'N', 'n')
                seed = scalar(d, 'seed', 'seed0')
                beta = scalar(d, 'beta')
                if Lt is None:
                    Lt = int(kshape[has_t][1][-1]) if has_t == 'S' else int(kshape[has_t][1][0])
                if N is None and has_t == 'S':
                    N = int(kshape[has_t][1][0])
                rows.append(dict(file=os.path.relpath(fp, ROOT), tseq=has_t, Ls=Ls, Lt=Lt,
                                 N=N, seed=seed, beta=beta, keys=kshape))

    witht = [r for r in rows if r.get('tseq')]
    nott = [r for r in rows if not r.get('tseq') and 'err' not in r]
    errs = [r for r in rows if 'err' in r]

    P("=" * 108)
    P(f"[A] npz files scanned: {nfile}   with t-sequence: {len(witht)}   without: {len(nott)}"
      f"   header errors: {len(errs)}")
    P("=" * 108)
    P(f"{'file':<52} {'key':>6} {'Ls':>4} {'Lt':>4} {'N':>8} {'seed':>10} {'beta':>6}"
      f" {'t (raw)':>7} {'t (p=0)':>7}")
    P("-" * 108)
    for r in sorted(witht, key=lambda x: (str(x['Ls']), str(x['Lt']), x['file'])):
        Lsr, Ltr = r['Ls'], r['Lt']
        raw = f"{Ltr:>3}/{Ltr:>3}" if Ltr else '?'
        p0 = f"{Ltr:>3}/{Ltr:>3}" if Ltr else '?'
        P(f"{r['file']:<52} {r['tseq']:>6} {str(Lsr):>4} {str(Ltr):>4} "
          f"{str(r['N']):>8} {str(r['seed']):>10} {str(r['beta']):>6} {raw:>7} {p0:>7}")

    P("")
    P("=" * 108)
    P("[B] npz WITHOUT a t-sequence (not usable for the two-point test)")
    P("=" * 108)
    for r in sorted(nott, key=lambda x: x['file']):
        P(f"   {r['file']:<62} keys={list(r['keys'])}")
    for r in errs:
        P(f"   HEADER-FAIL {r['file']}  {r['err']}")

    json.dump([{k: (v if not isinstance(v, dict) else {kk: list(vv) if isinstance(vv, tuple) else vv
                                                        for kk, vv in v.items()})
                for k, v in r.items()} for r in rows],
              open(os.path.join(HERE, 'hr4_inventory.json'), 'w'), indent=1, default=str)
    open(os.path.join(HERE, 'hr4_inventory.log'), 'w').write("\n".join(OUT) + "\n")
    P("\n[ok] wrote hr4_inventory.log / hr4_inventory.json")


if __name__ == '__main__':
    main()
