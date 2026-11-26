# rat3_a7.py -- exact count a_7 of rooted connected plaquette sets (4D Wilson plaquette graph, Delta=20)
#
# Object (pinned by cor:
#   V = plaquettes of Z^4 (torus T^4_L used as a finite stand-in, L = nx),
#   adjacency: two plaquettes adjacent iff they share a link  =>  deg = 4*(2(d-1)-1) = 20 (d=4)
#   a_n = #{ Gamma ⊂ V : p0 ∈ Gamma, Gamma connected (plaquette adjacency), |Gamma| = n }
#   p0 = plaquette at origin with orientation (0,1)  -> id 0
#
# Method (exact, no sampling, no hashing of the answer):
#   level-by-level incremental growth + global sort/dedup on exact 126-bit keys
#   (sorted list of 7 vertex ids, 18 bits each => injective encoding).
#   Level 7 is counted in `shards` sequential passes (memory-safe); each pass is
#   parallel over the level-6 list, then sorted and deduped.
#
# Why nx=14 gives the true Z^4 number for n<=7:
#   any connected n<=7 set is a walk of <=6 steps, so its integer coordinate
#   spread is <=6; torus wrap-adjacency needs a spread of nx-1 >= 13, which is
#   unreachable => the finite torus adds no edges among such sets, and the
#   rooted count is translation-free.  (nx=6, the project's usual harness, is an
#   UPPER bound for the same reason: it *can* wrap; both are reported.)
import os, subprocess, sys, time

RUST = r'''
use std::env;
use std::time::Instant;
use std::thread;

#[inline(always)]
fn wr(v: i64, nx: i64) -> usize { (((v % nx) + nx) % nx) as usize }

#[inline(always)]
fn ocode(a: usize, b: usize) -> usize {
    match (a, b) { (0,1)=>0, (0,2)=>1, (0,3)=>2, (1,2)=>3, (1,3)=>4, (2,3)=>5, _=>99 }
}

#[inline(always)]
fn pid(x: [i64; 4], a: usize, b: usize, nx: i64) -> u32 {
    let (a, b) = if a < b { (a, b) } else { (b, a) };
    let n = nx as usize;
    let pos = ((wr(x[3], nx) * n + wr(x[2], nx)) * n + wr(x[1], nx)) * n + wr(x[0], nx);
    (pos * 6 + ocode(a, b)) as u32
}

fn build_nb(nx: i64) -> (Vec<[u32; 20]>, usize) {
    let n = nx as usize;
    let nsites = n * n * n * n;
    let v = nsites * 6;
    let mut nb = vec![[u32::MAX; 20]; v];
    for id in 0..v {
        let o = id % 6;
        let pos = id / 6;
        let x0 = (pos % n) as i64; let p1 = pos / n;
        let x1 = (p1 % n) as i64; let p2 = p1 / n;
        let x2 = (p2 % n) as i64; let x3 = (p2 / n) as i64;
        let x = [x0, x1, x2, x3];
        let (mu, nu) = match o { 0=>(0usize,1usize), 1=>(0,2), 2=>(0,3), 3=>(1,2), 4=>(1,3), _=>(2,3) };
        let mut links = [([0i64; 4], 0usize); 4];
        links[0] = (x, mu);
        links[1] = (x, nu);
        let mut y = x; y[mu] += 1; links[2] = (y, nu);
        let mut y2 = x; y2[nu] += 1; links[3] = (y2, mu);
        let mut out = [0u32; 20];
        let mut k = 0usize;
        for &(yy, d) in links.iter() {
            for dd in 0..4usize {
                if dd == d { continue; }
                let c1 = pid(yy, d, dd, nx);
                let mut y3 = yy; y3[dd] -= 1;
                let c2 = pid(y3, d, dd, nx);
                for c in [c1, c2] {
                    if c as usize != id { out[k] = c; k += 1; }
                }
            }
        }
        if k != 20 { panic!("degree {}", k); }
        for i in 0..20 { for j in (i + 1)..20 { if out[i] == out[j] { panic!("dup nbr id {}", id); } } }
        nb[id] = out;
    }
    (nb, v)
}

#[inline(always)]
fn shard_of(k: u128, m: u64) -> u64 {
    let mut z = (k as u64) ^ ((k >> 64) as u64).wrapping_mul(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^= z >> 31;
    z % m
}

fn cand_keys(key: u128, nset: usize, nb: &[[u32; 20]], filt: Option<(u64, u64)>, out: &mut Vec<u128>) -> usize {
    let mut b = [0u32; 7];
    for i in 0..nset { b[i] = ((key >> (18 * i)) & 0x3FFFF) as u32; }
    let mut cd = [0u32; 140];
    let mut m = 0usize;
    for i in 0..nset {
        let vv = b[i] as usize;
        let row = &nb[vv];
        for j in 0..20 {
            let w = row[j];
            let mut inside = false;
            for t in 0..nset { if b[t] == w { inside = true; break; } }
            if inside { continue; }
            let mut seen = false;
            for t in 0..m { if cd[t] == w { seen = true; break; } }
            if !seen { cd[m] = w; m += 1; }
        }
    }
    for idx in 0..m {
        let w = cd[idx];
        let mut t = nset;
        while t > 0 && b[t - 1] > w { b[t] = b[t - 1]; t -= 1; }
        b[t] = w;
        let mut k2 = 0u128;
        for i in 0..(nset + 1) { k2 |= (b[i] as u128) << (18 * i); }
        let keep = match filt { None => true, Some((s, sh)) => shard_of(k2, sh) == s };
        if keep { out.push(k2); }
        let mut i = t;
        while i < nset { b[i] = b[i + 1]; i += 1; }
    }
    m
}

fn gen_level(prev: &Vec<u128>, nset: usize, nb: &Vec<[u32; 20]>, threads: usize,
             filt: Option<(u64, u64)>) -> (Vec<u128>, usize) {
    let plen = prev.len();
    let mut parts: Vec<Vec<u128>> = Vec::new();
    let mut mtot = 0usize;
    thread::scope(|sc| {
        let mut hs = Vec::new();
        for c in 0..threads {
            let lo = plen * c / threads;
            let hi = plen * (c + 1) / threads;
            let pv: &Vec<u128> = prev;
            let nbr: &Vec<[u32; 20]> = nb;
            hs.push(sc.spawn(move || {
                let mut out: Vec<u128> = Vec::with_capacity(1 << 16);
                let mut mt = 0usize;
                for i in lo..hi { mt += cand_keys(pv[i], nset, nbr, filt, &mut out); }
                (out, mt)
            }));
        }
        for h in hs {
            let (v, mt) = h.join().unwrap();
            mtot += mt;
            parts.push(v);
        }
    });
    let tot: usize = parts.iter().map(|p| p.len()).sum();
    let mut all: Vec<u128> = Vec::with_capacity(tot);
    for p in parts { all.extend_from_slice(&p); }
    (all, mtot)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let nx: i64 = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(14);
    let nmax: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(7);
    let shards: u64 = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(4);
    let threads: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(16);
    let t_all = Instant::now();
    let (nb, v) = build_nb(nx);
    println!("[graph] torus T^4_{} : V = {} plaquettes, deg = 20 (share-a-link), p0 = id 0", nx, v);
    let mut cur: Vec<u128> = vec![0u128];
    println!("[lv 1] a_1 = 1   (cum {:.2?})", t_all.elapsed());
    for lvl in 2..=nmax {
        let t0 = Instant::now();
        let plen = cur.len();
        if lvl < nmax {
            let (mut cand, mtot) = gen_level(&cur, lvl - 1, &nb, threads, None);
            let ncand = cand.len();
            cand.sort_unstable();
            cand.dedup();
            println!("[lv {}] a_{} = {}   cand = {} (avg bdry {:.3})  dedup {:.2}x   {:.2?}",
                     lvl, lvl, cand.len(), ncand, mtot as f64 / plen as f64,
                     ncand as f64 / cand.len() as f64, t0.elapsed());
            cur = cand;
        } else {
            // memory-safe sharded counting
            let mut total = 0usize;
            let mut mtot_all = 0usize;
            for s in 0..shards {
                let ts = Instant::now();
                let (mut part, mtot) = gen_level(&cur, lvl - 1, &nb, threads, Some((s, shards)));
                mtot_all += mtot;
                let nlen = part.len();
                part.sort_unstable();
                part.dedup();
                total += part.len();
                println!("   [shard {}/{}] uniq = {}  (raw {})  {:.2?}", s, shards, part.len(), nlen, ts.elapsed());
            }
            println!("[lv {}] a_{} = {}   cand = {} (avg bdry {:.3})  dedup {:.2}x   {:.2?}",
                     lvl, lvl, total, mtot_all, mtot_all as f64 / plen as f64,
                     mtot_all as f64 / total as f64, t0.elapsed());
        }
    }
    println!("[done] total {:.2?}", t_all.elapsed());
}
'''

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(os.environ.get("TEMP", "."), "rat3_build")
os.makedirs(BUILD, exist_ok=True)
src = os.path.join(BUILD, "a7.rs")
exe = os.path.join(BUILD, "a7.exe")
with open(src, "w", encoding="utf-8") as f:
    f.write(RUST)

t0 = time.time()
subprocess.run(["rustc", "-O", "-C", "target-cpu=native", "-o", exe, src], check=True)
print("[compile] %.1f s -> %s" % (time.time() - t0, exe), flush=True)

for args in ([("6", "6", "1", "16")] + [("14", "7", "4", "16")]):
    print("\n=== run nx=%s nmax=%s shards=%s threads=%s ===" % args, flush=True)
    t0 = time.time()
    subprocess.run([exe] + list(args), check=True)
    print("=== wall %.1f s ===" % (time.time() - t0), flush=True)
