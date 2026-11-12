"""Consecutive round counts with par_min.py, each starting at the previous minimum
(and at least the additive bound max_a B_a + B_{r-a} from the minima found so far).

Run: python seq_par.py jdiff 16 8 4 8 10 --jobs 6 --known 2:3,3:6,4:9,5:12,6:21,7:24
"""
import argparse, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("n"); ap.add_argument("b1"); ap.add_argument("b2")
    ap.add_argument("lo", type=int); ap.add_argument("hi", type=int)
    ap.add_argument("--jobs", default="6"); ap.add_argument("--known", default="")
    a = ap.parse_args()
    B = {1: 0}
    B.update(dict(tuple(map(int, kv.split(":"))) for kv in a.known.split(",") if kv))
    for r in range(a.lo, a.hi + 1):
        start = max([B.get(r - 1, 0)] + [B[x] + B[r - x] for x in range(1, r) if x in B and r - x in B])
        print(f"rounds={r}: start at {start}", flush=True)
        p = subprocess.run([sys.executable, os.path.join(HERE, "par_min.py"), a.mode, a.n, a.b1, a.b2, str(r),
                            "--start", str(start), "--jobs", a.jobs], capture_output=True, text=True)
        print(p.stdout, end="", flush=True)
        m = re.search(r"MIN .*: (\d+)", p.stdout)
        if not m:
            print("no minimum:", p.stderr[-500:], flush=True)
            return
        B[r] = int(m.group(1))


if __name__ == "__main__":
    main()
