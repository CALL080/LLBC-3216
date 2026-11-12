"""Minimum weights for consecutive round counts, starting each search at the previous minimum
(B_r >= B_{r-1}: the first r-1 rounds of an r-round characteristic form an (r-1)-round one).
For diff, --expect r:w,... runs W = w-1 (must be unsat) and W = w (must be sat) to confirm known values.

Run: python chain.py lin 16 8 4 --from 2 --to 16
     python chain.py diff 16 8 4 --from 2 --to 10 --expect 2:4,3:8,...
"""
import argparse, subprocess, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("n"); ap.add_argument("b1"); ap.add_argument("b2")
    ap.add_argument("--from", dest="lo", type=int, default=1); ap.add_argument("--to", dest="hi", type=int)
    ap.add_argument("--expect", default="")
    a = ap.parse_args()
    expect = dict(tuple(map(int, kv.split(":"))) for kv in a.expect.split(",") if kv)
    prev = 0
    for r in range(a.lo, a.hi + 1):
        start = expect[r] - 1 if r in expect else prev
        out = subprocess.run([sys.executable, os.path.join(HERE, "run_min.py"), a.mode, a.n, a.b1, a.b2, str(r),
                              "--start", str(start)], capture_output=True, text=True)
        print(out.stdout, end="", flush=True)
        if out.returncode != 0:
            print(out.stderr, flush=True)
            return
        last = out.stdout.strip().splitlines()[-1]
        prev = int(last.split("W<=")[1].split(":")[0])
        if r in expect and prev != expect[r]:
            print(f"MISMATCH rounds={r}: found {prev}, expected {expect[r]}", flush=True)


if __name__ == "__main__":
    main()
