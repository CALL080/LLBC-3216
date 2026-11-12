"""One bound W with a time cap: is there a characteristic of weight <= W?  (sat / unsat / timeout)

Run: python cap_run.py jlin 16 8 4 13 15 10800        log: results/<mode>_n<n>_r<r>.log
"""
import os, subprocess, sys, time

from gen_model import generate
from check_trail import parse, check
from run_min import wsl_path, sh, STP, HERE


def main():
    mode, n, b1, b2, r, W, cap = sys.argv[1], *map(int, sys.argv[2:8])
    f = os.path.join(HERE, "models", f"{mode}_n{n}_r{r}_w{W}.cvc")
    with open(f, "w") as fh:
        fh.write(generate(mode, n, b1, b2, r, W))
    t = time.time()
    p = subprocess.run(sh(f"timeout {cap} {STP} {wsl_path(f)}"), capture_output=True, text=True)
    dt = time.time() - t
    out = p.stdout
    open(f[:-4] + ".out", "w").write(out + p.stderr)
    head = f"{mode} n={n} ({b1},{b2}) rounds={r} W<={W}:"
    if "Invalid." in out:
        vals = parse(out)
        w = check(mode, vals, n, b1, b2, r)
        line = (f"{head} sat {dt:.1f}s  checked weight={w} total={vals.get('total')}  trail=" +
                " ".join(f"L{i}={vals.get(f'L{i}',0):x},R{i}={vals.get(f'R{i}',0):x}" for i in range(r + 1)))
    elif "Valid." in out:
        line = f"{head} unsat {dt:.1f}s"
    else:
        line = f"{head} timeout/stopped after {dt:.1f}s (cap {cap}s)"
    print(line, flush=True)
    with open(os.path.join(HERE, "results", f"{mode}_n{n}_r{r}.log"), "a") as fh:
        fh.write(line + "\n")


if __name__ == "__main__":
    main()
