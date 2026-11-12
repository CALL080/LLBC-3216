"""Find the minimum characteristic weight for one (mode, n, b1, b2, rounds) with STP.

For W = start, start+1, ...: generate the model, run STP (official static build in WSL);
stop at the first satisfiable W, check the characteristic with check_trail.py, and keep
the model/output of W-1 (unsatisfiable = proof that no lighter characteristic exists).
With --only W: run just that bound (used to confirm a known value: W-1 unsat, W sat).

Run: python run_min.py diff 16 8 4 7 --start 23
Log: results/<mode>_n<n>_r<r>.log (one line per W: result, seconds)
"""
import argparse, os, subprocess, time

from gen_model import generate
from check_trail import parse, check

HERE = os.path.dirname(os.path.abspath(__file__))
# STP binary (static 2.3.3 build); override with the environment variable LLBC_STP.
STP = os.environ.get("LLBC_STP", "~/stp-official/stp")
USE_WSL = os.name == "nt"          # on Windows run Linux tools through WSL, on Linux call them directly


def sh(cmd):
    """argv that runs a shell command line (inside WSL on Windows)."""
    return ["wsl", "-e", "bash", "-lc", cmd] if USE_WSL else ["bash", "-lc", cmd]


def wsl_path(p):
    p = os.path.abspath(p).replace("\\", "/")
    if not USE_WSL:
        return p
    return "/mnt/" + p[0].lower() + p[2:]


def run(mode, n, b1, b2, r, W):
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    f = os.path.join(HERE, "models", f"{mode}_n{n}_r{r}_w{W}.cvc")
    with open(f, "w") as fh:
        fh.write(generate(mode, n, b1, b2, r, W))
    t = time.time()
    p = subprocess.run(sh(f"{STP} {wsl_path(f)}"), capture_output=True, text=True)
    dt = time.time() - t
    out = p.stdout
    with open(f[:-4] + ".out", "w") as fh:
        fh.write(out + p.stderr)
    if "Valid." in out and "Invalid." not in out:
        return "unsat", dt, None
    if "Invalid." in out:
        return "sat", dt, out
    raise RuntimeError(f"unexpected STP output for {f}: {out[-500:]} {p.stderr[-500:]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("n", type=int); ap.add_argument("b1", type=int)
    ap.add_argument("b2", type=int); ap.add_argument("rounds", type=int)
    ap.add_argument("--start", type=int, default=0); ap.add_argument("--only", type=int)
    a = ap.parse_args()
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    log = os.path.join(HERE, "results", f"{a.mode}_n{a.n}_r{a.rounds}.log")
    Ws = [a.only] if a.only is not None else range(a.start, 4095)
    for W in Ws:
        res, dt, out = run(a.mode, a.n, a.b1, a.b2, a.rounds, W)
        line = f"{a.mode} n={a.n} ({a.b1},{a.b2}) rounds={a.rounds} W<={W}: {res} {dt:.1f}s"
        if res == "sat":
            vals = parse(out)
            w = check(a.mode, vals, a.n, a.b1, a.b2, a.rounds)
            line += f"  checked weight={w} total={vals.get('total')}"
            line += "  trail=" + " ".join(f"L{i}={vals.get(f'L{i}',0):x},R{i}={vals.get(f'R{i}',0):x}" for i in range(a.rounds + 1))
        print(line, flush=True)
        with open(log, "a") as fh:
            fh.write(line + "\n")
        if res == "sat":
            break


if __name__ == "__main__":
    main()
