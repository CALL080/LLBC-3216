"""Minimum weight for one round count, trying several bounds W at once (one STP per core).

Bounds start, start+1, ... are run with at most --jobs processes.  When some W is satisfiable and
every smaller W has finished unsatisfiable, W is the minimum; larger runs are stopped.
Each run's result is appended to results/<mode>_n<n>_r<r>.log like run_min.py.

Run: python par_min.py jlin 16 8 4 7 --start 7 --jobs 6
"""
import argparse, os, subprocess, time

from gen_model import generate
from check_trail import parse, check
from run_min import wsl_path, sh, STP, HERE, USE_WSL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("n", type=int); ap.add_argument("b1", type=int)
    ap.add_argument("b2", type=int); ap.add_argument("rounds", type=int)
    ap.add_argument("--start", type=int, required=True); ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    log = os.path.join(HERE, "results", f"{a.mode}_n{a.n}_r{a.rounds}.log")
    running, done, nextW, t0 = {}, {}, a.start, {}

    def launch(W):
        f = os.path.join(HERE, "models", f"{a.mode}_n{a.n}_r{a.rounds}_w{W}.cvc")
        with open(f, "w") as fh:
            fh.write(generate(a.mode, a.n, a.b1, a.b2, a.rounds, W))
        out = open(f[:-4] + ".out", "w")
        running[W] = (subprocess.Popen(sh(f"{STP} {wsl_path(f)}"), stdout=out,
                                       stderr=subprocess.STDOUT, text=True), f, out)
        t0[W] = time.time()

    while True:
        best = min((W for W, r in done.items() if r == "sat"), default=None)
        if best is not None and all(done.get(W) == "unsat" for W in range(a.start, best)):
            for W, (p, f, out) in running.items():
                subprocess.run((["wsl", "-e"] if USE_WSL else []) + ["pkill", "-f", os.path.basename(f)])
                p.wait(); out.close()
            print(f"MIN {a.mode} n={a.n} rounds={a.rounds}: {best}", flush=True)
            return
        while len(running) < a.jobs and (best is None or nextW < best):
            launch(nextW); nextW += 1
        time.sleep(2)
        for W in list(running):
            p, f, out = running[W]
            if p.poll() is None:
                continue
            out.close(); del running[W]
            text = open(f[:-4] + ".out").read()
            dt = time.time() - t0[W]
            if "Invalid." in text:
                vals = parse(text)
                w = check(a.mode, vals, a.n, a.b1, a.b2, a.rounds)
                res = "sat"
                line = (f"{a.mode} n={a.n} ({a.b1},{a.b2}) rounds={a.rounds} W<={W}: sat {dt:.1f}s  checked weight={w} "
                        f"total={vals.get('total')}  trail=" +
                        " ".join(f"L{i}={vals.get(f'L{i}',0):x},R{i}={vals.get(f'R{i}',0):x}" for i in range(a.rounds + 1)))
            elif "Valid." in text:
                res = "unsat"
                line = f"{a.mode} n={a.n} ({a.b1},{a.b2}) rounds={a.rounds} W<={W}: unsat {dt:.1f}s"
            else:
                raise RuntimeError(f"unexpected STP output for W={W}: {text[-300:]}")
            done[W] = res
            print(line, flush=True)
            with open(log, "a") as fh:
                fh.write(line + "\n")


if __name__ == "__main__":
    main()
