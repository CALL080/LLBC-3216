"""One case of impossible.py (to run the cases in parallel).
Run: python imp_one.py jlin 20 2      (mode, rounds, case index 0..3)
"""
import os, sys

import impossible as imp

mode, rounds, i = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
cases = [("in R0 bit0 (part0) -> out L bit4 (part1)", 0, 1, 1 << 4, 0),
         ("in R0 bit0 (part0) -> out L bit8 (part2)", 0, 1, 1 << 8, 0),
         ("in R0 bit0 (part0) -> out R bit12 (part3)", 0, 1, 0, 1 << 12),
         ("in L0 bit63 (part3) -> out R bit0 (part0)", 1 << 63, 0, 0, 1)]
desc, *fx = cases[i]
res, dt = imp.run(mode, rounds, tuple(fx), f"p{i}")
line = f"{mode} rounds={rounds} {desc}: {res}  {dt:.1f}s"
print(line, flush=True)
with open(os.path.join(imp.HERE, "results", f"impossible_r{rounds}.txt"), "a") as fh:
    fh.write(line + "\n")
