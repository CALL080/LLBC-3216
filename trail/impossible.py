"""Impossible differentials and zero-correlation linear hulls of the full 128-bit official LLBC,
found the usual automated way: fix an input and an output (difference or mask) and ask STP whether
any characteristic connects them (joint S-box model, any weight).  Unsatisfiable = impossible
differential (jdiff) / zero-correlation (jlin, trail-based) for that pair.

Pairs tried: input = one bit in part 0, output = one bit in part 1, 2 or 3 (expected impossible by
Theorem 1), and controls with the output bit in part 0 at a few rounds (expected possible).
Part k = bits i with (i mod 16)//4 == k.

Run: python impossible.py [rounds]       log: results/impossible_r<rounds>.txt
"""
import os, subprocess, sys, time

from gen_model import generate
from run_min import wsl_path, sh, STP, HERE

N, B1, B2 = 64, 32, 16


def run(mode, rounds, fixed, tag):
    f = os.path.join(HERE, "models", f"imp_{mode}_r{rounds}_{tag}.cvc")
    with open(f, "w") as fh:
        fh.write(generate(mode, N, B1, B2, rounds, None, fixed))
    t = time.time()
    p = subprocess.run(sh(f"{STP} {wsl_path(f)}"), capture_output=True, text=True)
    res = "IMPOSSIBLE (unsat)" if "Valid." in p.stdout and "Invalid." not in p.stdout else \
          "possible (sat)" if "Invalid." in p.stdout else "ERROR " + p.stdout[-200:] + p.stderr[-200:]
    return res, time.time() - t


def unpack(w16, k=0):
    """16-bit part word (nibble j = cell k+4j) -> 64-bit word with only part k set."""
    return sum(((w16 >> (4 * j)) & 0xF) << (4 * (k + 4 * j)) for j in range(4))


def control(mode, rounds):
    """Input/output of the best 32-bit characteristic for this round count, placed in part 0
    (must come out possible; shows the model is not vacuous)."""
    import re
    f = os.path.join(HERE, "results", f"{mode}_n16_r{rounds}.log")
    if not os.path.exists(f):
        return None
    sat = [l for l in open(f) if ": sat" in l]
    if not sat:
        return None
    v = dict((k, int(x, 16)) for k, x in re.findall(r"([LR]\d+)=([0-9a-f]+)", sat[-1]))
    return (f"control: best {rounds}-round characteristic placed in part 0",
            unpack(v["L0"]), unpack(v["R0"]), unpack(v[f"L{rounds}"]), unpack(v[f"R{rounds}"]))


def main():
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    log = open(os.path.join(HERE, "results", f"impossible_r{rounds}.txt"), "a")
    # (description, L0, R0, Lr, Rr): input bit 0 of R (part 0); output single bits
    cases = [("in R0 bit0 (part0) -> out L bit4 (part1)", 0, 1, 1 << 4, 0),
             ("in R0 bit0 (part0) -> out L bit8 (part2)", 0, 1, 1 << 8, 0),
             ("in R0 bit0 (part0) -> out R bit12 (part3)", 0, 1, 0, 1 << 12),
             ("in L0 bit63 (part3) -> out R bit0 (part0)", 1 << 63, 0, 0, 1)]
    for mode in ("jdiff", "jlin"):
        mc = list(cases)
        ctl = control(mode, rounds)
        if ctl:
            mc.append(ctl)
        for i, (desc, *fx) in enumerate(mc):
            res, dt = run(mode, rounds, tuple(fx), f"c{i}")
            line = f"{mode} rounds={rounds} {desc}: {res}  {dt:.1f}s"
            print(line, flush=True)
            log.write(line + "\n"); log.flush()


if __name__ == "__main__":
    main()
