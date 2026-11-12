"""Exhaustive output-difference distributions on one 32-bit part (measure.c ddist):
is there any differential (not just a characteristic) usable at 9-11 rounds?

For each round count: the input differences of the best joint characteristics found
(results/jdiff_n16_r<r>.log, up to 3 distinct) and some one-nibble input differences.
Reference: 20 rounds (python ddist.py 20, one-nibble input differences).  For a random permutation
the number of pairs per output difference is Poisson(1/2), so the maximum count over the 2^32
output differences is 18 or 20 ordered x with probability about 0.97 (22 with probability 0.03).
Keys: key schedule of the reference implementation, part 0, master key from random.Random(900 + r).

Run: MEASURE_THREADS=18 python ddist.py 9 10 20      log: results/ddist.txt
"""
import os, re, sys
import measure as m

HERE = os.path.dirname(os.path.abspath(__file__))
NIBBLE = [(0, 0x5), (0, 0x500), (0x5, 0), (0, 0xd)]


def inputs(r):
    seen = []
    f = os.path.join(HERE, "results", f"jdiff_n16_r{r}.log")
    for line in (open(f) if os.path.exists(f) else []):      # no log (e.g. r = 20): one-nibble differences only
        if ": sat" in line:
            v = dict((k, int(x, 16)) for k, x in re.findall(r"([LR]\d+)=([0-9a-f]+)", line))
            d = (v["L0"], v["R0"])
            if d not in seen:
                seen.append(d)
    return seen[:3] + [d for d in NIBBLE if d not in seen]


log = open(os.path.join(HERE, "results", "ddist.txt"), "a")
for r in map(int, sys.argv[1:]):
    K = m.keys(1, 900 + r)[0]
    for dl, dr in inputs(r):
        out, t = m.run(["ddist", f"{dl:x}", f"{dr:x}"], m.keyline(K, r))
        s = f"r={r} din=(L {dl:04x}, R {dr:04x}): " + out.replace("\n", "  ") + f"  [{t}]"
        print(s, flush=True)
        log.write(s + "\n"); log.flush()
